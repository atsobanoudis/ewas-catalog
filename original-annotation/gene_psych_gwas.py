#!/usr/bin/env python3
"""Summarize psychiatric GWAS Catalog associations per gene into `gene_psych_gwas.tsv`."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import List, Callable

import pandas as pd

# ---------------------------------------------------------------------
RESOURCES_DIR = Path("resources")
GWAS_RESOURCES_DIR = RESOURCES_DIR / "gwas"

GWAS_ASSOC_TSV = GWAS_RESOURCES_DIR / "gwas-catalog-associations_ontology-annotated.tsv"

GENE_PSYCH_OUT = GWAS_RESOURCES_DIR / "gene_psych_gwas.tsv"

# ---------------------------------------------------------------------
# Psychiatric keyword heuristics
# ---------------------------------------------------------------------

PSYCH_KEYWORDS: List[str] = [
    # core psych disorders
    "schizophrenia",
    "psychosis",
    "psychotic",
    "bipolar",
    "mania",
    "mood disorder",
    "depression",
    "depressive",
    "major depressive",
    "unipolar",
    "dysthymia",

    # neurodevelopmental / autism / adhd
    "autism",
    "asperger",
    "pervasive developmental",
    "adhd",
    "attention deficit",
    "hyperactivity",
    "conduct disorder",
    "oppositional defiant",
    "disruptive behavior",

    # anxiety / ocd / ptsd / stress
    "anxiety",
    "panic disorder",
    "agoraphobia",
    "social phobia",
    "obsessive-compulsive",
    "obsessive compulsive",
    "ocd",
    "post-traumatic stress",
    "posttraumatic stress",
    "ptsd",
    "stress-related",
    "stress related",

    # eating / substance
    "eating disorder",
    "anorexia nervosa",
    "bulimia",
    "binge eating",
    "substance use",
    "substance-use",
    "substance abuse",
    "drug dependence",
    "alcohol dependence",
    "alcohol-use disorder",
    "nicotine dependence",

    # other psych/neuropsychiatric baskets
    "personality disorder",
    "somatoform",
    "tic disorder",
    "tourette",
    "neuropsychiatric",
    "mental disorder",
    "mental or behavioural",
    "mental or behavioral",
]

# Column names in gwas-catalog-associations_ontology-annotated.tsv
# (per GWAS Catalog / Dipper docs)
TRAIT_COL = "MAPPED_TRAIT"
DISEASE_COL = "DISEASE/TRAIT"
EFO_URI_COL = "MAPPED_TRAIT_URI"
GENE_COL = "MAPPED_GENE"
PMID_COL = "PUBMEDID"
STUDY_ACCESSION_COL = "STUDY ACCESSION"


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def _ensure_dirs() -> None:
    """Create resources/gwas if needed."""
    for d in (RESOURCES_DIR, GWAS_RESOURCES_DIR):
        d.mkdir(parents=True, exist_ok=True)


def load_associations(path: Path = GWAS_ASSOC_TSV) -> pd.DataFrame:
    """Load GWAS associations table."""
    if not path.is_file():
        raise FileNotFoundError(
            f"GWAS associations file not found at {path}.\n"
            "Download it from the GWAS Catalog FTP 'latest' release:\n"
            "  ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
            "gwas-catalog-associations_ontology-annotated.tsv"
        )
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    return df


def is_psych_trait_text(text: str) -> bool:
    """Heuristic: does this free-text trait description look psychiatric?"""
    if not isinstance(text, str):
        return False
    lt = text.lower()
    return any(kw in lt for kw in PSYCH_KEYWORDS)


def _make_agg_unique() -> Callable[[pd.Series], str]:
    """Return a function that aggregates a Series into sorted, deduped '; ' string."""
    def _fn(series: pd.Series) -> str:
        vals = {
            str(x).strip()
            for x in series
            if pd.notna(x) and str(x).strip()
        }
        return "; ".join(sorted(vals)) if vals else ""
    return _fn


# ---------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------

def build_gene_psych_table(
    assoc_path: Path = GWAS_ASSOC_TSV,
) -> pd.DataFrame:
    """
    Build gene-level psychiatric GWAS association summary by filtering traits,
    exploding multi-gene rows, and aggregating per symbol.
    """
    df = load_associations(assoc_path)

    missing_cols = [c for c in [TRAIT_COL, DISEASE_COL, EFO_URI_COL,
                                GENE_COL, PMID_COL, STUDY_ACCESSION_COL]
                    if c not in df.columns]
    if missing_cols:
        raise ValueError(
            f"Expected columns not found in associations file: {missing_cols}\n"
            "Check that you downloaded the correct 'gwas-catalog-associations_ontology-annotated.tsv' "
            "and that the header has not changed."
        )

    logging.info("Associations loaded: %d rows", len(df))

    # Combine trait text fields
    logging.info("Flagging psychiatric traits using keyword heuristics...")
    def _row_trait_text(row) -> str:
        parts = []
        for col in (TRAIT_COL, DISEASE_COL):
            val = row.get(col)
            if pd.notna(val):
                parts.append(str(val))
        return " ".join(parts)

    trait_text_series = df.apply(_row_trait_text, axis=1)
    psych_mask = trait_text_series.apply(is_psych_trait_text)
    psych_df = df[psych_mask].copy()

    logging.info(
        "Psychiatric-like associations: %d (%.2f%%)",
        len(psych_df),
        100.0 * len(psych_df) / max(len(df), 1),
    )

    # Handle genes
    # MAPPED_GENE is often comma-separated list of HGNC symbols
    psych_df[GENE_COL] = psych_df[GENE_COL].fillna("")
    psych_df = psych_df[psych_df[GENE_COL].str.strip() != ""].copy()

    if psych_df.empty:
        logging.warning("No psychiatric associations with non-empty MAPPED_GENE found.")
        return pd.DataFrame(
            columns=[
                "gene",
                "n_psych_associations",
                "psych_mapped_traits",
                "disease_trait_labels",
                "efo_uris",
                "pubmed_ids",
                "study_accessions",
            ]
        )

    logging.info("Psych associations with mapped genes: %d rows", len(psych_df))

    psych_df["__gene_list"] = psych_df[GENE_COL].str.split(r"\s*,\s*")

    exploded = psych_df.explode("__gene_list").rename(columns={"__gene_list": "gene"})
    exploded["gene"] = exploded["gene"].str.strip()
    exploded = exploded[exploded["gene"] != ""].copy()

    logging.info("After exploding multi-gene rows: %d gene-association rows", len(exploded))

    agg_unique = _make_agg_unique()

    grouped = (
        exploded.groupby("gene", as_index=False)
        .agg(
            n_psych_associations=("gene", "size"),
            psych_mapped_traits=(TRAIT_COL, agg_unique),
            disease_trait_labels=(DISEASE_COL, agg_unique),
            efo_uris=(EFO_URI_COL, agg_unique),
            pubmed_ids=(PMID_COL, agg_unique),
            study_accessions=(STUDY_ACCESSION_COL, agg_unique),
        )
    )

    logging.info("Gene-level psych table: %d genes", len(grouped))

    return grouped


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build a gene-level table of psychiatric GWAS associations "
            "from gwas-catalog-associations_ontology-annotated.tsv"
        )
    )
    parser.add_argument(
        "--assoc-path",
        type=str,
        default=str(GWAS_ASSOC_TSV),
        help=f"Path to gwas-catalog-associations_ontology-annotated.tsv "
             f"(default: {GWAS_ASSOC_TSV})",
    )
    parser.add_argument(
        "--out",
        type=str,
        default=str(GENE_PSYCH_OUT),
        help=f"Output TSV path (default: {GENE_PSYCH_OUT})",
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO).",
    )
    parser.add_argument(
        "--gene",
        type=str,
        default=None,
        help="Optional: after building the full table, print the row for this gene symbol.",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(levelname)s: %(message)s",
    )

    _ensure_dirs()

    assoc_path = Path(args.assoc_path)
    out_path = Path(args.out)

    logging.info("Building gene-level psychiatric GWAS table from: %s", assoc_path)
    gene_df = build_gene_psych_table(assoc_path=assoc_path)

    if gene_df.empty:
        logging.warning(
            "No psychiatric gene associations found; NOT writing an empty file."
        )
    else:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        gene_df.to_csv(out_path, sep="\t", index=False)
        logging.info("Wrote gene-level psych GWAS table to: %s", out_path)

    if args.gene:
        gene = args.gene.strip()
        if gene:
            hit = gene_df[gene_df["gene"].str.upper() == gene.upper()]
            if hit.empty:
                logging.info("No psychiatric GWAS associations found for gene: %s", gene)
            else:
                # Pretty-print to stdout
                print("\n==== Psychiatric GWAS associations for gene:", gene, "====")
                print(hit.to_string(index=False))


if __name__ == "__main__":
    main()
