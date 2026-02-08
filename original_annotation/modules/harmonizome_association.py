#!/usr/bin/env python3
from __future__ import annotations

"""Pull psychiatric-leaning Harmonizome associations and attach them to gene tables."""

import logging
from http.client import IncompleteRead
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd
import requests
from requests.exceptions import ChunkedEncodingError, ContentDecodingError, ConnectionError

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# API configuration
# ---------------------------------------------------------------------------

HARMONIZOME_API_BASE = "https://maayanlab.cloud/Harmonizome/api/1.0"

# Datasets that represent gene–disease / gene–phenotype associations.
HARMONIZOME_DISEASE_DATASETS = {
    "CTD Gene-Disease Associations",
    "CTD Gene-Disease Associations 2025",
    "DisGeNET Gene-Disease Associations",
    "DisGeNET Gene-Phenotype Associations",
    "GAD Gene-Disease Associations",
    "GAD High Level Gene-Disease Associations",
    "GWASdb SNP-Disease Associations",
    "DISEASES Experimental Gene-Disease Association Evidence Scores",
    "DISEASES Experimental Gene-Disease Association Evidence Scores 2025",
    "DISEASES Text-mining Gene-Disease Association Evidence Scores",
    "DISEASES Text-mining Gene-Disease Association Evidence Scores 2025",
}

# Psychiatric keyword heuristic, aligned roughly with your GWAS psych list.
PSYCH_KEYWORDS: List[str] = [
    # core psychotic/mood
    "schizophrenia",
    "psychosis",
    "psychotic",
    "bipolar",
    "mania",
    "manic",
    "mood disorder",
    "depression",
    "depressive",
    "major depressive",
    "unipolar",

    # anxiety / ocd / ptsd / stress
    "anxiety",
    "panic disorder",
    "agoraphobia",
    "social phobia",
    "phobia",
    "obsessive-compulsive",
    "obsessive compulsive",
    "ocd",
    "post-traumatic stress",
    "posttraumatic stress",
    "ptsd",
    "stress-related",
    "stress related",

    # neurodevelopmental / disruptive
    "autism",
    "asperger",
    "pervasive developmental",
    "adhd",
    "attention deficit",
    "hyperactivity",
    "conduct disorder",
    "oppositional defiant",
    "disruptive behavior",

    # eating / substance
    "eating disorder",
    "anorexia nervosa",
    "bulimia",
    "binge eating",
    "substance use",
    "substance-use",
    "substance abuse",
    "addiction",
    "alcohol dependence",
    "alcohol-use disorder",
    "nicotine dependence",

    # personality & broad psych terms
    "personality disorder",
    "somatoform",
    "tic disorder",
    "tourette",
    "neuropsychiatric",
    "mental disorder",
    "mental or behavioural",
    "mental or behavioral",
    "psychiatric",
    "suicide",
    "suicidal",
    "self-harm",
    "self harm",
    "neuroticism",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _normalize_symbol(symbol: object) -> str:
    """Uppercase, stripped string representation of a gene symbol."""
    if symbol is None or (isinstance(symbol, float) and pd.isna(symbol)):
        return ""
    return str(symbol).strip().upper()


def _is_psych_label(text: object) -> bool:
    """Return True if the label text looks psychiatric based on keywords."""
    if not isinstance(text, str):
        return False
    lt = text.lower()
    return any(kw in lt for kw in PSYCH_KEYWORDS)


# ---------------------------------------------------------------------------
# Harmonizome API interaction
# ---------------------------------------------------------------------------

def fetch_harmonizome_associations_for_gene(
    gene_symbol: str,
    timeout: float = 20.0,
) -> List[Tuple[str, str, Optional[float]]]:
    """
    Query Harmonizome for all associations of a gene and return
    (attribute_name, dataset_name, score) tuples for *disease datasets only*.

    Uses the "download associations" endpoint:

        GET /download/associations?gene=GENE_SYMBOL

    The endpoint returns a plain-text table (tab- or space-separated).
    We are robust to slight format variation; any parsing failures are logged
    and skipped, rather than causing the entire gene to fail.
    """
    gene_symbol = gene_symbol.strip()
    if not gene_symbol:
        return []

    url = f"{HARMONIZOME_API_BASE}/download/associations"
    logger.info("[HARMONIZOME] GET %s?gene=%s", url, gene_symbol)
    text_payload: Optional[str] = None
    resp: Optional[requests.Response] = None

    try:
        resp = requests.get(url, params={"gene": gene_symbol}, timeout=timeout)
        text_payload = resp.text
    except (ChunkedEncodingError, ContentDecodingError, ConnectionError) as exc:
        # Harmonizome sometimes closes chunked responses early (e.g., gene NTM),
        # causing requests to raise before the payload can be parsed. Fall back
        # to urllib and salvage any partial body returned by the server.
        logger.warning(
            "[HARMONIZOME] Chunked download issue for gene %s: %s; retrying with urllib.",
            gene_symbol,
            exc,
        )
        try:
            req = Request(f"{url}?{urlencode({'gene': gene_symbol})}")
            with urlopen(req, timeout=timeout) as fh:  # type: ignore[arg-type]
                try:
                    data = fh.read()
                except IncompleteRead as ir_exc:
                    data = ir_exc.partial
                text_payload = data.decode("utf-8", errors="replace")
        except Exception as exc2:
            logger.warning(
                "[HARMONIZOME] Fallback request failed for gene %s: %s",
                gene_symbol,
                exc2,
            )
            return []
    except Exception as exc:
        logger.warning(
            "[HARMONIZOME] Request failed for gene %s: %s",
            gene_symbol,
            exc,
        )
        return []

    if resp is not None and resp.status_code == 404:
        # No associations for this gene
        logger.debug("[HARMONIZOME] No associations found for gene %s (404).", gene_symbol)
        return []

    if resp is not None and resp.status_code != 200:
        logger.warning(
            "[HARMONIZOME] Unexpected status %s for gene %s",
            resp.status_code,
            gene_symbol,
        )
        return []

    if text_payload is None:
        return []

    lines = text_payload.splitlines()
    results: List[Tuple[str, str, Optional[float]]] = []

    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            # Skip header / comment lines if present
            continue

        parts = line.split("\t")
        if len(parts) < 3:
            # Fallback: try whitespace split if tabs not used
            parts = line.split()
        if len(parts) < 3:
            # Can't reliably parse
            continue

        # Heuristic: try attr, dataset, score (last column)
        attr = parts[0].strip()
        dataset = parts[1].strip()
        score_val = parts[-1].strip()
        try:
            score = float(score_val) if score_val and score_val.lower() != "null" else None
        except ValueError:
            score = None

        if dataset not in HARMONIZOME_DISEASE_DATASETS:
            continue

        results.append((attr, dataset, score))

    return results


# ---------------------------------------------------------------------------
# Core summarisation logic
# ---------------------------------------------------------------------------

def build_harmonizome_summary(
    gene_symbols: Sequence[str],
) -> pd.DataFrame:
    """
    Build a per-gene psychiatric summary table from Harmonizome.

    Parameters
    ----------
    gene_symbols:
        Sequence of gene symbols (HGNC, case-insensitive).

    Returns
    -------
    DataFrame with columns:
        gene_symbol                (uppercased)
        harmonizome_count    (Int64)
        harmonizome_terms    (semicolon-separated labels)
        harmonizome_datasets (semicolon-separated dataset names)
    """
    # Normalise and deduplicate
    symbols_norm = sorted(
        { _normalize_symbol(g) for g in gene_symbols if _normalize_symbol(g) }
    )
    if not symbols_norm:
        return pd.DataFrame(
            columns=[
                "gene_symbol",
                "harmonizome_count",
                "harmonizome_terms",
                "harmonizome_datasets",
            ]
        )

    logger.info("[HARMONIZOME] Building psychiatric summary for %d genes...", len(symbols_norm))

    summary: Dict[str, Dict[str, object]] = {}

    total = len(symbols_norm)
    for idx, sym in enumerate(symbols_norm, start=1):
        logger.info("[HARMONIZOME] (%d/%d) Fetching associations for %s", idx, total, sym)
        assoc = fetch_harmonizome_associations_for_gene(sym)
        if not assoc:
            continue

        psych_terms: set[str] = set()
        psych_datasets: set[str] = set()

        for attr, dataset, score in assoc:
            if _is_psych_label(attr):
                psych_terms.add(attr)
                psych_datasets.add(dataset)

        if psych_terms:
            summary[sym] = {
                "gene_symbol": sym,
                "harmonizome_count": len(psych_terms),
            "harmonizome_terms": "; ".join(sorted(psych_terms)),
            "harmonizome_datasets": "; ".join(sorted(psych_datasets)),
        }

    if not summary:
        return pd.DataFrame(
            columns=[
                "gene_symbol",
                "harmonizome_count",
                "harmonizome_terms",
                "harmonizome_datasets",
            ]
        )

    df = pd.DataFrame.from_records(list(summary.values()))
    df["harmonizome_count"] = pd.to_numeric(
        df["harmonizome_count"], errors="coerce"
    ).astype("Int64")
    logger.info("[HARMONIZOME] Completed summary with %d genes that had psychiatric signals.", len(df))
    return df


def attach_harmonizome(
    gene_db: pd.DataFrame,
    symbol_col: str = "approved_symbol",
) -> pd.DataFrame:
    """
    Attach Harmonizome psychiatric summary columns to a working gene DataFrame.

    Parameters
    ----------
    gene_db:
        Your working gene table (from HGNC + function + GWAS pipeline).
    symbol_col:
        Column in `gene_db` containing HGNC symbols. If not present,
        falls back to 'symbol' if available.

    Returns
    -------
    A copy of `gene_db` with added columns:

        harmonizome_count
        harmonizome_terms
        harmonizome_datasets
    """
    df = gene_db.copy()

    if symbol_col not in df.columns:
        if "symbol" in df.columns:
            symbol_col = "symbol"
        else:
            logger.warning(
                "[HARMONIZOME] No '%s' or 'symbol' column in gene_db; "
                "skipping Harmonizome attachment.",
                symbol_col,
            )
            df["harmonizome_count"] = pd.NA
            df["harmonizome_terms"] = pd.NA
            df["harmonizome_datasets"] = pd.NA
            return df

    logger.info("[HARMONIZOME] Attaching associations using column '%s' (%d genes).", symbol_col, len(df))
    symbols = (
        df[symbol_col]
        .dropna()
        .astype(str)
        .str.strip()
        .tolist()
    )

    if not symbols:
        logger.info("[HARMONIZOME] No valid gene symbols in '%s'; skipping.", symbol_col)
        df["harmonizome_count"] = pd.NA
        df["harmonizome_terms"] = pd.NA
        df["harmonizome_datasets"] = pd.NA
        return df

    hm_df = build_harmonizome_summary(symbols)

    if hm_df.empty:
        logger.info("[HARMONIZOME] No psychiatric associations found for input genes.")
        df["harmonizome_count"] = pd.NA
        df["harmonizome_terms"] = pd.NA
        df["harmonizome_datasets"] = pd.NA
        return df

    # Prepare merge keys
    hm_df["gene_symbol_upper"] = hm_df["gene_symbol"].apply(_normalize_symbol)
    df["__symbol_upper"] = df[symbol_col].apply(_normalize_symbol)

    df = df.merge(
        hm_df[["gene_symbol_upper", "harmonizome_count",
               "harmonizome_terms", "harmonizome_datasets"]],
        how="left",
        left_on="__symbol_upper",
        right_on="gene_symbol_upper",
    )

    df = df.drop(columns=["__symbol_upper", "gene_symbol_upper"], errors="ignore")

    # Ensure Int64 dtype for count column
    if "harmonizome_count" in df.columns:
        df["harmonizome_count"] = pd.to_numeric(
            df["harmonizome_count"], errors="coerce"
        ).astype("Int64")

    attached = df["harmonizome_count"].notna().sum()
    logger.info("[HARMONIZOME] Attached psychiatric signals for %d/%d genes.", attached, len(df))

    return df


# ---------------------------------------------------------------------------
# Minimal CLI for testing
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    import sys

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    parser = argparse.ArgumentParser(
        description="Query Harmonizome for psychiatric-like disease associations."
    )
    parser.add_argument(
        "genes",
        nargs="*",
        help="Gene symbols (HGNC), e.g. NTM COMT DRD2. "
             "If omitted, use --gene-file or --tsv instead.",
    )
    parser.add_argument(
        "--gene-file",
        type=str,
        help="Optional text file with one gene symbol per line.",
    )
    parser.add_argument(
        "--tsv",
        type=str,
        help="Optional TSV with a gene column; prints TSV with Harmonizome columns attached.",
    )
    parser.add_argument(
        "--gene-col",
        type=str,
        default="approved_symbol",
        help="Gene symbol column name in TSV (default: approved_symbol).",
    )

    args = parser.parse_args()

    # Case 1: attach to an existing TSV
    if args.tsv:
        tsv_path = Path(args.tsv)
        if not tsv_path.is_file():
            logger.error("Input TSV not found: %s", tsv_path)
            sys.exit(1)
        df_in = pd.read_csv(tsv_path, sep="\t", dtype=str)
        df_out = attach_harmonizome(df_in, symbol_col=args.gene_col)
        out_path = tsv_path.with_suffix(".harmonizome.tsv")
        df_out.to_csv(out_path, sep="\t", index=False)
        logger.info("Wrote TSV with Harmonizome columns to %s", out_path)
        sys.exit(0)

    # Case 2: simple gene list and print summary
    genes: List[str] = []
    if args.gene_file:
        path = Path(args.gene_file)
        if not path.is_file():
            logger.error("Gene file not found: %s", path)
            sys.exit(1)
        genes = (
            pd.read_csv(path, header=None, names=["gene"], dtype=str)
            ["gene"].dropna().astype(str).str.strip().tolist()
        )

    if args.genes:
        genes.extend(args.genes)

    if not genes:
        logger.error("No genes provided. Use positional genes, --gene-file, or --tsv.")
        sys.exit(1)

    summary_df = build_harmonizome_summary(genes)
    if summary_df.empty:
        print("No psychiatric-like Harmonizome associations found for input genes.")
    else:
        print(summary_df.to_string(index=False))




