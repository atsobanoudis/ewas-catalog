#!/usr/bin/env python3
from __future__ import annotations

"""Attach precomputed psychiatric GWAS associations to a gene DataFrame."""

import logging
from pathlib import Path
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)

# Default path to psychiatric GWAS summary
DEFAULT_PSYCH_GWAS_PATH = Path("resources") / "gwas" / "gene_psych_gwas.tsv"

# Name of the log column in the working DataFrame
DEFAULT_LOG_COL = "log"

# Canonical GWAS column names and their accepted legacy aliases
COLUMN_ALIASES = {
    "gwas_assoc_count": ["gwas_assoc_count", "n_psych_associations"],
    "gwas_traits": ["gwas_traits", "psych_mapped_traits"],
    "gwas_labels": ["gwas_labels", "disease_trait_labels"],
    "gwas_pmids": ["gwas_pmids", "pubmed_ids"],
    "gwas_efo_uris": ["gwas_efo_uris", "efo_uris"],
    "gwas_study_accessions": ["gwas_study_accessions", "study_accessions"],
}

REQUIRED_CANONICAL = (
    "gwas_assoc_count",
    "gwas_traits",
    "gwas_labels",
    "gwas_pmids",
)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _ensure_log_column(df: pd.DataFrame, log_col: str = DEFAULT_LOG_COL) -> pd.DataFrame:
    """Ensure the working DataFrame has a string log column."""
    if log_col not in df.columns:
        df[log_col] = ""
    else:
        # Force to string (but keep empty where NA)
        df[log_col] = df[log_col].fillna("").astype(str)
    return df


def _add_empty_gwas_columns(df: pd.DataFrame) -> None:
    """Add empty GWAS columns with consistent dtypes."""
    df["gwas_assoc_count"] = pd.Series([pd.NA] * len(df), dtype="Int64")
    df["gwas_traits"] = ""
    df["gwas_labels"] = ""
    df["gwas_pmids"] = ""
    df["gwas_efo_uris"] = ""
    df["gwas_study_accessions"] = ""
    df["gwas_traits_example"] = ""


def _append_log(df: pd.DataFrame, message: str, mask: Optional[pd.Series] = None,
                log_col: str = DEFAULT_LOG_COL) -> None:
    """
    Append a message to the log column.

    If mask is None, apply the message to all rows; otherwise, only rows
    where mask is True get the message.
    """
    if mask is None:
        mask = pd.Series([True] * len(df), index=df.index)

    if df.empty or not mask.any():
        return

    # Only update rows where mask is True
    current = df.loc[mask, log_col].fillna("").astype(str)

    def _merge_msg(old: str) -> str:
        if not old:
            return message
        # Separate multiple messages with " | "
        return f"{old} | {message}"

    df.loc[mask, log_col] = current.apply(_merge_msg)


def _normalize_symbol(series: pd.Series) -> pd.Series:
    """Uppercase and strip whitespace for gene symbols."""
    return (
        series.fillna("")
        .astype(str)
        .str.strip()
        .str.upper()
    )


def _resolve_canonical_columns(psych_df: pd.DataFrame) -> tuple[dict, list[str]]:
    """
    Map canonical GWAS column names to the first available alias in the GWAS table.

    Returns a tuple: (resolved_map, missing_required)
      - resolved_map: dict of canonical_name -> column_found_in_psych_df
      - missing_required: list of canonical names that were not found
    """
    resolved: dict = {}
    for canonical, aliases in COLUMN_ALIASES.items():
        for col in aliases:
            if col in psych_df.columns:
                resolved[canonical] = col
                break

    missing_required = [col for col in REQUIRED_CANONICAL if col not in resolved]
    return resolved, missing_required


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def attach_psychiatric_gwas(
    df: pd.DataFrame,
    gene_col: str = "gene_symbol",
    psych_path: Optional[Path] = None,
    log_col: str = DEFAULT_LOG_COL,
) -> pd.DataFrame:
    """
    Attach psychiatric GWAS associations to a working gene DataFrame.

    Parameters
    ----------
    df:
        Working DataFrame. Must have one row per gene (or at least a gene
        symbol column).
    gene_col:
        Name of the column in `df` containing HGNC gene symbols.
        Example: "approved_symbol" or "gene_symbol".
    psych_path:
        Path to resources/gwas/gene_psych_gwas.tsv.
        If None, DEFAULT_PSYCH_GWAS_PATH is used.
    log_col:
        Name of the log column in `df` to store textual error/info messages.

    Returns
    -------
    A new DataFrame with additional columns:

        gwas_assoc_count     (Int64; count of psych GWAS hits)
        gwas_traits          (string; GWAS trait names)
        gwas_labels          (string; DISEASE/TRAIT labels)
        gwas_pmids           (string; PMIDs, semicolon-separated)
        gwas_efo_uris        (string; optional EFO URIs)
        gwas_study_accessions (string; optional study accessions)
        gwas_traits_example  (string; first trait for quick inspection)

    Any major problems are recorded in the `log` column rather than raising.
    """
    df = df.copy()
    df = _ensure_log_column(df, log_col=log_col)

    logger.info("[GWAS] Starting attachment using gene column '%s'.", gene_col)

    # ------------------------------------------------------------------
    # 1. Sanity-check gene column
    # ------------------------------------------------------------------
    if gene_col not in df.columns:
        msg = (
            f"[GWAS] Gene column '{gene_col}' not found in DataFrame; "
            "cannot attach psychiatric GWAS associations."
        )
        logger.error(msg)
        _append_log(df, msg, log_col=log_col)
        _add_empty_gwas_columns(df)
        return df

    # ------------------------------------------------------------------
    # 2. Load psychiatric GWAS summary table
    # ------------------------------------------------------------------
    if psych_path is None:
        psych_path = DEFAULT_PSYCH_GWAS_PATH
    else:
        psych_path = Path(psych_path)

    logger.info("[GWAS] Reading psychiatric GWAS summary from %s", psych_path)
    if not psych_path.is_file():
        msg = (
            f"[GWAS] Psychiatric GWAS summary file not found at '{psych_path}'. "
            "Run gene_psych_gwas.py first to generate it."
        )
        logger.error(msg)
        _append_log(df, msg, log_col=log_col)
        # Add empty columns and return
        _add_empty_gwas_columns(df)
        return df

    try:
        psych_df = pd.read_csv(psych_path, sep="\t", dtype=str)
    except Exception as exc:
        msg = f"[GWAS] Failed to read psychiatric GWAS file '{psych_path}': {exc}"
        logger.error(msg)
        _append_log(df, msg, log_col=log_col)
        _add_empty_gwas_columns(df)
        return df
    else:
        logger.info("[GWAS] Loaded %d GWAS records.", len(psych_df))

    if "gene" not in psych_df.columns:
        msg = (
            f"[GWAS] Psychiatric GWAS file '{psych_path}' is missing the 'gene' column; "
            "cannot attach associations."
        )
        logger.error(msg)
        _append_log(df, msg, log_col=log_col)
        _add_empty_gwas_columns(df)
        return df

    resolved_cols, missing_required = _resolve_canonical_columns(psych_df)
    if missing_required:
        msg = (
            f"[GWAS] Psychiatric GWAS file '{psych_path}' is missing expected "
            f"columns: {sorted(missing_required)}; cannot attach associations."
        )
        logger.error(msg)
        _append_log(df, msg, log_col=log_col)
        _add_empty_gwas_columns(df)
        return df

    # ------------------------------------------------------------------
    # 3. Normalize join keys (HGNC symbols)
    # ------------------------------------------------------------------
    df["__gwas_gene_key"] = _normalize_symbol(df[gene_col])
    psych_df["__gwas_gene_key"] = _normalize_symbol(psych_df["gene"])

    # Rows in df with missing or blank gene symbols
    missing_gene_mask = df["__gwas_gene_key"] == ""
    if missing_gene_mask.any():
        msg = (
            "[GWAS] Some rows have empty or missing gene symbols in "
            f"'{gene_col}'; GWAS associations cannot be attached for these."
        )
        logger.warning(msg)
        _append_log(df, msg, mask=missing_gene_mask, log_col=log_col)

    # ------------------------------------------------------------------
    # 4. Merge
    # ------------------------------------------------------------------
    rename_map = {src: canonical for canonical, src in resolved_cols.items()}
    psych_df = psych_df.rename(columns=rename_map)

    merge_cols = ["__gwas_gene_key"] + list(rename_map.values())
    psych_subset = psych_df[[c for c in merge_cols if c in psych_df.columns]]

    merged = df.merge(
        psych_subset,
        on="__gwas_gene_key",
        how="left",
        suffixes=("", "_gwas"),
    )

    # Ensure optional columns exist even if absent upstream
    for col in COLUMN_ALIASES:
        if col not in merged.columns:
            if col == "gwas_assoc_count":
                merged[col] = pd.Series([pd.NA] * len(merged), dtype="Int64")
            else:
                merged[col] = ""

    # ------------------------------------------------------------------
    # 5. Create final GWAS columns
    # ------------------------------------------------------------------
    # Count of associations: Int64 so we can represent missing as <NA>
    assoc_series = merged["gwas_assoc_count"]
    merged["gwas_assoc_count"] = (
        pd.to_numeric(assoc_series, errors="coerce").astype("Int64")
    )

    # Textual fields: fill NaN with "" for clarity
    text_cols = [
        "gwas_traits",
        "gwas_labels",
        "gwas_pmids",
        "gwas_efo_uris",
        "gwas_study_accessions",
    ]
    for col in text_cols:
        merged[col] = merged[col].fillna("").astype(str)

    merged["gwas_traits_example"] = (
        merged["gwas_traits"]
        .fillna("")
        .astype(str)
        .apply(lambda s: s.split(";")[0].strip() if s else "")
    )

    # ------------------------------------------------------------------
    # 6. Clean up temp columns and return
    # ------------------------------------------------------------------
    merged = merged.drop(columns=["__gwas_gene_key"], errors="ignore")
    legacy_cols = [
        "n_psych_associations",
        "psych_mapped_traits",
        "disease_trait_labels",
        "pubmed_ids",
        "gwas_associations",
        "gwas_trait",
        "gwas_label",
    ]
    merged = merged.drop(columns=[c for c in legacy_cols if c in merged.columns], errors="ignore")

    hit_count = merged["gwas_assoc_count"].notna().sum()
    logger.info("[GWAS] Attached psychiatric GWAS data for %d/%d genes.", hit_count, len(merged))

    return merged


# ---------------------------------------------------------------------------
# Minimal CLI for testing
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    import sys

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    parser = argparse.ArgumentParser(
        description="Attach psychiatric GWAS associations to a gene table TSV."
    )
    parser.add_argument(
        "input_tsv",
        help="Path to input TSV with at least a gene column.",
    )
    parser.add_argument(
        "output_tsv",
        help="Path to output TSV with GWAS columns added.",
    )
    parser.add_argument(
        "--gene-col",
        default="gene_symbol",
        help="Name of the gene symbol column in the input TSV (default: gene_symbol).",
    )
    parser.add_argument(
        "--psych-gwas",
        default=str(DEFAULT_PSYCH_GWAS_PATH),
        help="Path to resources/gwas/gene_psych_gwas.tsv (default: that path).",
    )
    parser.add_argument(
        "--log-col",
        default=DEFAULT_LOG_COL,
        help="Name of the log column in the output (default: 'log').",
    )

    args = parser.parse_args()

    try:
        df_in = pd.read_csv(args.input_tsv, sep="\t", dtype=str)
    except Exception as exc:
        logger.error("Failed to read input TSV '%s': %s", args.input_tsv, exc)
        sys.exit(1)

    df_out = attach_psychiatric_gwas(
        df_in,
        gene_col=args.gene_col,
        psych_path=Path(args.psych_gwas),
        log_col=args.log_col,
    )

    df_out.to_csv(args.output_tsv, sep="\t", index=False)
    logger.info("Wrote output with GWAS columns to %s", args.output_tsv)
