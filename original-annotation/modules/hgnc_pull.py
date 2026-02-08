"""
Fetch core gene metadata from NCBI Datasets and fill gaps with HGNC where
possible. Handles symbol or LOC-style inputs, records simple logs, and returns
one row per requested gene with IDs and summaries needed downstream.
"""

from __future__ import annotations

import json
import time
from typing import Dict, List, Optional

import pandas as pd
import requests


NCBI_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2/gene"
HGNC_BASE = "https://rest.genenames.org"


def _http_get_with_retries(
    url: str,
    headers: Optional[Dict[str, str]] = None,
    max_retries: int = 3,
    delay: float = 1.0,
) -> Optional[dict]:
    """Simple GET with basic retry logic and print-based logging."""
    label = "HTTP"
    if "ncbi.nlm.nih.gov" in url:
        label = "NCBI"
    elif "genenames.org" in url:
        label = "HGNC"

    for attempt in range(1, max_retries + 1):
        try:
            print(f"[{label}] GET {url} (attempt {attempt})")
            resp = requests.get(url, headers=headers or {}, timeout=10)
            if resp.status_code == 200:
                return resp.json()
            else:
                print(f"[{label}] Non-200 status {resp.status_code}: {resp.text[:200]}")
        except requests.RequestException as e:
            print(f"[{label}] Request error on attempt {attempt}: {e}")
        time.sleep(delay)
    print(f"[{label}] Failed after {max_retries} attempts for URL: {url}")
    return None


def _parse_ncbi_gene_report(raw_input: str, data: dict) -> Dict[str, object]:
    """Extracts relevant fields from NCBI Datasets gene JSON for a single raw_input."""
    row: Dict[str, object] = {
        "raw_input": raw_input,
        "approved_symbol": None,
        "hgnc_id": None,
        "entrez_id": None,
        "ensembl_id": None,
        "gene_type_raw": None,
        "synonyms": None,
        "description": None,
        "summary_text": None,
        # "gene_ontology": None,
        "uniprot_primary_accession": None,
        "annotation_status": "unknown",
        "log": "",
    }

    try:
        reports = data.get("reports") or []
        if not reports:
            row["annotation_status"] = "ncbi_not_found"
            row["log"] = "No reports returned by NCBI"
            return row

        gene = reports[0].get("gene", {})
        # Keep as numeric (later coerced to pandas Int64 for NA support)
        row["entrez_id"] = gene.get("gene_id")
        row["approved_symbol"] = gene.get("symbol")
        row["description"] = gene.get("description")

        # HGNC
        na = gene.get("nomenclature_authority") or {}
        row["hgnc_id"] = na.get("identifier")

        # IDs
        ensembl_ids = gene.get("ensembl_gene_ids") or []
        row["ensembl_id"] = ensembl_ids[0] if ensembl_ids else None

        # gene type
        gene_type = gene.get("type")
        row["gene_type_raw"] = gene_type

        # synonyms (as JSON list string)
        synonyms = gene.get("synonyms") or []
        row["synonyms"] = json.dumps(synonyms)

        # UniProt primary accession (Swiss-Prot)
        swiss = gene.get("swiss_prot_accessions") or []
        row["uniprot_primary_accession"] = swiss[0] if swiss else None

        # summary text
        summaries = gene.get("summary") or []
        row["summary_text"] = summaries[0].get("description") if summaries else None

        # gene ontology (store whole dict as JSON string)
        # go = gene.get("gene_ontology")
        # row["gene_ontology"] = json.dumps(go) if go is not None else None

        row["annotation_status"] = "ok"
    except Exception as e:  # be defensive; never crash the pipeline for one gene
        row["annotation_status"] = "parse_error"
        row["log"] = f"Error parsing NCBI response: {e}"

    return row


def _append_log(existing: object, message: str) -> str:
    """Append a log message, preserving any existing text."""
    base = str(existing).strip() if isinstance(existing, str) else ""
    return f"{base}; {message}" if base else message


def _fetch_ensembl_from_hgnc(hgnc_id: str) -> Optional[str]:
    """
    Fallback lookup for Ensembl ID via HGNC REST API.
    Returns string Ensembl ID or None if not present.
    """
    if not hgnc_id:
        return None

    url = f"{HGNC_BASE}/fetch/hgnc_id/{hgnc_id}"
    headers = {"Accept": "application/json"}
    data = _http_get_with_retries(url, headers=headers)
    if not data:
        return None

    response = data.get("response") or {}
    docs = response.get("docs") or []
    if not docs:
        return None

    doc = docs[0]
    ensembl_id = doc.get("ensembl_gene_id") or doc.get("ensembl_id")

    if isinstance(ensembl_id, list):
        return ensembl_id[0] if ensembl_id else None

    return ensembl_id


def _build_ncbi_url_for_gene(raw_input: str) -> str:
    """
    Decide which NCBI endpoint to use:
    - If raw_input starts with LOC + digits: use gene ID endpoint with numeric part.
    - Else: use symbol + human taxon endpoint.
    """
    s = raw_input.strip()
    if s.upper().startswith("LOC") and s[3:].isdigit():
        entrez_id = s[3:]
        return f"{NCBI_BASE}/id/{entrez_id}"
    # default: treat as symbol, restrict to human
    symbol = s
    return f"{NCBI_BASE}/symbol/{symbol}/taxon/9606"


def run_hgnc_pull(
    gene_list_path: str,
    existing_db: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """
    NCBI-based gene annotation module.
    Reads a list of gene identifiers and annotates them using NCBI Datasets API.
    Keeps HGNC ID from NCBI when available.
    """

    # Read gene list
    with open(gene_list_path, "r") as f:
        raw_genes = [line.strip() for line in f if line.strip() and not line.strip().startswith("#")]

    # Deduplicate while preserving order
    seen = set()
    genes: List[str] = []
    for g in raw_genes:
        if g not in seen:
            seen.add(g)
            genes.append(g)

    print(f"[NCBI] Total unique genes to process: {len(genes)}")

    rows: List[Dict[str, object]] = []

    for raw_input in genes:
        url = _build_ncbi_url_for_gene(raw_input)
        data = _http_get_with_retries(url)
        if data is None:
            rows.append(
                {
                    "raw_input": raw_input,
                    "approved_symbol": None,
                    "hgnc_id": None,
                    "entrez_id": None,
                    "ensembl_id": None,
                    "gene_type_raw": None,
                    "synonyms": None,
                    "description": None,
                    "summary_text": None,
                    # "gene_ontology": None,
                    "uniprot_primary_accession": None,
                    "annotation_status": "network_error",
                    "log": "NCBI request failed after retries",
                }
            )
            continue

        row = _parse_ncbi_gene_report(raw_input, data)

        # Fallback: if NCBI did not provide an Ensembl ID, try HGNC
        if row.get("annotation_status") == "ok" and (not row.get("ensembl_id") or pd.isna(row.get("ensembl_id"))):
            if row.get("hgnc_id"):
                fallback_ensembl = _fetch_ensembl_from_hgnc(str(row["hgnc_id"]))
                if fallback_ensembl:
                    row["ensembl_id"] = fallback_ensembl
                else:
                    row["log"] = _append_log(row.get("log", ""), f"Missing Ensembl ID from NCBI; HGNC lookup for {row['hgnc_id']} also missing")
            else:
                row["log"] = _append_log(row.get("log", ""), "Missing Ensembl ID from NCBI and no HGNC ID available for fallback")

        rows.append(row)

    new_df = pd.DataFrame(rows)
    if "entrez_id" in new_df.columns:
        # Ensure Entrez IDs are numeric for cleaner downstream CSV/Excel
        new_df["entrez_id"] = pd.to_numeric(new_df["entrez_id"], errors="coerce").astype("Int64")

    if existing_db is not None and not existing_db.empty:
        # Drop deprecated columns if they exist in an incoming database
        existing_db = existing_db.drop(columns=["gene_type_snake", "locus_type_simplified"], errors="ignore")

        # Merge on raw_input, with new annotations overwriting old
        existing_db = existing_db.set_index("raw_input")
        new_df = new_df.set_index("raw_input")
        combined = existing_db.combine_first(new_df)
        # Overwrite with any new values where present
        combined.update(new_df)
        combined = combined.reset_index()
        return combined

    return new_df
