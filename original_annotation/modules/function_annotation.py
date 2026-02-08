"""
Add function annotations to a gene table using NCBI summaries plus UniProt
and HGNC lookups. Produces per-source fields and a combined function block
that labels duplicates across sources.
"""

from __future__ import annotations

import time
from typing import Dict, Optional

import pandas as pd
import requests


UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
HGNC_BASE = "https://rest.genenames.org"  # Accept: application/json
NA_PLACEHOLDER = "n/a"


def _http_get_json_with_retries(
    url: str,
    headers: Optional[Dict[str, str]] = None,
    max_retries: int = 3,
    delay: float = 1.0,
) -> Optional[dict]:
    """Generic JSON GET with simple retry logic."""
    for attempt in range(1, max_retries + 1):
        try:
            print(f"[HTTP] GET {url} (attempt {attempt})")
            resp = requests.get(url, headers=headers or {}, timeout=15)
            if resp.status_code == 200:
                try:
                    return resp.json()
                except ValueError:
                    print(f"[HTTP] Failed to decode JSON for {url}")
                    return None
            else:
                print(f"[HTTP] Non-200 status {resp.status_code}: {resp.text[:200]}")
                if resp.status_code in (400, 404):
                    # Invalid ID; don't bother retrying
                    return None
        except requests.RequestException as e:
            print(f"[HTTP] Request error on attempt {attempt}: {e}")
        time.sleep(delay)
    print(f"[HTTP] Failed after {max_retries} attempts for URL: {url}")
    return None


# ------------------------ UniProt helpers ------------------------------------


def _parse_uniprot_info(entry: dict) -> Dict[str, Optional[str]]:
    """
    Extract function, entry_status, and protein_existence from a UniProt entry.

    Pulls FUNCTION comments, entryType, and proteinExistence fields.
    """
    info: Dict[str, Optional[str]] = {
        "function": None,
        "entry_status": None,
        "protein_existence": None,
    }

    if not entry:
        return info

    # Entry status (reviewed vs unreviewed)
    entry_type = entry.get("entryType")
    if entry_type == "UniProtKB reviewed (Swiss-Prot)":
        info["entry_status"] = "Reviewed"
    elif entry_type == "UniProtKB unreviewed (TrEMBL)":
        info["entry_status"] = "Unreviewed"
    else:
        info["entry_status"] = entry_type

    # Protein existence (e.g., 'Evidence at protein level', etc.)
    pe = entry.get("proteinExistence")
    if isinstance(pe, str):
        info["protein_existence"] = pe
    elif pe is not None:
        info["protein_existence"] = str(pe)

    # Function text from comments
    # NOTE: UniProt JSON uses 'commentType', not 'type'
    comments = entry.get("comments") or []
    function_texts = []
    for comment in comments:
        if comment.get("commentType") == "FUNCTION":
            for text_obj in comment.get("texts", []):
                val = text_obj.get("value")
                if val:
                    function_texts.append(val.strip())

    if function_texts:
        info["function"] = " ".join(function_texts)

    return info


def fetch_uniprot_info(accession: str) -> Dict[str, Optional[str]]:
    """
    Fetch UniProt info for a given accession:
      - function
      - entry_status
      - protein_existence

    Uses /uniprotkb/{accession}.json endpoint.
    """
    empty = {"function": None, "entry_status": None, "protein_existence": None}

    if not accession:
        return empty

    url = f"{UNIPROT_BASE}/{accession}.json"
    data = _http_get_json_with_retries(url)
    if not data:
        return empty

    return _parse_uniprot_info(data)


# ------------------------ HGNC helpers ---------------------------------------


def _parse_hgnc_function(json_data: dict) -> Optional[str]:
    """
    Extract a function-like text field from HGNC JSON if available.
    HGNC rarely has a narrative function; this is best-effort.
    """
    if not json_data:
        return None

    response = json_data.get("response") or {}
    docs = response.get("docs") or []
    if not docs:
        return None

    doc = docs[0]

    # Pragmatic heuristic: try a few plausible text-ish fields.
    for key in ("curator_summary", "gene_group", "name"):
        value = doc.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip()
        if isinstance(value, list) and value:
            joined = "; ".join(str(v) for v in value if v)
            if joined.strip():
                return joined.strip()

    return None


def fetch_hgnc_function(hgnc_id: str) -> Optional[str]:
    """
    Fetch function-like text for an HGNC ID.
    hgnc_id is expected to look like 'HGNC:1234'.
    """
    if not hgnc_id:
        return None

    url = f"{HGNC_BASE}/fetch/hgnc_id/{hgnc_id}"
    headers = {"Accept": "application/json"}
    data = _http_get_json_with_retries(url, headers=headers)
    if not data:
        return None

    return _parse_hgnc_function(data)


# ------------------------ Aggregation logic ----------------------------------


def _normalize_text_for_dedup(text: Optional[str]) -> Optional[str]:
    """Normalize text for duplicate detection: lowercase, collapse whitespace."""
    if not text:
        return None
    norm = " ".join(str(text).split()).lower()
    return norm or None


def _value_or_na(value: Optional[str]) -> str:
    """Return a string value or a consistent placeholder when missing."""
    if value is None:
        return NA_PLACEHOLDER
    text = str(value).strip()
    return text if text else NA_PLACEHOLDER


def _aggregate_functions_for_row(
    ncbi_text: Optional[str],
    uniprot_text: Optional[str],
    hgnc_text: Optional[str],
) -> str:
    """
    Build the combined 'function' block with duplicate detection and labeling.

    Example:

    [SOURCE: NCBI] ...
    [SOURCE: UniProt] same as NCBI
    [SOURCE: HGNC] ...
    """
    sources = [
        ("NCBI", ncbi_text),
        ("UniProt", uniprot_text),
        ("HGNC", hgnc_text),
    ]

    norm_to_source: Dict[str, str] = {}
    lines = []

    for source_name, text in sources:
        if not text or _normalize_text_for_dedup(text) == NA_PLACEHOLDER:
            lines.append(f"[SOURCE: {source_name}] {NA_PLACEHOLDER}")
            continue

        norm = _normalize_text_for_dedup(text)
        if norm and norm in norm_to_source:
            original_source = norm_to_source[norm]
            lines.append(f"[SOURCE: {source_name}] same as {original_source}")
        else:
            if norm:
                norm_to_source[norm] = source_name
            lines.append(f"[SOURCE: {source_name}] {text}")

    return "\n".join(lines)


def run_function_annotation(gene_db: pd.DataFrame) -> pd.DataFrame:
    """
    Enrich the gene_db with:
      - function_ncbi
      - function_uniprot
      - function_hgnc
      - uniprot_entry_status
      - uniprot_protein_existence
      - function  (combined, multi-line string with duplicate labeling)

    Expects gene_db to contain:
      - summary_text                 (NCBI narrative summary)
      - uniprot_primary_accession    (UniProt accession from NCBI)
      - hgnc_id                      (HGNC ID from NCBI)
    """

    df = gene_db.copy()

    # Initialize columns if missing
    for col in [
        "function_ncbi",
        "function_uniprot",
        "function_hgnc",
        "uniprot_entry_status",
        "uniprot_protein_existence",
        "function",
    ]:
        if col not in df.columns:
            df[col] = None

    # Cache UniProt / HGNC requests so repeated accessions/IDs are cheap
    uniprot_cache: Dict[str, Dict[str, Optional[str]]] = {}
    hgnc_cache: Dict[str, Optional[str]] = {}

    for idx, row in df.iterrows():
        ncbi_text = row.get("summary_text")

        # ---- UniProt via accession ----
        uniprot_acc = row.get("uniprot_primary_accession")
        uni_function = None
        uni_status = None
        uni_pe = None

        if uniprot_acc:
            if uniprot_acc in uniprot_cache:
                uni_info = uniprot_cache[uniprot_acc]
            else:
                uni_info = fetch_uniprot_info(uniprot_acc)
                uniprot_cache[uniprot_acc] = uni_info

            uni_function = uni_info.get("function")
            uni_status = uni_info.get("entry_status")
            uni_pe = uni_info.get("protein_existence")

        # ---- HGNC best-effort ----
        hgnc_id = row.get("hgnc_id")
        hgnc_text = None
        if hgnc_id:
            if hgnc_id in hgnc_cache:
                hgnc_text = hgnc_cache[hgnc_id]
            else:
                hgnc_text = fetch_hgnc_function(hgnc_id)
                hgnc_cache[hgnc_id] = hgnc_text

        # Fill per-source columns
        ncbi_value = _value_or_na(ncbi_text)
        uni_value = _value_or_na(uni_function)
        hgnc_value = _value_or_na(hgnc_text)
        status_value = _value_or_na(uni_status)
        pe_value = _value_or_na(uni_pe)

        df.at[idx, "function_ncbi"] = ncbi_value
        df.at[idx, "function_uniprot"] = uni_value
        df.at[idx, "function_hgnc"] = hgnc_value
        df.at[idx, "uniprot_entry_status"] = status_value
        df.at[idx, "uniprot_protein_existence"] = pe_value

        # Build combined function block
        combined = _aggregate_functions_for_row(ncbi_value, uni_value, hgnc_value)
        df.at[idx, "function"] = combined

    return df
