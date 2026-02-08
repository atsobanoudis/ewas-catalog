"""Annotate genes with mental-health PubMed hits using Entrez links or symbol searches."""

from __future__ import annotations

import json
import logging
import re
import time
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import requests
import pandas as pd
from xml.etree import ElementTree as ET

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# CONFIG / CONSTANTS
# ---------------------------------------------------------------------------

NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Default mental health MeSH headings to look for.
# These should be exact MeSH Heading or Qualifier names as they appear in PubMed XML.
DEFAULT_MENTAL_HEALTH_MESH_TERMS: List[str] = [
    # core mood disorders
    "Depressive Disorder",
    "Depressive Disorder, Major",
    "Bipolar Disorder",
    "Cyclothymic Disorder",
    # psychotic disorders
    "Schizophrenia",
    "Psychotic Disorders",
    "Schizoaffective Disorder",
    # anxiety, trauma, stress
    "Anxiety Disorders",
    "Panic Disorder",
    "Phobic Disorders",
    "Obsessive-Compulsive Disorder",
    "Stress Disorders, Post-Traumatic",
    # neurodevelopmental
    "Autism Spectrum Disorder",
    "Attention Deficit Disorder with Hyperactivity",
    "Intellectual Disability",
    # substance use
    "Substance-Related Disorders",
    "Alcohol-Related Disorders",
    "Opioid-Related Disorders",
    "Cocaine-Related Disorders",
    # self-harm/suicide
    "Suicide",
    "Suicidal Ideation",
    # broad
    "Mental Disorders",
    "Mental Health",
]

# Text keywords to catch mental-health context even if MeSH is incomplete.
DEFAULT_MENTAL_HEALTH_TEXT_TERMS: List[str] = [
    "schizophrenia",
    "bipolar",
    "depressi",      # depression, depressive
    "mania",
    "psychosis",
    "psychotic",
    "autism",
    "asperger",
    "adhd",
    "attention-deficit",
    "anxiety",
    "panic disorder",
    "ptsd",
    "post-traumatic",
    "suicide",
    "suicidal",
    "substance use",
    "substance-use",
    "addiction",
    "alcohol use",
    "alcohol-use",
    "mental disorder",
    "mental health",
]

# Patterns to label a hit as "genetic/variant" related.
DEFAULT_GENETIC_TEXT_PATTERNS: List[re.Pattern] = [
    re.compile(r"\bpolymorphism\b", re.I),
    re.compile(r"\bpolymorphisms\b", re.I),
    re.compile(r"\bmutation\b", re.I),
    re.compile(r"\bmutations\b", re.I),
    re.compile(r"\bvariant\b", re.I),
    re.compile(r"\bvariants\b", re.I),
    re.compile(r"\bgenetic\b", re.I),
    re.compile(r"\bgenome[- ]wide association\b", re.I),
    re.compile(r"\bgwas\b", re.I),
    re.compile(r"\bsnp\b", re.I),
    re.compile(r"\bcopy number variation\b", re.I),
    re.compile(r"\bcnv\b", re.I),
]

# MeSH terms that strongly suggest a genetic/variant study.
DEFAULT_GENETIC_MESH_TERMS: List[str] = [
    "Polymorphism, Genetic",
    "Genetic Variation",
    "Genome-Wide Association Study",
    "Genotype",
    "Genetic Predisposition to Disease",
]


# ---------------------------------------------------------------------------
# DATA STRUCTURES
# ---------------------------------------------------------------------------

@dataclass
class PubMedRecord:
    pmid: str
    title: str
    year: Optional[int]
    mesh_terms: List[str]
    is_mental_health: bool
    is_genetic: bool


# ---------------------------------------------------------------------------
# LOW-LEVEL EUTILS HELPERS
# ---------------------------------------------------------------------------

def _eutils_get(
    endpoint: str,
    params: Dict[str, str],
    email: str,
    api_key: Optional[str] = None,
    pause: float = 0.34,
    max_retries: int = 3,
    backoff_factor: float = 1.5,
) -> str:
    """
    Helper for GET requests to NCBI E-utilities with rate limiting and
    sensible defaults.

    Parameters
    ----------
    endpoint : str
        e.g. "elink.fcgi", "esummary.fcgi", "efetch.fcgi".
    params : dict
        Query parameters excluding 'email' and 'api_key', which are added here.
    email : str
        Your email address (NCBI requirement).
    api_key : str, optional
        NCBI API key to increase rate limits.
    pause : float
        Seconds to sleep after each request (NCBI suggests <= 3 req/sec without key).

    Returns
    -------
    str
        Response text (XML or similar), raise for HTTP errors.
    """
    url = f"{NCBI_EUTILS_BASE}/{endpoint}"
    merged_params = dict(params)
    merged_params["email"] = email
    if api_key:
        merged_params["api_key"] = api_key

    safe_params = dict(merged_params)
    if "api_key" in safe_params:
        safe_params["api_key"] = "***"
    logger.info("[PUBMED] GET %s with params %s", endpoint, safe_params)

    headers = {
        "Connection": "close",           # avoid keep-alive issues on large payloads
        "Accept-Encoding": "identity",   # disable gzip/chunked compression to reduce truncation errors
    }

    last_exc: Optional[Exception] = None
    for attempt in range(1, max_retries + 1):
        try:
            resp = requests.get(
                url,
                params=merged_params,
                headers=headers,
                timeout=45,
            )
            logger.info(
                "[PUBMED] %s -> status %s (%d chars)",
                endpoint,
                resp.status_code,
                len(resp.text),
            )
            resp.raise_for_status()
            time.sleep(pause)
            return resp.text
        except (
            requests.exceptions.ChunkedEncodingError,
            requests.exceptions.ProtocolError,
        ) as exc:
            last_exc = exc
            logger.warning(
                "[PUBMED] %s response ended prematurely (attempt %d/%d): %s",
                endpoint,
                attempt,
                max_retries,
                exc,
            )
        except requests.exceptions.RequestException as exc:
            last_exc = exc
            logger.warning(
                "[PUBMED] %s request failed (attempt %d/%d): %s",
                endpoint,
                attempt,
                max_retries,
                exc,
            )

        sleep_for = backoff_factor * attempt
        logger.info("[PUBMED] sleeping %.1fs before retrying %s", sleep_for, endpoint)
        time.sleep(sleep_for)

    raise RuntimeError(
        f"Failed to GET {endpoint} after {max_retries} attempts"
    ) from last_exc


def elink_gene_to_pubmed(
    entrez_gene_id: int,
    email: str,
    api_key: Optional[str] = None,
    pause: float = 0.34,
    max_ids: Optional[int] = None,
) -> List[str]:
    """
    Use ELink to retrieve PubMed IDs linked to an Entrez Gene ID.

    This leverages NCBI's gene2pubmed mapping and avoids the
    NTM/tuberculosis ambiguity you ran into with symbol-only searches.

    Parameters
    ----------
    entrez_gene_id : int
        Entrez Gene ID.
    email : str
        For NCBI.
    api_key : str, optional
        NCBI API key.
    pause : float
        Seconds to sleep after request.
    max_ids : int, optional
        Optional cap on number of PMIDs returned.

    Returns
    -------
    List[str]
        List of PubMed IDs as strings.
    """
    xml_text = _eutils_get(
        "elink.fcgi",
        params={
            "dbfrom": "gene",
            "db": "pubmed",
            "id": str(entrez_gene_id),
            "retmode": "xml",
        },
        email=email,
        api_key=api_key,
        pause=pause,
    )

    root = ET.fromstring(xml_text)
    pmids: List[str] = []

    for linkset in root.findall(".//LinkSetDb"):
        db = linkset.findtext("DbTo")
        if db != "pubmed":
            continue
        for link in linkset.findall("Link"):
            pmid = link.findtext("Id")
            if pmid:
                pmids.append(pmid)

    if max_ids is not None:
        pmids = pmids[:max_ids]

    logger.debug("ELink gene %s -> %d PMIDs", entrez_gene_id, len(pmids))
    return pmids


def esearch_pubmed_by_symbol(
    gene_symbol: str,
    email: str,
    api_key: Optional[str] = None,
    pause: float = 0.34,
    max_ids: int = 200,
) -> List[str]:
    """
    Fallback: PubMed ESearch using gene symbol. This is inherently noisier
    than gene->pubmed links, but useful when Entrez ID is missing.

    We keep this deliberately simple and *do not* bake mental health or
    genetic filters into the query; those are handled later based on
    metadata to keep things transparent.

    Parameters
    ----------
    gene_symbol : str
        Official gene symbol.
    email : str
        For NCBI.
    api_key : str, optional
        NCBI API key.
    pause : float
        Seconds to sleep after request.
    max_ids : int
        Maximum PMIDs to retrieve.

    Returns
    -------
    List[str]
        PubMed IDs as strings.
    """
    # Simple query: gene symbol in title/abstract OR MeSH.
    # You can customize this if needed.
    query = f'"{gene_symbol}"[Title/Abstract] OR "{gene_symbol}"[MeSH Terms]'

    xml_text = _eutils_get(
        "esearch.fcgi",
        params={
            "db": "pubmed",
            "term": query,
            "retmode": "xml",
            "retmax": str(max_ids),
        },
        email=email,
        api_key=api_key,
        pause=pause,
    )

    root = ET.fromstring(xml_text)
    pmids = [e.text for e in root.findall(".//IdList/Id") if e.text]
    logger.debug("ESearch symbol %s -> %d PMIDs", gene_symbol, len(pmids))
    return pmids


def efetch_pubmed_records(
    pmids: Sequence[str],
    email: str,
    api_key: Optional[str] = None,
    pause: float = 0.34,
    batch_size: int = 200,
) -> List[PubMedRecord]:
    """
    Fetch PubMed records as XML and parse key fields.

    Parameters
    ----------
    pmids : sequence of str
        PubMed IDs to fetch.
    email : str
        For NCBI.
    api_key : str, optional
        NCBI API key.
    pause : float
        Seconds to sleep between batches.
    batch_size : int
        Number of PMIDs per EFetch call.

    Returns
    -------
    List[PubMedRecord]
        Parsed PubMedRecord instances.
    """
    def _fetch_batches(batch_pmids: Sequence[str], current_batch_size: int) -> List[PubMedRecord]:
        records: List[PubMedRecord] = []
        batch_pmids = list(batch_pmids)

        for i in range(0, len(batch_pmids), current_batch_size):
            batch = batch_pmids[i : i + current_batch_size]
            if not batch:
                continue

            logger.info(
                "[PUBMED] Fetching batch %d-%d of %d PMIDs (size=%d)",
                i + 1,
                min(i + current_batch_size, len(batch_pmids)),
                len(batch_pmids),
                current_batch_size,
            )

            try:
                xml_text = _eutils_get(
                    "efetch.fcgi",
                    params={
                        "db": "pubmed",
                        "id": ",".join(batch),
                        "retmode": "xml",
                    },
                    email=email,
                    api_key=api_key,
                    pause=pause,
                )
            except Exception as exc:
                # Large batches can occasionally be truncated; split and retry smaller chunks.
                if len(batch) > 1 and current_batch_size > 25:
                    smaller = max(25, current_batch_size // 2)
                    logger.warning(
                        "[PUBMED] Batch %d-%d failed (%s); retrying with smaller chunk size %d",
                        i + 1,
                        min(i + current_batch_size, len(batch_pmids)),
                        exc,
                        smaller,
                    )
                    records.extend(_fetch_batches(batch, smaller))
                    continue

                logger.error(
                    "[PUBMED] Failed to fetch PMIDs %s (batch size=%d): %s",
                    ",".join(batch),
                    current_batch_size,
                    exc,
                )
                continue

            root = ET.fromstring(xml_text)
            for article in root.findall(".//PubmedArticle"):
                try:
                    record = _parse_pubmed_article(article)
                    records.append(record)
                except Exception as exc:  # be defensive
                    logger.warning("Failed to parse PubMed article: %s", exc)

        return records

    return _fetch_batches(pmids, batch_size)


def _parse_pubmed_article(node: ET.Element) -> PubMedRecord:
    """Parse a single <PubmedArticle> node into PubMedRecord."""
    pmid = node.findtext(".//MedlineCitation/PMID") or ""

    # Title
    title = node.findtext(".//Article/ArticleTitle") or ""
    title = title.strip()

    # Year (try ArticleDate, then PubDate Year, then MedlineDate)
    year: Optional[int] = None
    yr_text = None

    # Try <ArticleDate>
    article_date = node.find(".//Article/ArticleDate/Year")
    if article_date is not None and article_date.text:
        yr_text = article_date.text
    else:
        # Try <PubDate><Year>
        pub_date_year = node.find(".//Article/Journal/JournalIssue/PubDate/Year")
        if pub_date_year is not None and pub_date_year.text:
            yr_text = pub_date_year.text
        else:
            # Fallback: parse MedlineDate like '2007 Jan-Feb'
            medline_date = node.findtext(
                ".//Article/Journal/JournalIssue/PubDate/MedlineDate"
            )
            if medline_date:
                match = re.search(r"(\d{4})", medline_date)
                if match:
                    yr_text = match.group(1)

    if yr_text:
        try:
            year = int(yr_text)
        except ValueError:
            year = None

    # MeSH headings
    mesh_terms: List[str] = []
    for mh in node.findall(".//MeshHeadingList/MeshHeading/DescriptorName"):
        if mh.text:
            mesh_terms.append(mh.text)

    # We'll compute is_mental_health and is_genetic later to keep this
    # parser generic. Placeholders here:
    return PubMedRecord(
        pmid=pmid,
        title=title,
        year=year,
        mesh_terms=mesh_terms,
        is_mental_health=False,
        is_genetic=False,
    )


# ---------------------------------------------------------------------------
# CLASSIFICATION HELPERS
# ---------------------------------------------------------------------------

def is_mental_health_record(
    record: PubMedRecord,
    mesh_terms: Optional[Iterable[str]] = None,
    text_terms: Optional[Iterable[str]] = None,
) -> Tuple[bool, List[str]]:
    """
    Determine whether a PubMedRecord is mental-health-related, and
    return (is_mental_health, matched_terms).

    Parameters
    ----------
    record : PubMedRecord
    mesh_terms : iterable of str, optional
        MeSH terms to match against DescriptorNames.
    text_terms : iterable of str, optional
        Lowercased substrings to search in the title.

    Returns
    -------
    (bool, List[str])
        is_mental_health, matched_terms
    """
    mesh_terms = list(mesh_terms or DEFAULT_MENTAL_HEALTH_MESH_TERMS)
    text_terms = list(text_terms or DEFAULT_MENTAL_HEALTH_TEXT_TERMS)

    matched: List[str] = []

    # Check MeSH
    mesh_lower = {m.lower() for m in record.mesh_terms}
    for mt in mesh_terms:
        if mt.lower() in mesh_lower:
            matched.append(mt)

    # Check title text
    title_lower = (record.title or "").lower()
    for tt in text_terms:
        if tt in title_lower:
            matched.append(tt)

    return (len(matched) > 0), matched


def is_genetic_record(
    record: PubMedRecord,
    genetic_mesh_terms: Optional[Iterable[str]] = None,
    genetic_patterns: Optional[Iterable[re.Pattern]] = None,
) -> bool:
    """
    Tag whether a PubMedRecord looks "genetic/variant" related.

    This is *not* used to filter the mental-health hits by default; it's
    an annotation that lets you later subset or stratify.

    Parameters
    ----------
    record : PubMedRecord
    genetic_mesh_terms : iterable of str, optional
        MeSH terms signalling genetic studies.
    genetic_patterns : iterable of compiled regex, optional
        Patterns to scan in the title.

    Returns
    -------
    bool
        True if likely a genetic/variant paper.
    """
    genetic_mesh_terms = list(genetic_mesh_terms or DEFAULT_GENETIC_MESH_TERMS)
    genetic_patterns = list(genetic_patterns or DEFAULT_GENETIC_TEXT_PATTERNS)

    # MeSH-based
    mesh_lower = {m.lower() for m in record.mesh_terms}
    for gm in genetic_mesh_terms:
        if gm.lower() in mesh_lower:
            return True

    # Title-based
    title = record.title or ""
    for pattern in genetic_patterns:
        if pattern.search(title):
            return True

    return False


# ---------------------------------------------------------------------------
# GENE-LEVEL ANNOTATION
# ---------------------------------------------------------------------------

def fetch_psych_pubmed_for_gene(
    gene_symbol: Optional[str],
    entrez_gene_id: Optional[int],
    email: str,
    api_key: Optional[str] = None,
    pause: float = 0.34,
    max_pmids_per_gene: int = 500,
) -> List[PubMedRecord]:
    """
    Fetch PubMed records for a single gene and tag mental health / genetic.

    Priority:
      1. Use gene->pubmed link if Entrez ID available.
      2. Otherwise fall back to text search by gene symbol.

    Parameters
    ----------
    gene_symbol : str or None
        Approved gene symbol.
    entrez_gene_id : int or None
        Entrez Gene ID.
    email : str
        NCBI email.
    api_key : str, optional
        NCBI API key.
    pause : float
        Rate limiting pause.
    max_pmids_per_gene : int
        Maximum PMIDs to fetch per gene (after link/search).

    Returns
    -------
    List[PubMedRecord]
        All PubMed records (not yet filtered), with mental/genetic flags set.
    """
    pmids: List[str] = []

    logger.info(
        "[PUBMED] Processing gene %s (Entrez: %s)",
        gene_symbol or "<no symbol>",
        entrez_gene_id if entrez_gene_id is not None else "n/a",
    )

    if entrez_gene_id is not None:
        try:
            pmids = elink_gene_to_pubmed(
                entrez_gene_id=entrez_gene_id,
                email=email,
                api_key=api_key,
                pause=pause,
                max_ids=max_pmids_per_gene,
            )
        except Exception as exc:
            logger.warning(
                "ELink gene->pubmed failed for Entrez %s: %s", entrez_gene_id, exc
            )

    # Fallback if no PMIDs found via Entrez or Entrez missing
    if not pmids and gene_symbol:
        try:
            pmids = esearch_pubmed_by_symbol(
                gene_symbol=gene_symbol,
                email=email,
                api_key=api_key,
                pause=pause,
                max_ids=max_pmids_per_gene,
            )
        except Exception as exc:
            logger.warning(
                "ESearch symbol->pubmed failed for %s: %s", gene_symbol, exc
            )

    if not pmids:
        logger.info("[PUBMED] No PMIDs returned for %s.", gene_symbol or "<no symbol>")
        return []

    logger.info("[PUBMED] Retrieved %d PMIDs for %s; fetching details...", len(pmids), gene_symbol or "<no symbol>")
    records = efetch_pubmed_records(
        pmids=pmids,
        email=email,
        api_key=api_key,
        pause=pause,
    )

    # Tag mental health and genetic flags
    for i, rec in enumerate(records):
        is_mh, _ = is_mental_health_record(rec)
        is_gen = is_genetic_record(rec)
        records[i] = PubMedRecord(
            pmid=rec.pmid,
            title=rec.title,
            year=rec.year,
            mesh_terms=rec.mesh_terms,
            is_mental_health=is_mh,
            is_genetic=is_gen,
        )

    return records


def annotate_df_with_psych_literature(
    df: pd.DataFrame,
    *,
    gene_symbol_col: str = "approved_symbol",
    entrez_col: str = "entrez_id",
    email: str,
    api_key: Optional[str] = None,
    pause: float = 0.34,
    max_pmids_per_gene: int = 500,
    mental_health_mesh_terms: Optional[Iterable[str]] = None,
    mental_health_text_terms: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    """
    Main entry point: annotate a DataFrame of genes with columns describing
    their mental-health-related PubMed hits.

    New columns added
    -----------------
    'pubmed_count' : int
        Number of mental-health-related PubMed hits for the gene.
    'pubmed_genetic_count' : int
        Subset of hits that look "genetic/variant" related.
    'pubmed_pmids' : str
        Semicolon-delimited PMIDs for mental-health-related hits.
    'pubmed_terms' : str
        Semicolon-delimited mental-health MeSH/text terms that matched across hits.
    'pubmed_brief' : str
        Semicolon-delimited summaries per hit in the order: title, year, PMID, genetic label.

    Parameters
    ----------
    df : pd.DataFrame
        Must have columns for gene symbol and/or Entrez Gene ID.
    gene_symbol_col : str
        Column name for approved gene symbol.
    entrez_col : str
        Column name for Entrez Gene ID (int or str convertible).
    email : str
        email required by NCBI.
    api_key : str, optional
        NCBI API key for higher rate limits.
    pause : float
        Seconds to sleep between e-utils requests.
    max_pmids_per_gene : int
        Hard cap per gene for PMIDs.
    mental_health_mesh_terms : iterable of str, optional
        Override default MeSH terms for mental health.
    mental_health_text_terms : iterable of str, optional
        Override default text terms.

    Returns
    -------
    pd.DataFrame
        Copy of df with appended annotation columns.
    """
    mesh_terms = list(mental_health_mesh_terms or DEFAULT_MENTAL_HEALTH_MESH_TERMS)
    text_terms = list(mental_health_text_terms or DEFAULT_MENTAL_HEALTH_TEXT_TERMS)

    df = df.copy()

    # Initialize columns
    df["pubmed_count"] = 0
    df["pubmed_genetic_count"] = 0
    df["pubmed_pmids"] = ["" for _ in range(len(df))]
    df["pubmed_terms"] = ["" for _ in range(len(df))]
    df["pubmed_brief"] = ["" for _ in range(len(df))]

    total_genes = len(df)
    for idx, row in df.iterrows():
        gene_symbol = row.get(gene_symbol_col)
        entrez_val = row.get(entrez_col)
        logger.info("[PUBMED] (%d/%d) Annotating %s", idx + 1, total_genes, gene_symbol or "<no symbol>")

        # normalize Entrez ID
        entrez_id: Optional[int] = None
        if pd.notna(entrez_val):
            try:
                entrez_id = int(entrez_val)
            except Exception:
                entrez_id = None

        records = fetch_psych_pubmed_for_gene(
            gene_symbol=gene_symbol,
            entrez_gene_id=entrez_id,
            email=email,
            api_key=api_key,
            pause=pause,
            max_pmids_per_gene=max_pmids_per_gene,
        )

        # Filter to mental health hits using current meshes/texts
        mh_records: List[PubMedRecord] = []
        combined_terms: List[str] = []

        for rec in records:
            is_mh, matched_terms = is_mental_health_record(
                rec,
                mesh_terms=mesh_terms,
                text_terms=text_terms,
            )
            if is_mh:
                mh_records.append(rec)
                combined_terms.extend(matched_terms)

        if not mh_records:
            logger.info("[PUBMED] -> No mental-health hits for %s", gene_symbol or "<no symbol>")
            continue

        # Aggregate
        pmids = [r.pmid for r in mh_records]
        # Unique, sorted terms for readability
        unique_terms = sorted(set(combined_terms))

        # Human-readable strings; easier to scan in spreadsheets.
        genetic_hits = [r for r in mh_records if r.is_genetic]
        df.at[idx, "pubmed_count"] = len(mh_records)
        df.at[idx, "pubmed_genetic_count"] = len(genetic_hits)
        df.at[idx, "pubmed_pmids"] = "; ".join(pmids)
        df.at[idx, "pubmed_terms"] = "; ".join(unique_terms)
        brief_entries = []
        for r in mh_records:
            genetic_label = "Genetic" if r.is_genetic else "Not Genetic"
            year_part = str(r.year) if r.year is not None else ""
            parts = [p for p in (r.title, year_part, f"PMID: {r.pmid}", genetic_label) if p]
            brief_entries.append(" ".join(parts))
        df.at[idx, "pubmed_brief"] = "; ".join(brief_entries)
        logger.info(
            "[PUBMED] -> %d mental-health hits (%d genetic) for %s",
            len(mh_records),
            len(genetic_hits),
            gene_symbol or "<no symbol>",
        )

    return df

