#!/usr/bin/env python3
"""Run the gene annotation pipeline and produce spreadsheet-friendly outputs."""

import json
import sys
from importlib import util as importlib_util
from pathlib import Path
from typing import Optional
import re
import logging

import pandas as pd

from modules.hgnc_pull import run_hgnc_pull
from modules.function_annotation import run_function_annotation
from modules.gwas_association import attach_psychiatric_gwas
from modules.pubmed_association import annotate_df_with_psych_literature
from modules.harmonizome_association import attach_harmonizome
from modules.cpg_integration import load_pi_cpg_mappings, attach_ewas_atlas_traits
from modules.disgenet_utils import annotate_with_disgenet


def load_ncbi_config(config_path: str = "config.json") -> tuple[Optional[str], Optional[str]]:
    """Load NCBI_EMAIL (required) and NCBI_API_KEY (optional) from config.json."""
    path = Path(config_path)
    if not path.is_file():
        print(f"[MAIN] Warning: config.json not found at {config_path}.")
        return None, None

    try:
        with path.open() as f:
            cfg = json.load(f)
    except Exception as exc:
        print(f"[MAIN] Error reading config.json: {exc}")
        return None, None

    email = cfg.get("NCBI_EMAIL")
    api_key = cfg.get("NCBI_API_KEY")

    if not email:
        print("[MAIN] Warning: NCBI_EMAIL missing in config.json.")
    return email, api_key


def _format_synonyms(val: object) -> Optional[str]:
    """Normalize synonyms to a ';'-separated string for readability."""
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return None

    if isinstance(val, list):
        cleaned = [str(v) for v in val if pd.notna(v)]
        return "; ".join(cleaned) if cleaned else None

    if isinstance(val, str):
        stripped = val.strip()
        try:
            parsed = json.loads(stripped)
            if isinstance(parsed, list):
                cleaned = [str(v) for v in parsed if pd.notna(v)]
                return "; ".join(cleaned) if cleaned else None
        except json.JSONDecodeError:
            pass
        return stripped

    return str(val)


def prompt_path_with_completion(prompt: str) -> str:
    """Prompt for a filesystem path with tab completion; degrade gracefully if unavailable."""
    # Try prompt_toolkit first (works cross-platform).
    try:
        from prompt_toolkit import prompt as pt_prompt
        from prompt_toolkit.completion import PathCompleter

        return pt_prompt(
            prompt,
            completer=PathCompleter(expanduser=True),
            complete_while_typing=True,
        )
    except Exception:
        pass

    try:
        import glob
        import os
        import readline
    except ImportError:
        return input(prompt)
    
    def complete_path(text: str, state: int) -> Optional[str]:
        expanded = os.path.expanduser(text)
        matches = glob.glob(expanded + '*')
        matches = [m + '/' if os.path.isdir(m) else m for m in matches]
        return matches[state] if state < len(matches) else None
    
    binding = "tab: complete"
    if getattr(readline, "__doc__", "") and "libedit" in readline.__doc__:
        binding = "bind ^I rl_complete"
    
    old_completer = readline.get_completer()
    old_delims = readline.get_completer_delims()
    readline.set_completer_delims(' \t\n;')
    readline.set_completer(complete_path)
    readline.parse_and_bind(binding)
    
    try:
        return input(prompt)
    finally:
        readline.set_completer(old_completer)
        readline.set_completer_delims(old_delims)


def export_readable_excel(gene_db: pd.DataFrame, output_path: str = "annotated_genes.xlsx") -> None:
    """Write a readable Excel version of the annotations with flattened lists and pruned columns."""
    df = gene_db.copy()

    if "synonyms" in df.columns:
        df["synonyms"] = df["synonyms"].apply(_format_synonyms)

    if "entrez" in df.columns:
        df["entrez"] = pd.to_numeric(df["entrez"], errors="coerce").astype("Int64")

    # Ensure GWAS association counts are numeric and default blanks to 0
    for col in ("gwas_assoc_count", "gwas_association_count"):
        if col in df.columns:
            df[col] = (
                pd.to_numeric(df[col], errors="coerce")
                .fillna(0)
                .astype("Int64")
            )

    # Present delimited lists as multi-line text for readability in Excel cells.
    def _to_multiline(text: object) -> object:
        if not isinstance(text, str):
            return text
        parts = [p.strip() for p in re.split(r";|\n", text) if p.strip()]
        return "\n".join(parts)

    list_cols = [
        "gwas_traits",
        "gwas_labels",
        "pubmed_terms",
        "pubmed_brief",
        "harmonizome_terms",
        "harmonizome_datasets",
    ]
    for col in list_cols:
        if col in df.columns:
            df[col] = df[col].apply(_to_multiline)

    # Add spacing between source blocks in the combined function field for readability.
    def _space_function_blocks(text: object) -> object:
        if not isinstance(text, str):
            return text
        lines = text.splitlines()
        spaced: list[str] = []
        for line in lines:
            stripped = line.strip()
            spaced.append(line.rstrip())
            # Treat any line that is exactly "n/a" or ends with " n/a" (case-insensitive) as empty
            is_na = stripped.lower() == "n/a" or stripped.lower().endswith(" n/a")
            if stripped and not is_na:
                spaced.append("")  # blank line after non-empty/non-NA entries
        # Remove trailing blank lines
        while spaced and spaced[-1] == "":
            spaced.pop()
        return "\n".join(spaced)

    if "function" in df.columns:
        df["function"] = df["function"].apply(_space_function_blocks)

    drop_cols = [
        "summary_text",
        "summary",
        "function_ncbi",
        "function_uniprot",
        "function_hgnc",
        "annotation_status",
        "annotation",
        "gwas_traits_example"
    ]
    df = df.drop(columns=[c for c in drop_cols if c in df.columns])

    if "name" in df.columns and "symbol" in df.columns:
        col_order = list(df.columns)
        col_order.remove("name")
        insert_at = col_order.index("symbol") + 1
        col_order.insert(insert_at, "name")
        df = df[col_order]

    # Pick an available Excel engine (cross-platform).
    engine = None
    if importlib_util.find_spec("openpyxl"):
        engine = "openpyxl"
    elif importlib_util.find_spec("xlsxwriter"):
        engine = "xlsxwriter"

    if engine is None:
        print("[MAIN] Skipping Excel export: no Excel writer engine found. Install `openpyxl` (recommended) or `XlsxWriter`, then rerun to produce the .xlsx file.")
        return

    try:
        df.to_excel(output_path, index=False, engine=engine)
        print(f"[MAIN] Readable Excel written to: {output_path} (engine: {engine})")
    except ModuleNotFoundError:
        # Rare race: engine discovered but import failed; tell user to install.
        print(f"[MAIN] Skipping Excel export: the {engine} package is not available at runtime. Install it and rerun to produce the .xlsx file.")


def build_output_df(gene_db: pd.DataFrame) -> pd.DataFrame:
    """Apply final column renames/order for deliverables and drop noisy intermediate fields."""
    df = gene_db.copy()
    rename_map = {
        "raw_input": "input",
        "approved_symbol": "symbol",
        "description": "name",
        "hgnc_id": "hgnc",
        "entrez_id": "entrez",
        "ensembl_id": "ensembl",
        "gene_type_raw": "gene_type",
        "uniprot_primary_accession": "uniprot",
        "uniprot_entry_status": "uniprot_review",
        "uniprot_protein_existence": "uniprot_existence",
        "Start_hg38": "cpg_start",
        "End_hg38": "cpg_end",
        "chr": "cpg_chr"
    }
    df = df.rename(columns=rename_map)

    if "synonyms" in df.columns:
        df["synonyms"] = df["synonyms"].apply(_format_synonyms)

    drop_cols = [
        "summary_text",
        "summary",
        "function_ncbi",
        "function_uniprot",
        "function_hgnc",
        "annotation_status",
        "annotation",
        "unique_gene_name"
    ]
    df = df.drop(columns=[c for c in drop_cols if c in df.columns])

    if "entrez" in df.columns:
        df["entrez"] = pd.to_numeric(df["entrez"], errors="coerce").astype("Int64")

    # Reorder to put CpG info near the front
    cols = list(df.columns)
    cpg_cols = ['cpg', 'cpg_chr', 'cpg_start', 'cpg_end']
    for c in reversed(cpg_cols):
        if c in cols:
            cols.remove(c)
            cols.insert(0, c)

    if "name" in cols and "symbol" in cols:
        cols.remove("name")
        insert_at = cols.index("symbol") + 1
        cols.insert(insert_at, "name")
    
    df = df[cols]

    return df



def main() -> None:
    """Run the annotation pipeline end to end."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    
    # Get gene list file path from command line or prompt user
    # Default to the PI's CpG mapping file
    default_input = "ewas_res_groupsig_128.xlsx"
    if len(sys.argv) > 1:
        input_path = sys.argv[1]
    else:
        input_path = prompt_path_with_completion(f"Enter path to input file (default: {default_input}): ").strip()
        if not input_path:
            input_path = default_input
    
    # Validate file exists
    input_path = str(Path(input_path).expanduser())
    if not Path(input_path).is_file():
        print(f"Error: File not found: {input_path}", file=sys.stderr)
        sys.exit(1)
    
    print(f"\n{'='*60}")
    print("Gene Annotation Pipeline")
    print(f"{'='*60}")
    print(f"Input file: {input_path}\n")
    
    # Determine if input is the CpG mapping Excel
    if input_path.endswith('.xlsx') and 'ewas_res' in input_path:
        print("[MAIN] Detected CpG mapping file. Loading PI mappings...")
        cpg_db = load_pi_cpg_mappings(input_path)
        # We need to annotate these genes
        unique_genes = cpg_db['gene'].unique().tolist()
        print(f"[MAIN] Pulling metadata for {len(unique_genes)} unique genes from CpG mapping...")
        gene_db_raw = run_hgnc_pull(unique_genes)
        # Merge back
        gene_db = pd.merge(cpg_db, gene_db_raw, left_on='gene', right_on='raw_input', how='left')
    else:
        # Run HGNC/NCBI annotation module from a standard gene list
        print("[MAIN] Starting HGNC/NCBI annotation module...")
        gene_db: pd.DataFrame = run_hgnc_pull(input_path)


    # Run function aggregation module (NCBI + UniProt + HGNC)
    print("\n[MAIN] Starting function annotation (NCBI + UniProt + HGNC)...")
    gene_db = run_function_annotation(gene_db)

    # ------------------------------------------------------------------
    # CpG-level EWAS Atlas associations
    # ------------------------------------------------------------------
    if 'cpg' in gene_db.columns:
        print("\n[MAIN] Attaching EWAS Atlas CpG-trait associations...")
        atlas_path = "data/ewas_atlas.csv"
        if Path(atlas_path).is_file():
            atlas_df = pd.read_csv(atlas_path, low_memory=False)
            gene_db = attach_ewas_atlas_traits(gene_db, atlas_df)
        else:
            print(f"[MAIN] Warning: {atlas_path} not found; skipping Atlas step.")

    # ------------------------------------------------------------------
    # Psychiatric GWAS Catalog associations (gene-level)
    # ------------------------------------------------------------------
    print("\n[MAIN] Attaching GWAS psychiatric associations (gene-level)...")

    # Use approved_symbol if available; otherwise fall back to symbol.
    symbol_col = None
    if "approved_symbol" in gene_db.columns:
        symbol_col = "approved_symbol"
    elif "symbol" in gene_db.columns:
        symbol_col = "symbol"

    if symbol_col is not None:
        gene_db = attach_psychiatric_gwas(gene_db, gene_col=symbol_col)
    else:
        print("[MAIN] No symbol / approved_symbol column found; skipping GWAS step.")

    # Harmonizome psychiatric disease associations
    print("\n[MAIN] Attaching Harmonizome-based psychiatric disease evidence...")
    gene_db = attach_harmonizome(gene_db, symbol_col="approved_symbol")
    
    # ------------------------------------------------------------------
    # DisGeNET Broad-Spectrum associations
    # ------------------------------------------------------------------
    print("\n[MAIN] Attaching DisGeNET associations (Broad-Spectrum & Psych)...")
    gea_path = 'data/disgenet_gea.csv'
    if Path(gea_path).is_file():
        disgenet_df = pd.read_csv(gea_path)
        # 1. Broad spectrum (all traits)
        gene_db = annotate_with_disgenet(gene_db, disgenet_df, psych_only=False)
        # 2. Psych specific (to maintain compatibility with previous outputs)
        gene_db = annotate_with_disgenet(gene_db, disgenet_df, psych_only=True)
    else:
        print(f"[MAIN] Warning: {gea_path} not found; skipping DisGeNET step.")

    # ------------------------------------------------------------------
    # PubMed mental-health literature (MeSH/text-based, gene-level)
    # ------------------------------------------------------------------
    print("\n[MAIN] Attaching PubMed mental-health literature associations...")

    ncbi_email, ncbi_api_key = load_ncbi_config()

    if not ncbi_email:
        print("[MAIN] NCBI_EMAIL missing; skipping PubMed psych literature step.")
    else:
        gene_db = annotate_df_with_psych_literature(
            gene_db,
            gene_symbol_col="approved_symbol" if "approved_symbol" in gene_db.columns else "symbol",
            entrez_col="entrez_id",
            email=ncbi_email,
            api_key=ncbi_api_key,
        )

    # Print summary
    print(f"\n{'='*60}")
    print("Annotation Summary")
    print(f"{'='*60}")
    print(f"Total genes processed: {len(gene_db)}")
    
    # Status breakdown
    if 'annotation_status' in gene_db.columns:
        print("\nAnnotation Status Breakdown:")
        status_counts = gene_db['annotation_status'].value_counts()
        for status, count in status_counts.items():
            print(f"  {status}: {count}")

    # Write full, un-pruned results to CSV for debugging/traceability
    raw_output_path = "annotated_genes_raw.csv"
    gene_db.to_csv(raw_output_path, index=False)
    print(f"\n[MAIN] Raw results written to: {raw_output_path} (all columns retained)")
    
    # Prepare deliverable with final column names/order
    output_df = build_output_df(gene_db)

    # Show first few rows (deliverable view)
    print(f"\n{'='*60}")
    print("First 5 Annotated Genes (deliverable columns):")
    print(f"{'='*60}")
    display_cols = [
        'cpg',
        'symbol',
        'name',
        'ewas_atlas_traits',
        'disgenet_diseases',
        'pubmed_count',
        'gwas_assoc_count',
    ]
    available_cols = [col for col in display_cols if col in output_df.columns]
    print(output_df[available_cols].head().to_string(index=False))

    
    export_readable_excel(output_df)
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
