import pandas as pd
import requests
import json
import time
import sys

# --- CONFIGURATION ---
GENOME_BUILD = "hg38"
GWAS_WINDOW_SIZE = 5000  # Look 5kb upstream and downstream for GWAS hits
INPUT_FILE = "ewas_res_groupsig_128.xlsx"
OUTPUT_FILE = "ucsc/ewas_ucsc_annotated.xlsx"

# Define your tracks separated by search strategy
TRACKS = {
    # STRATEGY 1: Exact Overlap (Did a mutation break my CpG?)
    "direct_overlap": [
        "clinvarMain",      # Clinical significance
        # "alphaMissense",    # Pathogenicity prediction
        # "snp155"            # Common SNPs (dbSNP)
    ],
    # STRATEGY 2: Neighborhood (Is this CpG near a disease driver?)
    "neighborhood": [
        "gwasCatalog"       # Disease Associations
    ]
}

def query_ucsc_track(chrom, start, end, track_name):
    """
    Queries the UCSC API for a specific track at specific coordinates.
    """
    base_url = "https://api.genome.ucsc.edu/getData/track"
    params = {
        "genome": GENOME_BUILD,
        "track": track_name,
        "chrom": chrom,
        "start": start,
        "end": end
    }
    
    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()
        
        if track_name in data:
            return data[track_name]
        return []
        
    except requests.exceptions.RequestException as e:
        print(f"  [!] Error querying {track_name}: {e}")
        return []

def verify_cpg(cpg_id, chrom, start, end):
    """
    Verifies that the coordinates resolve to the given CpG ID using the snpArrayIllumina850k track.
    """
    track_name = "snpArrayIllumina850k"
    hits = query_ucsc_track(chrom, start, end, track_name)
    
    for hit in hits:
        if hit.get('name') == cpg_id:
            return True
    return False

def flatten_result(cpg_id, original_row_data, hit_data):
    """
    Combines metadata with the hit data.
    """
    # Start with a base record identifying the CpG
    flat = {
        'cpg_id': cpg_id,
        'query_chrom': original_row_data['chr'],
        'query_start': original_row_data['Start_hg38'],
        'query_end': original_row_data['End_hg38']
    }
    # Update with the actual API result
    flat.update(hit_data)
    return flat

def main():
    print(f"Reading input file: {INPUT_FILE}...")
    try:
        df = pd.read_excel(INPUT_FILE)
    except FileNotFoundError:
        print(f"Error: Could not find {INPUT_FILE}")
        sys.exit(1)

    # Initialize storage for results (one list per track)
    all_tracks = [t for cat in TRACKS.values() for t in cat]
    results_storage = {track: [] for track in all_tracks}
    
    # Check if required columns exist
    required_cols = ['cpg', 'chr', 'Start_hg38', 'End_hg38']
    if not all(col in df.columns for col in required_cols):
        print(f"Error: Input file missing one of the required columns: {required_cols}")
        print(f"Found columns: {df.columns.tolist()}")
        sys.exit(1)

    print(f"Processing {len(df)} rows...")
    
    for index, row in df.iterrows():
        cpg_id = row['cpg']
        chrom = row['chr']
        start = row['Start_hg38']
        end = row['End_hg38']
        
        print(f"[{index+1}/{len(df)}] Processing {cpg_id}...", end="", flush=True)

        # 1. Verification
        if not verify_cpg(cpg_id, chrom, start, end):
            print(" [X] Verification Failed (Coordinates do not match CpG ID)")
            continue
        
        print(" [V] Verified. Running tracks...", end="", flush=True)
        
        # 2. Run Direct Overlap Tracks
        for track in TRACKS["direct_overlap"]:
            hits = query_ucsc_track(chrom, start, end, track)
            if hits:
                for hit in hits:
                    flat_hit = flatten_result(cpg_id, row, hit)
                    results_storage[track].append(flat_hit)
            else:
                # Add a blank record if no hits found to keep the CpG in the sheet
                results_storage[track].append(flatten_result(cpg_id, row, {}))

        # 3. Run Neighborhood Tracks
        win_start = max(0, start - GWAS_WINDOW_SIZE)
        win_end = end + GWAS_WINDOW_SIZE
        
        for track in TRACKS["neighborhood"]:
            hits = query_ucsc_track(chrom, win_start, win_end, track)
            if hits:
                for hit in hits:
                    flat_hit = flatten_result(cpg_id, row, hit)
                    results_storage[track].append(flat_hit)
            else:
                # Add a blank record if no hits found to keep the CpG in the sheet
                results_storage[track].append(flatten_result(cpg_id, row, {}))

        print(" Done.")
        # Be polite to the API
        time.sleep(0.1)

    # 4. Write Output
    print(f"Writing results to {OUTPUT_FILE}...")
    with pd.ExcelWriter(OUTPUT_FILE, engine='openpyxl') as writer:
        files_written = False
        for track, data in results_storage.items():
            if data:
                print(f"  - Writing sheet '{track}' with {len(data)} rows.")
                # Convert list of dicts to DataFrame
                track_df = pd.DataFrame(data)
                track_df.to_excel(writer, sheet_name=track, index=False)
                files_written = True
        
        if not files_written:
            print("  Warning: No results found for any track. Creating an empty file.")
            pd.DataFrame({'info': ['No hits found']}).to_excel(writer, sheet_name='Summary', index=False)

    print("Process complete!")

if __name__ == "__main__":
    main()