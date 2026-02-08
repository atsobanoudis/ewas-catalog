import csv
import requests
import time
import json
import os

def get_publication_data(pmid):
    """
    Queries the EWAS Atlas API for a given PMID and returns the data as a JSON object.
    """
    url = f"https://ngdc.cncb.ac.cn/ewas/rest/publication?pmid={pmid}"
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        # Fix for Unicode apostrophe and other potential issues
        response_text = response.text.replace(r"\u2019", "'")
        return json.loads(response_text)
    except Exception as e:
        print(f"Error querying API for PMID {pmid}: {e}")
        return None

def main():
    input_file = "ewas_atlas.csv"
    
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return

    # Read the existing data
    with open(input_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)

    if not fieldnames:
        print("Error: CSV file has no headers.")
        return

    # Prepare new fieldnames
    # 'p' between 'correlation' and 'rank'
    # 'total_associations' after 'pmid'
    new_fieldnames = []
    for field in fieldnames:
        new_fieldnames.append(field)
        if field == "correlation":
            new_fieldnames.append("p")
    
    if "total_associations" not in new_fieldnames:
        # User said "after PMID", if PMID is last, it goes at the end.
        pmid_index = new_fieldnames.index("pmid") if "pmid" in new_fieldnames else -1
        if pmid_index != -1:
            new_fieldnames.insert(pmid_index + 1, "total_associations")
        else:
            new_fieldnames.append("total_associations")

    # Collect unique PMIDs to minimize API calls
    unique_pmids = sorted(list(set(row["pmid"] for row in rows if row.get("pmid") and row["pmid"].strip() and row["pmid"] != "NA")))
    
    print(f"Found {len(unique_pmids)} unique PMIDs to query.")
    
    pmid_cache = {}
    for pmid in unique_pmids:
        print(f"Fetching data for PMID {pmid}...")
        data = get_publication_data(pmid)
        if data and data.get("code") == 0:
            pmid_cache[pmid] = data.get("data")
        else:
            print(f"No data found for PMID {pmid}")
        time.sleep(0.5)  # Respectful delay

    # Process rows
    for row in rows:
        pmid = row.get("pmid")
        cpg = row.get("cpg")
        rank = row.get("rank")
        
        row["p"] = ""
        row["total_associations"] = ""
        
        if pmid and pmid in pmid_cache:
            pub_data = pmid_cache[pmid]
            study_list = pub_data.get("studyList", [])
            
            # 1. Total associations
            total = 0
            for study in study_list:
                # Handle potential typo in API field name 'assocaitionList'
                assoc_list = study.get("assocaitionList") or study.get("associationList") or []
                total += len(assoc_list)
            row["total_associations"] = total
            
            # 2. Find p-value
            found_p = None
            for study in study_list:
                assoc_list = study.get("assocaitionList") or study.get("associationList") or []
                for assoc in assoc_list:
                    # Match by CpG ID
                    if assoc.get("probeId") == cpg:
                        # If rank is provided in CSV, try to match it as well
                        if rank and rank.strip():
                            if str(assoc.get("rank")) == str(rank):
                                found_p = assoc.get("pvalue")
                                break
                        else:
                            # If no rank, take the first match
                            found_p = assoc.get("pvalue")
                            break
                if found_p is not None:
                    break
            
            row["p"] = found_p if found_p is not None else ""

    # Write updated data back to CSV
    with open(input_file, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=new_fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Successfully updated {input_file} with p-values and total associations.")

if __name__ == "__main__":
    main()
