import csv
import requests
import json
import time

def get_cpg_data(cpg_id):
    """
    Queries the EWAS Atlas API for a given CpG ID and returns the data as a JSON object.
    """
    url = f"https://ngdc.cncb.ac.cn/ewas/rest/probe?probeId={cpg_id}"
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes
        # Fix for Unicode apostrophe
        response_text = response.text.replace(r"\u2019", "'")
        return json.loads(response_text)
    except requests.exceptions.RequestException as e:
        print(f"Error querying API for {cpg_id}: {e}")
        return None
def write_ewas_csv(data, filename="ewas_atlas.csv"):
    """
    Writes EWAS data to a CSV file.
    """
    with open(filename, "w", newline="", encoding="utf-8") as csvfile:
        fieldnames = ["cpg", "genes", "cpg_island", "trait", "correlation", "rank", "pmid"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for cpg_data in data:
            cpg = cpg_data.get("probeId")
            cpg_island = cpg_data.get("cpgIsland")
            
            genes = ";".join(sorted(list(set(rt.get("geneName") for rt in cpg_data.get("relatedTranscription", []) if rt.get("geneName")))))

            if not cpg_data.get("associationList"):
                writer.writerow({
                    "cpg": cpg,
                    "genes": genes,
                    "cpg_island": cpg_island,
                    "trait": None,
                    "correlation": None,
                    "rank": None,
                    "pmid": None
                })
            else:
                for association in cpg_data.get("associationList", []):
                    writer.writerow({
                        "cpg": cpg,
                        "genes": genes,
                        "cpg_island": cpg_island,
                        "trait": association.get("trait"),
                        "correlation": association.get("correlation"),
                        "rank": association.get("rank"),
                        "pmid": association.get("pmid")
                    })

def main():
    """
    Reads CpG IDs from a file, queries the EWAS Atlas API for each CpG,
    and writes the collected data to a JSON file and summarizes in csv.
    """
    with open("cpgs.txt", "r") as f:
        cpg_ids = [line.strip() for line in f if line.strip()]

    all_data = []
    for cpg_id in cpg_ids:
        print(f"Querying data for {cpg_id}...")
        cpg_data = get_cpg_data(cpg_id)
        if cpg_data and cpg_data.get("code") == 0:
            all_data.append(cpg_data["data"])
        time.sleep(0.5)  # delay

    with open("ewas_atlas.json", "w") as f:
        json.dump(all_data, f, indent=2)

    print("Successfully wrote EWAS data to ewas_atlas.json")
    
    write_ewas_csv(all_data)
    print("Successfully wrote EWAS data to ewas_atlas.csv")

if __name__ == "__main__":
    main()
