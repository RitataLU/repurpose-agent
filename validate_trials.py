"""
validate_trials.py — Query ClinicalTrials.gov v2 for the top 20 orphan genes.
Usage: python3 validate_trials.py
Output: docs/trials.json
"""
import json, time, urllib.request, urllib.parse, urllib.error

GENES_JSON  = "docs/genes.json"
TRIALS_JSON = "docs/trials.json"

ACTIVE_STATUSES = {"RECRUITING", "ACTIVE_NOT_RECRUITING"}

PHASE_ORDER = {"PHASE4": 4, "PHASE3": 3, "PHASE2": 2, "PHASE1": 1, "NA": 0, "EARLY_PHASE1": 0}

def fetch_trials(gene):
    params = urllib.parse.urlencode({
        "query.intr": gene,
        "pageSize": 20,
        "fields": "NCTId,OverallStatus,Phase,StartDate",
        "format": "json",
    })
    url = f"https://clinicaltrials.gov/api/v2/studies?{params}"
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            return json.loads(resp.read().decode())
    except Exception as e:
        print(f"  [error] {gene}: {e}")
        return None

def parse_response(data):
    studies = data.get("studies", [])
    n_active = 0
    max_phase_rank = 0
    max_phase_label = None
    latest_start = None

    for study in studies:
        ps = study.get("protocolSection", {})
        status_mod = ps.get("statusModule", {})
        status = status_mod.get("overallStatus", "")
        if status in ACTIVE_STATUSES:
            n_active += 1

        # Phase
        design_mod = ps.get("designModule", {})
        phases = design_mod.get("phases", [])
        for ph in phases:
            rank = PHASE_ORDER.get(ph, 0)
            if rank > max_phase_rank:
                max_phase_rank = rank
                max_phase_label = ph

        # Start date
        start = status_mod.get("startDateStruct", {}).get("date", "")
        if start:
            if latest_start is None or start > latest_start:
                latest_start = start

    return n_active, max_phase_label, latest_start

def main():
    with open(GENES_JSON) as f:
        genes_data = json.load(f)

    # Top 20 by rank (already sorted, rank=1 first)
    top20 = sorted(genes_data, key=lambda g: g.get("rank", 9999))[:20]
    gene_names = [g["gene"] for g in top20]

    print(f"Querying ClinicalTrials.gov for {len(gene_names)} genes…\n")

    results = []
    for gene in gene_names:
        data = fetch_trials(gene)
        if data is None:
            entry = {"gene": gene, "n_active_trials": None, "max_phase": None, "latest_start": None}
            print(f"  {gene:<12} → ERROR (null)")
        else:
            n_active, max_phase, latest_start = parse_response(data)
            entry = {
                "gene": gene,
                "n_active_trials": n_active,
                "max_phase": max_phase,
                "latest_start": latest_start,
            }
            phase_str = max_phase or "—"
            print(f"  {gene:<12} → {n_active:>3} active trials  max_phase={phase_str:<8}  latest={latest_start or '—'}")
        results.append(entry)
        time.sleep(0.4)

    with open(TRIALS_JSON, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nWrote {TRIALS_JSON} ({len(results)} genes)")

if __name__ == "__main__":
    main()
