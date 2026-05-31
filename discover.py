"""
discover.py — Data pipeline for RepurposeAgent.
Phase 1: Gold standard from Open Targets
Phase 2: Orphan genes from GWAS Catalog
Phase 3: Enrich orphan genes → writes data.py

Usage:
  python discover.py            # live APIs
  python discover.py --fallback # hardcoded fallback if APIs blocked
"""

import urllib.request
import urllib.parse
import json
import time
import sys
import math
import argparse
from datetime import datetime

# ── Constants ──────────────────────────────────────────────────────────────

DISEASE_SCOPE = [
    {"name": "Alzheimer's disease",        "efo": "MONDO_0004975", "gwas_id": "MONDO_0004975", "area": "neurological"},
    {"name": "Parkinson's disease",        "efo": "MONDO_0005180", "gwas_id": "MONDO_0005180", "area": "neurological"},
    {"name": "coronary artery disease",    "efo": "MONDO_0021661", "gwas_id": "MONDO_0005010", "area": "cardiovascular"},
    {"name": "type 2 diabetes",            "efo": "MONDO_0005148", "gwas_id": "MONDO_0005148", "area": "metabolic"},
    {"name": "schizophrenia",              "efo": "MONDO_0005090", "gwas_id": "MONDO_0005090", "area": "psychiatric"},
    {"name": "inflammatory bowel disease", "efo": "EFO_0003767",   "gwas_id": "MONDO_0005265", "area": "immune"},
]

GOLD_PER_DISEASE     = 3
ORPHAN_MAX_PER_DISEASE = 30

BURDEN_MAP = {
    "neurological":   28.0,
    "cardiovascular": 182.0,
    "metabolic":      42.0,
    "psychiatric":    13.0,
    "immune":         3.0,
}

TISSUE_SPECIFICITY_MAP = {
    "neurological":   0.78,
    "cardiovascular": 0.52,
    "metabolic":      0.44,
    "psychiatric":    0.71,
    "immune":         0.38,
}

OT_GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"

# ── HTTP helpers ────────────────────────────────────────────────────────────

def _post_json(url, payload, timeout=30):
    data = json.dumps(payload).encode("utf-8")
    req  = urllib.request.Request(
        url, data=data,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8"))


def _get_json(url, timeout=30):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8"))


# ── Phase 1: Gold standard ──────────────────────────────────────────────────

GOLD_QUERY = """
query GoldStandard($efoId: String!) {
  disease(efoId: $efoId) {
    associatedTargets(
      page: { index: 0, size: 100 }
    ) {
      rows {
        target {
          approvedSymbol
          id
          drugAndClinicalCandidates {
            rows {
              drug { name }
              maxClinicalStage
            }
          }
        }
        score
        datatypeScores {
          id
          score
        }
      }
    }
  }
}
"""


def _genetics_score(datatype_scores):
    genetic_ids = {"genetic_association", "genetic_literature"}
    best = 0.0
    for ds in datatype_scores:
        if ds.get("id") in genetic_ids:
            best = max(best, ds.get("score", 0.0))
    return best


def phase1_gold_standard():
    print("\n=== Phase 1: Building gold standard from Open Targets ===")
    gold_standard = {}
    gold_gene_set  = set()

    for disease in DISEASE_SCOPE:
        efo = disease["efo"]
        print(f"  Querying {disease['name']} ({efo}) ...", end=" ", flush=True)
        try:
            resp = _post_json(OT_GRAPHQL_URL, {
                "query":     GOLD_QUERY,
                "variables": {"efoId": efo},
            })
            rows = (resp.get("data", {})
                       .get("disease", {})
                       .get("associatedTargets", {})
                       .get("rows", []))
        except Exception as e:
            print(f"ERROR: {e}")
            rows = []

        candidates = []
        for row in rows:
            target = row.get("target", {})
            gene   = target.get("approvedSymbol", "")
            if not gene:
                continue

            gen_score = _genetics_score(row.get("datatypeScores", []))
            if gen_score == 0:
                continue

            drugs = target.get("drugAndClinicalCandidates", {}).get("rows", [])
            approved = [d for d in drugs if d.get("maxClinicalStage") == "APPROVAL"]
            if not approved:
                continue

            best_drug = approved[0].get("drug", {}).get("name", "unknown")

            candidates.append({
                "gene":           gene,
                "disease":        disease["name"],
                "area":           disease["area"],
                "drug":           best_drug,
                "ot_score":       round(row.get("score", 0.0), 4),
                "genetics_score": round(gen_score, 4),
            })

        candidates.sort(key=lambda x: x["genetics_score"], reverse=True)
        added = 0
        for c in candidates:
            if added >= GOLD_PER_DISEASE:
                break
            if c["gene"] not in gold_gene_set:
                gold_standard[c["gene"]] = c
                gold_gene_set.add(c["gene"])
                added += 1

        print(f"found {added} gold genes")
        time.sleep(0.8)

    # Print table
    print(f"\n  {'Disease':<35} {'Gene':<10} {'Drug':<25} {'Gen.Score'}")
    print("  " + "-" * 80)
    for gene, info in gold_standard.items():
        print(f"  {info['disease']:<35} {gene:<10} {info['drug']:<25} {info['genetics_score']}")

    return gold_standard, gold_gene_set


# ── Phase 2: Orphan genes ───────────────────────────────────────────────────

def _parse_gene_from_allele(allele_name):
    """Extract gene symbol from risk allele name like 'GENE-rs...' or 'GENE'."""
    if not allele_name:
        return None
    part = allele_name.split("-")[0].split("_")[0].strip()
    if part and part[0].isalpha() and len(part) >= 2:
        return part.upper()
    return None


def _fetch_gwas_associations(gwas_id, page_size=500):
    """Fetch all genome-wide significant associations for a GWAS catalog trait ID."""
    url = (
        f"https://www.ebi.ac.uk/gwas/rest/api/efoTraits/{gwas_id}/associations"
        f"?size={page_size}&page=0"
    )
    try:
        data = _get_json(url, timeout=60)
        return data.get("_embedded", {}).get("associations", [])
    except Exception as e:
        print(f"ERROR({gwas_id}): {e}")
        return []


def phase2_orphan_genes(gold_gene_set):
    print("\n=== Phase 2: Finding orphan genes from GWAS Catalog ===")
    orphan_genes = {}

    for disease in DISEASE_SCOPE:
        gwas_id = disease["gwas_id"]
        print(f"  Querying {disease['name']} ({gwas_id}) ...", end=" ", flush=True)

        assocs = _fetch_gwas_associations(gwas_id)

        gene_gwas = {}
        for assoc in assocs:
            pval = assoc.get("pvalue", None)
            if pval is None or pval <= 0:
                continue
            # Apply genome-wide significance threshold
            if float(pval) > 5e-8:
                continue
            try:
                pval_log = min(50, max(8, -math.log10(float(pval))))
            except Exception:
                continue

            loci = assoc.get("loci", [])
            for locus in loci:
                # Gene names from authorReportedGenes (more reliable than parsing allele names)
                for gene_entry in locus.get("authorReportedGenes", []):
                    gene = gene_entry.get("geneName", "").strip().upper()
                    if not gene or gene in gold_gene_set:
                        continue
                    if gene in gene_gwas:
                        gene_gwas[gene]["pval_log"] = max(gene_gwas[gene]["pval_log"], pval_log)
                        gene_gwas[gene]["n_studies"] += 1
                    else:
                        gene_gwas[gene] = {"pval_log": pval_log, "n_studies": 1}

        top_genes = sorted(gene_gwas.keys(), key=lambda g: gene_gwas[g]["pval_log"], reverse=True)
        top_genes = top_genes[:ORPHAN_MAX_PER_DISEASE]

        for gene in top_genes:
            if gene not in orphan_genes:
                orphan_genes[gene] = {
                    "gwas_pval_log10": round(gene_gwas[gene]["pval_log"], 2),
                    "n_gwas_studies":  gene_gwas[gene]["n_studies"],
                    "disease_area":    disease["area"],
                    "disease_name":    disease["name"],
                    "burden_daly_m":   BURDEN_MAP[disease["area"]],
                }

        print(f"found {len(top_genes)} orphan genes")
        time.sleep(1.0)

    return orphan_genes


# ── Phase 3: Enrich orphan genes ────────────────────────────────────────────

SEARCH_QUERY = """
query SearchGene($symbol: String!) {
  search(queryString: $symbol, entityNames: ["target"], page: {index: 0, size: 1}) {
    hits { id name }
  }
}
"""

ENRICH_QUERY = """
query EnrichGene($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    approvedSymbol
    tractability {
      label
      modality
      value
    }
    associatedDiseases(page: {index: 0, size: 3}) {
      rows { score }
    }
  }
}
"""


def _fetch_pubmed_count(gene, years=5):
    year_start = datetime.now().year - years
    term = urllib.parse.quote(f"{gene}[Gene Name] AND {year_start}:{datetime.now().year}[PDAT]")
    url  = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={term}&retmode=json&retmax=0"
    try:
        data = _get_json(url, timeout=20)
        return int(data.get("esearchresult", {}).get("count", 0))
    except Exception:
        return 0


def phase3_enrich(orphan_genes):
    print("\n=== Phase 3: Enriching orphan genes with biological features ===")
    total = len(orphan_genes)

    for i, (gene, feat) in enumerate(orphan_genes.items(), 1):
        print(f"  [{i}/{total}] {gene} ...", end=" ", flush=True)

        # Open Targets — step 1: resolve symbol → ensembl ID
        ensembl_id = None
        try:
            s_resp = _post_json(OT_GRAPHQL_URL, {
                "query":     SEARCH_QUERY,
                "variables": {"symbol": gene},
            })
            hits = s_resp.get("data", {}).get("search", {}).get("hits", [])
            if hits and hits[0].get("name", "").upper() == gene.upper():
                ensembl_id = hits[0]["id"]
        except Exception:
            pass

        # Open Targets — step 2: enrich with ensembl ID
        ot_row = {}
        if ensembl_id:
            try:
                e_resp = _post_json(OT_GRAPHQL_URL, {
                    "query":     ENRICH_QUERY,
                    "variables": {"ensemblId": ensembl_id},
                })
                ot_row = e_resp.get("data", {}).get("target", {}) or {}
            except Exception:
                pass

        assoc_rows = ot_row.get("associatedDiseases", {}).get("rows", [])
        ot_score   = max((r.get("score", 0) for r in assoc_rows), default=0.3)
        ot_score   = max(ot_score, 0.3)

        # Tractability: SM values are boolean (Approved Drug > Advanced Clinical > Phase 1)
        tractability = ot_row.get("tractability", [])
        sm_labels = {t.get("label"): t.get("value") for t in tractability if t.get("modality") == "SM"}
        if sm_labels.get("Approved Drug"):
            druggability = 0.90
        elif sm_labels.get("Advanced Clinical"):
            druggability = 0.70
        elif sm_labels.get("Phase 1 Clinical"):
            druggability = 0.50
        else:
            druggability = 0.35
        druggability = max(0.0, min(1.0, druggability))

        time.sleep(0.4)

        # PubMed
        pubmed_count = _fetch_pubmed_count(gene)
        time.sleep(0.35)

        # Proxy features
        gwas_p = feat["gwas_pval_log10"]
        n_stud = feat["n_gwas_studies"]

        eqtl_effect = (gwas_p / 50) * 0.5 + ot_score * 0.5
        eqtl_effect = max(0.1, min(0.95, eqtl_effect))

        mr_z_score  = (gwas_p / 10) * (1 + n_stud / 20)
        mr_z_score  = max(0.5, min(8.0, mr_z_score))

        tissue_spec = TISSUE_SPECIFICITY_MAP[feat["disease_area"]]

        ppi_degree = (
            180 if pubmed_count > 2000 else
            110 if pubmed_count > 800  else
            65  if pubmed_count > 300  else 35
        )

        feat.update({
            "open_targets_score": round(ot_score, 3),
            "druggability_score": round(druggability, 3),
            "eqtl_effect":        round(eqtl_effect, 3),
            "mr_z_score":         round(mr_z_score, 2),
            "tissue_specificity": tissue_spec,
            "ppi_degree":         ppi_degree,
            "pubmed_count_5yr":   pubmed_count,
            "data_source":        "live_api",
        })

        print(f"ot={ot_score:.2f} drug={druggability:.2f} pubmed={pubmed_count}")

    return orphan_genes


# ── Write data.py ────────────────────────────────────────────────────────────

def write_data_py(gold_standard, orphan_genes):
    timestamp = datetime.now().isoformat(timespec="seconds")
    diseases  = list({v["disease"] for v in gold_standard.values()})

    lines = [
        '"""',
        "data.py — auto-generated by discover.py",
        f"Generated: {timestamp}",
        f"Gold standard genes: {len(gold_standard)}",
        f"Orphan genes: {len(orphan_genes)}",
        "DO NOT edit manually.",
        '"""',
        "",
        "DISCOVERY_META = " + json.dumps({
            "generated_at": timestamp,
            "n_gold":       len(gold_standard),
            "n_orphan":     len(orphan_genes),
            "diseases":     diseases,
        }, indent=4),
        "",
        "# Discovered from Open Targets — 2-3 approved drugs per disease",
        "# gene → {disease, area, drug, ot_score, genetics_score}",
        "GOLD_STANDARD = " + json.dumps(gold_standard, indent=4),
        "",
        "GOLD_GENE_SET = set(GOLD_STANDARD.keys())",
        "",
        "# Discovered from GWAS Catalog — no approved drug",
        "# gene → full feature dict",
        "ORPHAN_GENES = " + json.dumps(orphan_genes, indent=4),
        "",
        "FEATURE_KEYS = [",
        '    "gwas_pval_log10", "n_gwas_studies", "tissue_specificity",',
        '    "ppi_degree", "pubmed_count_5yr", "eqtl_effect",',
        '    "druggability_score", "mr_z_score", "open_targets_score"',
        "]",
        "",
    ]

    with open("data.py", "w") as f:
        f.write("\n".join(lines))
    print(f"\ndata.py written ({len(gold_standard)} gold, {len(orphan_genes)} orphan genes)")


# ── Fallback ─────────────────────────────────────────────────────────────────

FALLBACK_GOLD = {
    "HMGCR":  {"gene": "HMGCR",  "disease": "coronary artery disease",    "area": "cardiovascular", "drug": "atorvastatin",   "ot_score": 0.85, "genetics_score": 0.72},
    "PCSK9":  {"gene": "PCSK9",  "disease": "coronary artery disease",    "area": "cardiovascular", "drug": "evolocumab",     "ot_score": 0.82, "genetics_score": 0.68},
    "GLP1R":  {"gene": "GLP1R",  "disease": "type 2 diabetes",            "area": "metabolic",      "drug": "semaglutide",    "ot_score": 0.80, "genetics_score": 0.65},
    "KCNJ11": {"gene": "KCNJ11", "disease": "type 2 diabetes",            "area": "metabolic",      "drug": "glibenclamide",  "ot_score": 0.75, "genetics_score": 0.60},
    "DRD2":   {"gene": "DRD2",   "disease": "schizophrenia",              "area": "psychiatric",    "drug": "haloperidol",    "ot_score": 0.88, "genetics_score": 0.74},
    "HTR2A":  {"gene": "HTR2A",  "disease": "schizophrenia",              "area": "psychiatric",    "drug": "clozapine",      "ot_score": 0.79, "genetics_score": 0.61},
    "TNF":    {"gene": "TNF",    "disease": "inflammatory bowel disease", "area": "immune",         "drug": "infliximab",     "ot_score": 0.91, "genetics_score": 0.80},
    "IL23A":  {"gene": "IL23A",  "disease": "inflammatory bowel disease", "area": "immune",         "drug": "risankizumab",   "ot_score": 0.83, "genetics_score": 0.70},
    "LRRK2":  {"gene": "LRRK2",  "disease": "Parkinson's disease",        "area": "neurological",   "drug": "DNL201",         "ot_score": 0.77, "genetics_score": 0.63},
    "SNCA":   {"gene": "SNCA",   "disease": "Parkinson's disease",        "area": "neurological",   "drug": "prasinezumab",   "ot_score": 0.74, "genetics_score": 0.59},
    "PSEN1":  {"gene": "PSEN1",  "disease": "Alzheimer's disease",        "area": "neurological",   "drug": "semagacestat",   "ot_score": 0.76, "genetics_score": 0.62},
    "APP":    {"gene": "APP",    "disease": "Alzheimer's disease",        "area": "neurological",   "drug": "gantenerumab",   "ot_score": 0.73, "genetics_score": 0.58},
}

FALLBACK_ORPHAN = {
    "CLU":    {"gwas_pval_log10": 18.0, "n_gwas_studies": 8,  "disease_area": "neurological",   "disease_name": "Alzheimer's disease",        "burden_daly_m": 28.0,  "open_targets_score": 0.55, "druggability_score": 0.42, "eqtl_effect": 0.45, "mr_z_score": 2.1, "tissue_specificity": 0.78, "ppi_degree": 110, "pubmed_count_5yr": 950, "data_source": "fallback"},
    "BIN1":   {"gwas_pval_log10": 16.5, "n_gwas_studies": 7,  "disease_area": "neurological",   "disease_name": "Alzheimer's disease",        "burden_daly_m": 28.0,  "open_targets_score": 0.50, "druggability_score": 0.38, "eqtl_effect": 0.41, "mr_z_score": 1.9, "tissue_specificity": 0.78, "ppi_degree": 110, "pubmed_count_5yr": 820, "data_source": "fallback"},
    "PICALM": {"gwas_pval_log10": 15.0, "n_gwas_studies": 6,  "disease_area": "neurological",   "disease_name": "Alzheimer's disease",        "burden_daly_m": 28.0,  "open_targets_score": 0.48, "druggability_score": 0.35, "eqtl_effect": 0.39, "mr_z_score": 1.8, "tissue_specificity": 0.78, "ppi_degree": 65,  "pubmed_count_5yr": 600, "data_source": "fallback"},
    "GBA":    {"gwas_pval_log10": 14.0, "n_gwas_studies": 5,  "disease_area": "neurological",   "disease_name": "Parkinson's disease",        "burden_daly_m": 28.0,  "open_targets_score": 0.62, "druggability_score": 0.70, "eqtl_effect": 0.44, "mr_z_score": 2.0, "tissue_specificity": 0.78, "ppi_degree": 110, "pubmed_count_5yr": 1200, "data_source": "fallback"},
    "MAPT":   {"gwas_pval_log10": 20.0, "n_gwas_studies": 10, "disease_area": "neurological",   "disease_name": "Parkinson's disease",        "burden_daly_m": 28.0,  "open_targets_score": 0.65, "druggability_score": 0.45, "eqtl_effect": 0.55, "mr_z_score": 2.8, "tissue_specificity": 0.78, "ppi_degree": 180, "pubmed_count_5yr": 2500, "data_source": "fallback"},
    "LPA":    {"gwas_pval_log10": 22.0, "n_gwas_studies": 12, "disease_area": "cardiovascular", "disease_name": "coronary artery disease",    "burden_daly_m": 182.0, "open_targets_score": 0.72, "druggability_score": 0.60, "eqtl_effect": 0.56, "mr_z_score": 3.2, "tissue_specificity": 0.52, "ppi_degree": 110, "pubmed_count_5yr": 1800, "data_source": "fallback"},
    "LDLR":   {"gwas_pval_log10": 19.0, "n_gwas_studies": 9,  "disease_area": "cardiovascular", "disease_name": "coronary artery disease",    "burden_daly_m": 182.0, "open_targets_score": 0.70, "druggability_score": 0.65, "eqtl_effect": 0.52, "mr_z_score": 2.8, "tissue_specificity": 0.52, "ppi_degree": 180, "pubmed_count_5yr": 2100, "data_source": "fallback"},
    "APOB":   {"gwas_pval_log10": 17.5, "n_gwas_studies": 8,  "disease_area": "cardiovascular", "disease_name": "coronary artery disease",    "burden_daly_m": 182.0, "open_targets_score": 0.68, "druggability_score": 0.55, "eqtl_effect": 0.49, "mr_z_score": 2.5, "tissue_specificity": 0.52, "ppi_degree": 110, "pubmed_count_5yr": 1500, "data_source": "fallback"},
    "TCF7L2": {"gwas_pval_log10": 25.0, "n_gwas_studies": 15, "disease_area": "metabolic",      "disease_name": "type 2 diabetes",            "burden_daly_m": 42.0,  "open_targets_score": 0.75, "druggability_score": 0.50, "eqtl_effect": 0.60, "mr_z_score": 3.5, "tissue_specificity": 0.44, "ppi_degree": 110, "pubmed_count_5yr": 2000, "data_source": "fallback"},
    "PPARG":  {"gwas_pval_log10": 16.0, "n_gwas_studies": 7,  "disease_area": "metabolic",      "disease_name": "type 2 diabetes",            "burden_daly_m": 42.0,  "open_targets_score": 0.68, "druggability_score": 0.75, "eqtl_effect": 0.46, "mr_z_score": 2.2, "tissue_specificity": 0.44, "ppi_degree": 180, "pubmed_count_5yr": 2800, "data_source": "fallback"},
    "CACNA1C":{"gwas_pval_log10": 15.0, "n_gwas_studies": 6,  "disease_area": "psychiatric",    "disease_name": "schizophrenia",              "burden_daly_m": 13.0,  "open_targets_score": 0.55, "druggability_score": 0.70, "eqtl_effect": 0.42, "mr_z_score": 2.0, "tissue_specificity": 0.71, "ppi_degree": 110, "pubmed_count_5yr": 1100, "data_source": "fallback"},
    "NRXN1":  {"gwas_pval_log10": 14.5, "n_gwas_studies": 5,  "disease_area": "psychiatric",    "disease_name": "schizophrenia",              "burden_daly_m": 13.0,  "open_targets_score": 0.50, "druggability_score": 0.40, "eqtl_effect": 0.40, "mr_z_score": 1.8, "tissue_specificity": 0.71, "ppi_degree": 110, "pubmed_count_5yr": 900, "data_source": "fallback"},
    "NOD2":   {"gwas_pval_log10": 18.0, "n_gwas_studies": 9,  "disease_area": "immune",         "disease_name": "inflammatory bowel disease", "burden_daly_m": 3.0,   "open_targets_score": 0.70, "druggability_score": 0.55, "eqtl_effect": 0.50, "mr_z_score": 2.6, "tissue_specificity": 0.38, "ppi_degree": 110, "pubmed_count_5yr": 1600, "data_source": "fallback"},
    "IL10":   {"gwas_pval_log10": 16.0, "n_gwas_studies": 7,  "disease_area": "immune",         "disease_name": "inflammatory bowel disease", "burden_daly_m": 3.0,   "open_targets_score": 0.65, "druggability_score": 0.60, "eqtl_effect": 0.46, "mr_z_score": 2.2, "tissue_specificity": 0.38, "ppi_degree": 110, "pubmed_count_5yr": 1400, "data_source": "fallback"},
    "ATG16L1":{"gwas_pval_log10": 15.5, "n_gwas_studies": 6,  "disease_area": "immune",         "disease_name": "inflammatory bowel disease", "burden_daly_m": 3.0,   "open_targets_score": 0.58, "druggability_score": 0.42, "eqtl_effect": 0.43, "mr_z_score": 2.0, "tissue_specificity": 0.38, "ppi_degree": 65,  "pubmed_count_5yr": 700, "data_source": "fallback"},
    "ABCA7":  {"gwas_pval_log10": 14.0, "n_gwas_studies": 5,  "disease_area": "neurological",   "disease_name": "Alzheimer's disease",        "burden_daly_m": 28.0,  "open_targets_score": 0.46, "druggability_score": 0.55, "eqtl_effect": 0.38, "mr_z_score": 1.7, "tissue_specificity": 0.78, "ppi_degree": 65,  "pubmed_count_5yr": 450, "data_source": "fallback"},
    "CR1":    {"gwas_pval_log10": 13.5, "n_gwas_studies": 4,  "disease_area": "neurological",   "disease_name": "Alzheimer's disease",        "burden_daly_m": 28.0,  "open_targets_score": 0.44, "druggability_score": 0.38, "eqtl_effect": 0.36, "mr_z_score": 1.6, "tissue_specificity": 0.78, "ppi_degree": 35,  "pubmed_count_5yr": 350, "data_source": "fallback"},
    "PINK1":  {"gwas_pval_log10": 13.0, "n_gwas_studies": 4,  "disease_area": "neurological",   "disease_name": "Parkinson's disease",        "burden_daly_m": 28.0,  "open_targets_score": 0.58, "druggability_score": 0.50, "eqtl_effect": 0.40, "mr_z_score": 1.8, "tissue_specificity": 0.78, "ppi_degree": 65,  "pubmed_count_5yr": 800, "data_source": "fallback"},
    "CDKAL1": {"gwas_pval_log10": 14.5, "n_gwas_studies": 6,  "disease_area": "metabolic",      "disease_name": "type 2 diabetes",            "burden_daly_m": 42.0,  "open_targets_score": 0.45, "druggability_score": 0.35, "eqtl_effect": 0.38, "mr_z_score": 1.7, "tissue_specificity": 0.44, "ppi_degree": 35,  "pubmed_count_5yr": 280, "data_source": "fallback"},
    "DISC1":  {"gwas_pval_log10": 12.5, "n_gwas_studies": 3,  "disease_area": "psychiatric",    "disease_name": "schizophrenia",              "burden_daly_m": 13.0,  "open_targets_score": 0.40, "druggability_score": 0.30, "eqtl_effect": 0.34, "mr_z_score": 1.4, "tissue_specificity": 0.71, "ppi_degree": 35,  "pubmed_count_5yr": 200, "data_source": "fallback"},
    "IL6":    {"gwas_pval_log10": 17.0, "n_gwas_studies": 8,  "disease_area": "immune",         "disease_name": "inflammatory bowel disease", "burden_daly_m": 3.0,   "open_targets_score": 0.72, "druggability_score": 0.75, "eqtl_effect": 0.50, "mr_z_score": 2.5, "tissue_specificity": 0.38, "ppi_degree": 180, "pubmed_count_5yr": 3500, "data_source": "fallback"},
}


def run_fallback():
    print("\n=== FALLBACK MODE — using hardcoded gene-drug pairs ===")
    gold_standard = FALLBACK_GOLD
    orphan_genes  = FALLBACK_ORPHAN
    return gold_standard, orphan_genes


# ── Summary ──────────────────────────────────────────────────────────────────

def print_summary(gold_standard, orphan_genes):
    print("\n=== Discovery Complete ===")
    print(f"Gold standard: {len(gold_standard)} genes (2-3 per disease)")
    by_disease = {}
    for gene, info in gold_standard.items():
        by_disease.setdefault(info["disease"], []).append((gene, info["drug"]))
    for disease, pairs in by_disease.items():
        pairs_str = ", ".join(f"{g} ({d})" for g, d in pairs)
        print(f"  {disease}: {pairs_str}")

    print(f"\nOrphan genes: {len(orphan_genes)} total")
    by_area = {}
    for gene, info in orphan_genes.items():
        by_area.setdefault(info["disease_area"], []).append(gene)
    for area, genes in sorted(by_area.items()):
        print(f"  {area}: {len(genes)} genes")

    print("\ndata.py written. Next: python agent.py")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fallback", action="store_true",
                        help="Use hardcoded fallback data instead of live APIs")
    args = parser.parse_args()

    if args.fallback:
        gold_standard, orphan_genes = run_fallback()
    else:
        gold_standard, gold_gene_set = phase1_gold_standard()
        if not gold_standard:
            print("WARNING: No gold standard genes found — falling back to hardcoded set")
            gold_standard = FALLBACK_GOLD
            gold_gene_set = set(FALLBACK_GOLD.keys())

        orphan_genes = phase2_orphan_genes(gold_gene_set)
        if not orphan_genes:
            print("WARNING: No orphan genes found — falling back to hardcoded set")
            orphan_genes = FALLBACK_ORPHAN

        orphan_genes = phase3_enrich(orphan_genes)

    write_data_py(gold_standard, orphan_genes)
    print_summary(gold_standard, orphan_genes)


if __name__ == "__main__":
    main()
