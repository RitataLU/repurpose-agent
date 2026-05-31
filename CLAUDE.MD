# CLAUDE.md — RepurposeAgent
## Fully autonomous. No hardcoded gene lists. No hardcoded gold standard.
## Read every section before writing any code.

---

## What you are building

A drug repurposing research agent where **everything is discovered from APIs**:

1. `discover.py` — one script that does three things in sequence:
   - **Phase 1:** query Open Targets to find 2–3 approved drugs per disease → gold standard
   - **Phase 2:** query GWAS Catalog to find orphan genes per disease → candidate pool
   - **Phase 3:** enrich each orphan gene with biological features → writes `data.py`

2. `scorer.py` — the only file the agent edits (baseline equal-weight formula)

3. `evaluate.py` — fixed harness measuring recall@K against the discovered gold standard

4. `agent.py` — autoresearch loop: modify scorer → score genes → measure recall → keep/discard → git push

5. `docs/index.html` + `docs/showcase.html` — GitHub Pages site

**GitHub remote:** `https://github.com/RitataLU/repurpose-agent.git`
**GitHub Pages:** `https://ritataLU.github.io/repurpose-agent/`

---

## Disease scope — the only hardcoded thing

These 6 disease EFO IDs are the only values you hardcode.
Everything else — gold standard, orphan genes, features — comes from APIs.

```python
DISEASE_SCOPE = [
    {"name": "Alzheimer's disease",       "efo": "EFO_0000249", "area": "neurological"},
    {"name": "Parkinson's disease",       "efo": "EFO_0002508", "area": "neurological"},
    {"name": "coronary artery disease",   "efo": "EFO_0001645", "area": "cardiovascular"},
    {"name": "type 2 diabetes",           "efo": "EFO_0001360", "area": "metabolic"},
    {"name": "schizophrenia",             "efo": "EFO_0000692", "area": "psychiatric"},
    {"name": "inflammatory bowel disease","efo": "EFO_0003767", "area": "immune"},
]
```

---

## Execution order — do not deviate

```
Step 1  → scaffold repo + git init + first push
Step 2  → discover.py      (Phase 1 + 2 + 3 — the whole data pipeline)
Step 3  → scorer.py        (baseline formula)
Step 4  → evaluate.py      (fixed recall harness)
Step 5  → agent.py         (loop + git auto-push)
Step 6  → docs/index.html  (live dashboard)
Step 7  → docs/showcase.html (public page)
Step 8  → README.md
Step 9  → run: python discover.py
Step 10 → run: python agent.py
Step 11 → enable GitHub Pages
```

---

## Step 1 — Scaffold

```bash
mkdir -p repurpose_agent/docs
cd repurpose_agent
git init
git checkout -b main
git remote add origin https://github.com/RitataLU/repurpose-agent.git

cat > .gitignore << 'EOF'
__pycache__/
*.pyc
*.bak
scorer.py.bak
EOF

echo '[]' > docs/results.json
echo '[]' > docs/genes.json

git add .
git commit -m "init: RepurposeAgent scaffold"
git push -u origin main
```

---

## Step 2 — `discover.py`

This is the heart of the project. One script, three phases, no hardcoded biology.
Use only `urllib`, `json`, `time`, `sys`. No pip installs.

### Phase 1: Build gold standard from Open Targets

For each disease, query Open Targets for genes that have:
- An approved drug (`maximumClinicalTrialPhase == 4`)
- A known GWAS association for that disease (evidenceCount > 0 from genetics)

Take the **top 2–3 genes per disease** ranked by Open Targets genetics evidence score.
These become `GOLD_STANDARD` — the evaluation ground truth.

```python
DISEASE_SCOPE = [
    {"name": "Alzheimer's disease",       "efo": "EFO_0000249", "area": "neurological"},
    {"name": "Parkinson's disease",       "efo": "EFO_0002508", "area": "neurological"},
    {"name": "coronary artery disease",   "efo": "EFO_0001645", "area": "cardiovascular"},
    {"name": "type 2 diabetes",           "efo": "EFO_0001360", "area": "metabolic"},
    {"name": "schizophrenia",             "efo": "EFO_0000692", "area": "psychiatric"},
    {"name": "inflammatory bowel disease","efo": "EFO_0003767", "area": "immune"},
]

GOLD_PER_DISEASE = 3   # max gold standard genes to take per disease
ORPHAN_MAX_PER_DISEASE = 30  # cap orphan genes per disease for API efficiency
```

**Open Targets GraphQL query for Phase 1:**

```graphql
query GoldStandard($efoId: String!) {
  disease(efoId: $efoId) {
    associatedTargets(
      page: { index: 0, size: 50 }
      enableIndirect: false
    ) {
      rows {
        target {
          approvedSymbol
          id
        }
        score
        datatypeScores {
          id       # look for "genetics_literature" or "gwas"
          score
        }
        knownDrugs: target {
          knownDrugs(size: 5) {
            rows {
              maximumClinicalTrialPhase
              drug { name approvedYear }
            }
          }
        }
      }
    }
  }
}
```

**Phase 1 logic in pseudocode — implement this exactly:**

```
gold_standard = {}   # gene_symbol → {disease, drug_name, area, ot_score}
gold_gene_set = set()

for each disease in DISEASE_SCOPE:
    query Open Targets with disease EFO id
    get rows sorted by score descending
    
    candidates = []
    for each row in rows:
        gene = row.target.approvedSymbol
        
        # Must have GWAS/genetics evidence
        genetics_score = find datatypeScore where id in
            ["gwas", "genetics_literature", "affected_pathway"]
        if genetics_score == 0 or missing: skip
        
        # Must have approved drug (phase 4)
        drugs = row.target.knownDrugs.rows
        approved = [d for d in drugs if d.maximumClinicalTrialPhase >= 4]
        if not approved: skip
        
        # Take best approved drug (highest phase, then alphabetical)
        best_drug = approved[0].drug.name
        
        candidates.append({
            "gene": gene,
            "disease": disease.name,
            "area": disease.area,
            "drug": best_drug,
            "ot_score": row.score,
            "genetics_score": genetics_score,
        })
    
    # Take top GOLD_PER_DISEASE by genetics_score
    candidates.sort(by genetics_score, descending)
    for c in candidates[:GOLD_PER_DISEASE]:
        gold_standard[c.gene] = c
        gold_gene_set.add(c.gene)

    sleep(0.8)

print gold standard table:
    disease | gene | drug | genetics_score
```

### Phase 2: Find orphan genes from GWAS Catalog

For each disease, query GWAS Catalog for genome-wide significant associations (p < 5e-8).
Extract gene symbols. Remove any gene in `gold_gene_set`. These are the orphan candidates.

**GWAS Catalog REST endpoint:**
```
GET https://www.ebi.ac.uk/gwas/rest/api/associations/search
    ?q=efoTrait:{efo_id}
    &pvalueUpperLimit=0.00000005
    &size=500&page=0
```

**Phase 2 logic:**
```
orphan_genes = {}   # gene_symbol → feature dict

for each disease in DISEASE_SCOPE:
    fetch GWAS associations for disease EFO id
    
    gene_gwas = {}  # gene → {gwas_pval_log10, n_studies}
    
    for each association:
        extract gene symbols from loci → strongestRiskAlleles
        parse gene from riskAlleleName (format: "GENE-rs..." or "GENE")
        
        pval_log = -log10(association.pvalue), capped 8–50
        
        for each gene found:
            if gene in gold_gene_set: skip (already a gold standard)
            if gene already seen: update max pval, increment n_studies
            else: add with pval and n_studies=1
    
    # Take top ORPHAN_MAX_PER_DISEASE by gwas_pval_log10
    top_genes = sorted(gene_gwas, by pval_log, descending)[:ORPHAN_MAX_PER_DISEASE]
    
    for gene in top_genes:
        orphan_genes[gene] = {
            "gwas_pval_log10": gene_gwas[gene]["gwas_pval_log10"],
            "n_gwas_studies":  gene_gwas[gene]["n_studies"],
            "disease_area":    disease.area,
            "disease_name":    disease.name,
            "burden_daly_m":   BURDEN_MAP[disease.area],
        }
    
    sleep(1.0)
```

**BURDEN_MAP** (static reference values, DALYs in millions, GBD 2021):
```python
BURDEN_MAP = {
    "neurological":   28.0,
    "cardiovascular": 182.0,
    "metabolic":      42.0,
    "psychiatric":    13.0,
    "immune":         3.0,
}
```

### Phase 3: Enrich orphan genes with biological features

For each orphan gene, fetch additional features from Open Targets and PubMed.
Then compute proxy features for eQTL and MR (no GTEx/MR-Base call needed).

**Open Targets enrichment query per gene:**
```graphql
query EnrichGene($symbol: String!) {
  targets(queryString: $symbol, page: {index: 0, size: 1}) {
    rows {
      id
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
}
```

**Phase 3 logic:**
```
for each gene in orphan_genes:
    # Open Targets
    ot_data = query Open Targets for gene symbol
    
    ot_score = max association score across top 3 diseases, default 0.3
    
    tractability = ot_data.tractability
    sm_vals = [t.value for t in tractability if t.modality == "SM"]
    druggability_score = max(sm_vals) if sm_vals else 0.35
    clamp to [0.0, 1.0]
    
    sleep(0.4)
    
    # PubMed count (last 5 years)
    pubmed_count = fetch_pubmed_count(gene)  # eutils esearch
    sleep(0.35)
    
    # Proxy features (principled approximations)
    eqtl_effect = (gwas_pval_log10/50)*0.5 + ot_score*0.5
                  clamped [0.1, 0.95]
    
    mr_z_score = (gwas_pval_log10/10) * (1 + n_gwas_studies/20)
                 clamped [0.5, 8.0]
    
    tissue_specificity = {
        "neurological": 0.78, "cardiovascular": 0.52,
        "metabolic": 0.44,    "psychiatric": 0.71, "immune": 0.38
    }[disease_area]
    
    ppi_degree = (
        180 if pubmed_count > 2000 else
        110 if pubmed_count > 800  else
        65  if pubmed_count > 300  else 35
    )
    
    # Update orphan_genes entry with full feature set
    orphan_genes[gene].update({
        "open_targets_score": round(ot_score, 3),
        "druggability_score": round(druggability_score, 3),
        "eqtl_effect":        round(eqtl_effect, 3),
        "mr_z_score":         round(mr_z_score, 2),
        "tissue_specificity": tissue_specificity,
        "ppi_degree":         ppi_degree,
        "pubmed_count_5yr":   pubmed_count,
        "data_source":        "live_api",
    })
    
    print: gene | ot=X | drug=X | pubmed=X
```

### Write `data.py`

After all three phases, write `data.py` with this exact structure:

```python
"""
data.py — auto-generated by discover.py
Generated: {timestamp}
Gold standard genes: {n}
Orphan genes: {n}
DO NOT edit manually.
"""

DISCOVERY_META = {
    "generated_at": "...",
    "n_gold":       N,
    "n_orphan":     N,
    "diseases":     [...],
}

# Discovered from Open Targets — 2-3 approved drugs per disease
# gene → {disease, area, drug, ot_score, genetics_score}
GOLD_STANDARD = { ... }

GOLD_GENE_SET = set(GOLD_STANDARD.keys())

# Discovered from GWAS Catalog — no approved drug
# gene → full feature dict
ORPHAN_GENES = { ... }

FEATURE_KEYS = [
    "gwas_pval_log10", "n_gwas_studies", "tissue_specificity",
    "ppi_degree", "pubmed_count_5yr", "eqtl_effect",
    "druggability_score", "mr_z_score", "open_targets_score"
]
```

### Fallback mode (`--fallback` flag)

If APIs return nothing (network blocked, rate limited), write `data.py` with:
- `GOLD_STANDARD`: 12 real gene-drug pairs (2 per disease), hardcoded as fallback
- `ORPHAN_GENES`: 20 real orphan genes, hardcoded as fallback
- `data_source: "fallback"` on every entry

Use these as fallback gold standard (2 per disease, real pairs):
```
Alzheimer's:  BACE1→verubecestat(failed but Phase3 — use APOE4-targeting drug), 
              use: TNF→"thalidomide AD trial" — actually use real ones:
              Better fallback gold pairs:
              - "HMGCR" → "atorvastatin" (CAD)
              - "PCSK9"  → "evolocumab"  (CAD)
              - "GLP1R"  → "semaglutide" (T2D)
              - "KCNJ11" → "glibenclamide"(T2D)
              - "DRD2"   → "haloperidol" (schizophrenia)
              - "HTR2A"  → "clozapine"   (schizophrenia)
              - "TNF"    → "infliximab"  (IBD)
              - "IL23A"  → "risankizumab"(IBD)
              - "LRRK2"  → "DNL201"      (Parkinson's — Phase 2, use as example)
              - "SNCA"   → "prasinezumab"(Parkinson's — Phase 2)
              - "APOE"   → "lecanemab-related" — skip, use:
              - "PSEN1"  → "semagacestat"(Alzheimer's — Phase 3 trial)
              - "APP"    → "gantenerumab"(Alzheimer's — Phase 3)
```

Note: for fallback, use the most real-world accurate pairs you can find.
The fallback only runs if `python discover.py --fallback` is explicitly called.

### Print summary at the end of discover.py

```
=== Discovery Complete ===
Gold standard: N genes (2-3 per disease)
  Alzheimer's:  GENE1 (drug1), GENE2 (drug2)
  Parkinson's:  ...
  ...

Orphan genes: N total
  neurological: N genes
  cardiovascular: N genes
  ...

data.py written. Next: python agent.py
```

**Test:** `python discover.py`
Expected: data.py written with both GOLD_STANDARD and ORPHAN_GENES populated.
If APIs blocked: `python discover.py --fallback`

---

## Step 3 — `scorer.py`

```python
"""
scorer.py — Gene repurposing scorer.
THE ONLY FILE THE AGENT EDITS.
Baseline: equal weights.
"""

WEIGHTS = {
    "gwas_pval_log10":    1.0,
    "n_gwas_studies":     1.0,
    "tissue_specificity": 1.0,
    "ppi_degree":         1.0,
    "pubmed_count_5yr":   1.0,
    "eqtl_effect":        1.0,
    "druggability_score": 1.0,
    "mr_z_score":         1.0,
    "open_targets_score": 1.0,
}

BURDEN_MULTIPLIER = False
INTERACTION_TERMS = []  # list of (feature_a, feature_b, weight)

def score_gene(features: dict) -> float:
    score = sum(
        WEIGHTS.get(k, 0) * v
        for k, v in features.items()
        if k in WEIGHTS and isinstance(v, (int, float))
    )
    if BURDEN_MULTIPLIER:
        score *= (1 + features.get("burden_daly_m", 0) / 200)
    for a, b, w in INTERACTION_TERMS:
        score += w * features.get(a, 0) * features.get(b, 0)
    return score
```

---

## Step 4 — `evaluate.py`

```python
"""
evaluate.py — Fixed evaluation harness. Never modified.
Metric: recall@K where K = min(20, len(GOLD_GENE_SET))
"""
import importlib

def normalize(genes: dict) -> dict:
    skip = {"disease_area", "disease_name", "data_source", "ensembl_id"}
    sample = next(iter(genes.values()))
    keys = [k for k in sample if k not in skip
            and isinstance(sample.get(k), (int, float))]
    normed = {g: dict(f) for g, f in genes.items()}
    for k in keys:
        vals = [genes[g][k] for g in genes
                if isinstance(genes[g].get(k), (int, float))]
        mn, mx = min(vals), max(vals)
        rng = mx - mn if mx != mn else 1.0
        for g in normed:
            if isinstance(normed[g].get(k), (int, float)):
                normed[g][k] = (normed[g][k] - mn) / rng
    return normed

def run_evaluation() -> dict:
    import data as d
    import scorer as sc
    importlib.reload(d)
    importlib.reload(sc)

    K = min(20, len(d.GOLD_GENE_SET))
    normed = normalize(d.ORPHAN_GENES)
    scores = {g: sc.score_gene(f) for g, f in normed.items()}
    ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    top_K  = [g for g, _ in ranked[:K]]
    hits   = [g for g in top_K if g in d.GOLD_GENE_SET]

    return {
        "recall_at_k":       round(len(hits) / K, 4),
        "K":                 K,
        "n_hits":            len(hits),
        "n_gold":            len(d.GOLD_GENE_SET),
        "n_orphan":          len(d.ORPHAN_GENES),
        "top_k_genes":       top_K,
        "hit_genes":         hits,
        "all_scores":        dict(ranked),
        "weights_snapshot":  dict(sc.WEIGHTS),
        "interaction_snapshot": list(sc.INTERACTION_TERMS),
        "burden_mul":        sc.BURDEN_MULTIPLIER,
    }
```

**Test after discover.py has run:**
```bash
python -c "
from evaluate import run_evaluation
r = run_evaluation()
print(f'Orphan genes: {r[\"n_orphan\"]}')
print(f'Gold standard: {r[\"n_gold\"]} genes (K={r[\"K\"]})')
print(f'Baseline recall@{r[\"K\"]}: {r[\"recall_at_k\"]}')
print(f'Top 3 genes: {r[\"top_k_genes\"][:3]}')
"
```

---

## Step 5 — `agent.py`

```python
"""
agent.py — Autoresearch loop with auto git push.
Usage:
  python agent.py              # full 10-experiment run with git push
  python agent.py --no-push    # local only
  python agent.py --max 3      # limit experiments
"""
import os, sys, shutil, csv, json, subprocess, argparse
from datetime import datetime

RESULTS_TSV   = "results.tsv"
TOP_GENES_TSV = "top_genes.tsv"
DOCS_RESULTS  = "docs/results.json"
DOCS_GENES    = "docs/genes.json"
SCORER_FILE   = "scorer.py"
SCORER_BAK    = "scorer.py.bak"

EXPERIMENTS = [
    {
        "desc": "upweight causal features: MR + druggability",
        "weights": {
            "mr_z_score": 3.0, "druggability_score": 3.0,
            "open_targets_score": 2.0, "gwas_pval_log10": 1.5,
            "n_gwas_studies": 1.0, "tissue_specificity": 1.0,
            "ppi_degree": 0.5, "pubmed_count_5yr": 0.3, "eqtl_effect": 1.5,
        },
        "interaction": [], "burden": False,
    },
    {
        "desc": "MR × druggability interaction term",
        "weights": {
            "mr_z_score": 2.5, "druggability_score": 2.5,
            "open_targets_score": 2.0, "gwas_pval_log10": 1.5,
            "n_gwas_studies": 1.0, "tissue_specificity": 1.0,
            "ppi_degree": 0.5, "pubmed_count_5yr": 0.3, "eqtl_effect": 1.5,
        },
        "interaction": [("mr_z_score", "druggability_score", 4.0)],
        "burden": False,
    },
    {
        "desc": "zero pubmed — remove popularity bias",
        "weights": {
            "mr_z_score": 2.5, "druggability_score": 2.5,
            "open_targets_score": 2.0, "gwas_pval_log10": 1.5,
            "n_gwas_studies": 1.0, "tissue_specificity": 1.2,
            "ppi_degree": 0.5, "pubmed_count_5yr": 0.0, "eqtl_effect": 1.8,
        },
        "interaction": [("mr_z_score", "druggability_score", 4.0)],
        "burden": False,
    },
    {
        "desc": "add disease burden multiplier",
        "weights": {
            "mr_z_score": 2.5, "druggability_score": 2.5,
            "open_targets_score": 2.0, "gwas_pval_log10": 1.5,
            "n_gwas_studies": 1.0, "tissue_specificity": 1.2,
            "ppi_degree": 0.5, "pubmed_count_5yr": 0.0, "eqtl_effect": 1.8,
        },
        "interaction": [("mr_z_score", "druggability_score", 4.0)],
        "burden": True,
    },
    {
        "desc": "eQTL × tissue_specificity interaction",
        "weights": {
            "mr_z_score": 2.5, "druggability_score": 2.5,
            "open_targets_score": 2.0, "gwas_pval_log10": 1.5,
            "n_gwas_studies": 1.2, "tissue_specificity": 1.0,
            "ppi_degree": 0.5, "pubmed_count_5yr": 0.0, "eqtl_effect": 1.0,
        },
        "interaction": [
            ("mr_z_score", "druggability_score", 4.0),
            ("eqtl_effect", "tissue_specificity", 2.5),
        ],
        "burden": True,
    },
    {
        "desc": "zero PPI — network degree is noise",
        "weights": {
            "mr_z_score": 3.0, "druggability_score": 3.0,
            "open_targets_score": 2.5, "gwas_pval_log10": 1.0,
            "n_gwas_studies": 2.5, "tissue_specificity": 1.0,
            "ppi_degree": 0.0, "pubmed_count_5yr": 0.0, "eqtl_effect": 1.2,
        },
        "interaction": [
            ("mr_z_score", "druggability_score", 5.0),
            ("eqtl_effect", "tissue_specificity", 2.0),
        ],
        "burden": True,
    },
    {
        "desc": "dominant MR — pure causal genetics",
        "weights": {
            "mr_z_score": 6.0, "druggability_score": 2.0,
            "open_targets_score": 1.5, "gwas_pval_log10": 0.5,
            "n_gwas_studies": 1.5, "tissue_specificity": 0.8,
            "ppi_degree": 0.0, "pubmed_count_5yr": 0.0, "eqtl_effect": 1.0,
        },
        "interaction": [("mr_z_score", "druggability_score", 3.0)],
        "burden": True,
    },
    {
        "desc": "causal triad: MR + eQTL + OT",
        "weights": {
            "mr_z_score": 3.0, "druggability_score": 2.5,
            "open_targets_score": 3.0, "gwas_pval_log10": 0.8,
            "n_gwas_studies": 2.0, "tissue_specificity": 1.2,
            "ppi_degree": 0.0, "pubmed_count_5yr": 0.0, "eqtl_effect": 2.5,
        },
        "interaction": [
            ("mr_z_score", "druggability_score", 3.5),
            ("eqtl_effect", "open_targets_score", 2.0),
        ],
        "burden": True,
    },
    {
        "desc": "boost GWAS replication signal",
        "weights": {
            "mr_z_score": 3.0, "druggability_score": 2.5,
            "open_targets_score": 3.5, "gwas_pval_log10": 0.8,
            "n_gwas_studies": 3.0, "tissue_specificity": 1.2,
            "ppi_degree": 0.0, "pubmed_count_5yr": 0.0, "eqtl_effect": 2.5,
        },
        "interaction": [
            ("mr_z_score", "druggability_score", 3.5),
            ("eqtl_effect", "open_targets_score", 2.0),
        ],
        "burden": True,
    },
    {
        "desc": "fine-tune: OT score boost",
        "weights": {
            "mr_z_score": 3.0, "druggability_score": 2.5,
            "open_targets_score": 4.0, "gwas_pval_log10": 0.8,
            "n_gwas_studies": 3.0, "tissue_specificity": 1.2,
            "ppi_degree": 0.0, "pubmed_count_5yr": 0.0, "eqtl_effect": 2.5,
        },
        "interaction": [
            ("mr_z_score", "druggability_score", 3.5),
            ("eqtl_effect", "open_targets_score", 2.5),
        ],
        "burden": True,
    },
]

def write_scorer(exp: dict):
    w = json.dumps(exp["weights"], indent=4)
    code = f'''"""scorer.py — {exp["desc"]} — auto-written by agent.py"""

WEIGHTS = {w}
BURDEN_MULTIPLIER = {repr(exp.get("burden", False))}
INTERACTION_TERMS = {repr(exp.get("interaction", []))}

def score_gene(features: dict) -> float:
    score = sum(WEIGHTS.get(k,0)*v for k,v in features.items()
                if k in WEIGHTS and isinstance(v,(int,float)))
    if BURDEN_MULTIPLIER:
        score *= (1 + features.get("burden_daly_m",0)/200)
    for a,b,w in INTERACTION_TERMS:
        score += w * features.get(a,0) * features.get(b,0)
    return score
'''
    with open(SCORER_FILE, "w") as f:
        f.write(code)

def git_push(exp_num: int, recall: float, no_push: bool):
    files = [RESULTS_TSV, TOP_GENES_TSV, DOCS_RESULTS, DOCS_GENES, SCORER_FILE]
    try:
        subprocess.run(["git","add"] + files, check=True, capture_output=True)
        msg = f"exp {exp_num}: recall@K={recall:.4f} [{datetime.now().strftime('%Y-%m-%d %H:%M')}]"
        r = subprocess.run(["git","commit","-m",msg], capture_output=True, text=True)
        if "nothing to commit" in r.stdout:
            return
        if not no_push:
            subprocess.run(["git","push","origin","main"], check=True, capture_output=True)
            print("  → pushed to GitHub")
    except subprocess.CalledProcessError as e:
        print(f"  → git skip: {str(e)[:60]}")

def save_outputs(rows: list, best: dict):
    import data as d
    os.makedirs("docs", exist_ok=True)

    with open(DOCS_RESULTS,"w") as f:
        json.dump(rows, f, indent=2, default=str)

    genes_out = [
        {
            "rank": i,
            "gene": g,
            "score": round(s, 4),
            "hit":  g in best.get("hit_genes", []),
            "disease_area": d.ORPHAN_GENES.get(g,{}).get("disease_area",""),
            "disease_name": d.ORPHAN_GENES.get(g,{}).get("disease_name",""),
            "druggability": round(d.ORPHAN_GENES.get(g,{}).get("druggability_score",0),2),
            "mr_z":    round(d.ORPHAN_GENES.get(g,{}).get("mr_z_score",0),2),
            "ot_score":round(d.ORPHAN_GENES.get(g,{}).get("open_targets_score",0),2),
            "gwas_p":  round(d.ORPHAN_GENES.get(g,{}).get("gwas_pval_log10",0),1),
            "pubmed":  d.ORPHAN_GENES.get(g,{}).get("pubmed_count_5yr",0),
            "source":  d.ORPHAN_GENES.get(g,{}).get("data_source",""),
        }
        for i,(g,s) in enumerate(best["all_scores"].items(),1)
    ]
    with open(DOCS_GENES,"w") as f:
        json.dump(genes_out, f, indent=2)

    # also write gold standard info for the dashboard
    gold_out = [
        {"gene": g, **v}
        for g,v in d.GOLD_STANDARD.items()
    ]
    with open("docs/gold.json","w") as f:
        json.dump(gold_out, f, indent=2)

    with open(TOP_GENES_TSV,"w",newline="") as f:
        import csv as _c
        w2 = _c.DictWriter(f, delimiter="\t",
             fieldnames=["rank","gene","score","hit","disease_area","druggability","mr_z"],
             extrasaction="ignore")
        w2.writeheader()
        w2.writerows(genes_out)

def log_row(rows, exp_num, result, desc, status):
    row = {
        "experiment":  exp_num,
        "timestamp":   datetime.now().isoformat(timespec="seconds"),
        "recall_at_k": result["recall_at_k"],
        "K":           result["K"],
        "n_hits":      result["n_hits"],
        "n_gold":      result["n_gold"],
        "n_orphan":    result["n_orphan"],
        "top_gene":    result["top_k_genes"][0] if result["top_k_genes"] else "",
        "hit_genes":   "|".join(result.get("hit_genes",[])),
        "status":      status,
        "description": desc,
    }
    rows.append(row)
    header = not os.path.exists(RESULTS_TSV) or os.path.getsize(RESULTS_TSV)==0
    with open(RESULTS_TSV,"a",newline="") as f:
        import csv as _c
        w = _c.DictWriter(f, fieldnames=row.keys(), delimiter="\t")
        if header: w.writeheader()
        w.writerow(row)
    return row

def main():
    from evaluate import run_evaluation

    parser = argparse.ArgumentParser()
    parser.add_argument("--max",     type=int, default=len(EXPERIMENTS))
    parser.add_argument("--no-push", action="store_true")
    args = parser.parse_args()

    import data as d
    K = min(20, len(d.GOLD_GENE_SET))
    print("=" * 60)
    print("RepurposeAgent — Autoresearch Loop")
    print("=" * 60)
    print(f"Gold standard (discovered): {len(d.GOLD_STANDARD)} genes")
    print(f"Orphan candidates:          {len(d.ORPHAN_GENES)} genes")
    print(f"Metric: recall@{K}")
    print(f"Experiments: {min(args.max, len(EXPERIMENTS))}")
    print(f"Auto-push: {'off' if args.no_push else 'on → github.com/RitataLU/repurpose-agent'}\n")

    # Print gold standard table
    print("Gold standard discovered by API:")
    for gene, info in d.GOLD_STANDARD.items():
        print(f"  {gene:<10} {info['disease']:<30} → {info['drug']}")
    print()

    rows = []
    shutil.copy(SCORER_FILE, SCORER_BAK)

    # Baseline
    base = run_evaluation()
    best_recall = base["recall_at_k"]
    best_result = base
    log_row(rows, 0, base, "baseline: equal weights", "keep")
    save_outputs(rows, best_result)
    git_push(0, best_recall, args.no_push)
    print(f"[0] Baseline  recall@{K}={best_recall:.4f}  hits={base['n_hits']}/{K}")

    for i, exp in enumerate(EXPERIMENTS[:args.max], 1):
        shutil.copy(SCORER_FILE, SCORER_BAK)
        write_scorer(exp)
        try:
            result = run_evaluation()
            recall = result["recall_at_k"]
            if recall > best_recall:
                best_recall, best_result = recall, result
                status, marker = "keep", "✓ IMPROVED"
            else:
                shutil.copy(SCORER_BAK, SCORER_FILE)
                status, marker = "discard", "✗ discard"

            log_row(rows, i, result, exp["desc"], status)
            save_outputs(rows, best_result)
            git_push(i, recall, args.no_push)
            print(f"[{i}] {marker}  recall@{K}={recall:.4f}  "
                  f"hits={result['n_hits']}/{K}  | {exp['desc'][:50]}")
        except Exception as e:
            shutil.copy(SCORER_BAK, SCORER_FILE)
            empty = {"recall_at_k":0,"K":K,"n_hits":0,"n_gold":len(d.GOLD_GENE_SET),
                     "n_orphan":len(d.ORPHAN_GENES),"top_k_genes":[],"hit_genes":[],"all_scores":{}}
            log_row(rows, i, empty, f"crash: {str(e)[:60]}", "crash")
            print(f"[{i}] CRASH: {e}")

    print(f"\n{'='*60}")
    print(f"Best recall@{K}: {best_recall:.4f}  ({int(best_recall*K)}/{K} hits)")
    print(f"Top gene: {list(best_result['all_scores'].keys())[0]}")
    print(f"Dashboard: https://ritataLU.github.io/repurpose-agent/")

if __name__ == "__main__":
    main()
```

**Test:** `python agent.py --max 2 --no-push`

---

## Step 6 — `docs/index.html`

Dark scientific dashboard. Reads `results.json`, `genes.json`, `gold.json`.

Requirements:
- IBM Plex Mono + IBM Plex Sans, Chart.js 4.4.1 from cdnjs
- `#080c18` bg, `#00d4aa` teal, `#cdd6f4` text
- Header: sticky, GitHub link to `https://github.com/RitataLU/repurpose-agent`
- Stat bar: best recall@K | experiments | top gene | orphan genes count
- **Gold standard panel**: table showing what `discover.py` found —
  columns: Gene | Disease | Drug | Genetics score
  This shows the audience the gold standard was discovered, not hardcoded
- Experiment trace chart: line, green dots = keep, red = discard
- Ranked gene table: sortable, gold-standard badge, disease area pill
- Methodology accordion
- Empty-state when JSONs are `[]`
- Relative paths only (`./results.json`), works on file:// and GitHub Pages

---

## Step 7 — `docs/showcase.html`

Public editorial page. Reads same JSON files.

Requirements:
- Fraunces + DM Sans, light theme `#f7f4ef` bg, `#00966d` green
- Sections: hero / gap SVG / how it works (3 cards) / loop SVG /
  **gold standard discovery section** (show the table of what the API returned) /
  results panel / implications / CTA
- GitHub link: `https://github.com/RitataLU/repurpose-agent`
- Mobile responsive, all data from fetch('./results.json') etc.

---

## Step 8 — `README.md`

```markdown
# RepurposeAgent

Autonomous drug repurposing agent — discovers orphan genes and gold standard
entirely from public APIs, then iterates a scoring formula to rank candidates.

**Live results:** https://ritataLU.github.io/repurpose-agent/

## How it works
1. `discover.py` queries Open Targets for approved drugs per disease (gold standard),
   then GWAS Catalog for orphan genes (no drug), enriches with biological features
2. `agent.py` runs autoresearch loop: edit scorer → measure recall@K → keep/discard
3. Results auto-push to GitHub Pages after every experiment

## Quick start
python discover.py    # ~8 min, builds data.py from real APIs
python agent.py       # 10 experiments, auto-pushes to GitHub

## Disease scope
Alzheimer's · Parkinson's · coronary artery disease ·
type 2 diabetes · schizophrenia · inflammatory bowel disease

## Data sources
GWAS Catalog (EMBL-EBI) · Open Targets Platform · PubMed E-utilities
```

---

## Step 9 — Run

```bash
python discover.py                    # builds data.py (~8 min)
python discover.py --fallback         # if APIs blocked

python -c "from evaluate import run_evaluation; r=run_evaluation(); \
  print('gold:', r['n_gold'], 'orphan:', r['n_orphan'], 'recall:', r['recall_at_k'])"

python agent.py                       # 10 experiments + auto push
python agent.py --no-push --max 3    # local test
```

---

## Step 10 — GitHub Pages

```
github.com/RitataLU/repurpose-agent
→ Settings → Pages → Branch: main → Folder: /docs → Save
```

Live at: `https://ritataLU.github.io/repurpose-agent/`

---

## File checklist

```
repurpose_agent/
├── CLAUDE.md
├── README.md
├── .gitignore
├── discover.py          ← Phase1 + Phase2 + Phase3 → writes data.py
├── data.py              ← written by discover.py ONLY
├── scorer.py            ← agent edits this
├── evaluate.py          ← never modified
├── agent.py             ← loop + git push
└── docs/
    ├── index.html       ← live dashboard
    ├── showcase.html    ← public page
    ├── results.json     ← written by agent
    ├── genes.json       ← written by agent
    └── gold.json        ← written by agent (gold standard table)
```

## Token efficiency

- Build Python files in one pass: discover → scorer → evaluate → agent
- Test each immediately after: `python -c "import X; print('ok')"`
- HTML files last, one complete file per message
- Fix failures by patching lines only — never regenerate a whole file
- `data.py` is never created manually — only discover.py writes it
