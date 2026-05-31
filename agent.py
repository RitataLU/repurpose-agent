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


def git_push(exp_num: int, spearman: float, no_push: bool):
    files = [RESULTS_TSV, TOP_GENES_TSV, DOCS_RESULTS, DOCS_GENES, SCORER_FILE]
    try:
        subprocess.run(["git", "add"] + files, check=True, capture_output=True)
        msg = f"exp {exp_num}: spearman_r={spearman:.4f} [{datetime.now().strftime('%Y-%m-%d %H:%M')}]"
        r = subprocess.run(["git", "commit", "-m", msg], capture_output=True, text=True)
        if "nothing to commit" in r.stdout:
            return
        if not no_push:
            subprocess.run(["git", "push", "origin", "main"], check=True, capture_output=True)
            print("  → pushed to GitHub")
    except subprocess.CalledProcessError as e:
        print(f"  → git skip: {str(e)[:60]}")


def save_outputs(rows: list, best: dict):
    import data as d
    os.makedirs("docs", exist_ok=True)

    with open(DOCS_RESULTS, "w") as f:
        json.dump(rows, f, indent=2, default=str)

    all_sims = best.get("all_sims", {})
    genes_out = [
        {
            "rank":         i,
            "gene":         g,
            "score":        round(s, 4),
            "gold_sim":     round(all_sims.get(g, 0.0), 4),
            "hit":          g in d.GOLD_GENE_SET,
            "disease_area": d.ORPHAN_GENES.get(g, {}).get("disease_area", ""),
            "disease_name": d.ORPHAN_GENES.get(g, {}).get("disease_name", ""),
            "druggability": round(d.ORPHAN_GENES.get(g, {}).get("druggability_score", 0), 2),
            "mr_z":         round(d.ORPHAN_GENES.get(g, {}).get("mr_z_score", 0), 2),
            "ot_score":     round(d.ORPHAN_GENES.get(g, {}).get("open_targets_score", 0), 2),
            "gwas_p":       round(d.ORPHAN_GENES.get(g, {}).get("gwas_pval_log10", 0), 1),
            "pubmed":       d.ORPHAN_GENES.get(g, {}).get("pubmed_count_5yr", 0),
            "source":       d.ORPHAN_GENES.get(g, {}).get("data_source", ""),
        }
        for i, (g, s) in enumerate(best["all_scores"].items(), 1)
        if g in d.ORPHAN_GENES
    ]
    with open(DOCS_GENES, "w") as f:
        json.dump(genes_out, f, indent=2)

    gold_out = [{"gene": g, **v} for g, v in d.GOLD_STANDARD.items()]
    with open("docs/gold.json", "w") as f:
        json.dump(gold_out, f, indent=2)

    with open(TOP_GENES_TSV, "w", newline="") as f:
        import csv as _c
        w2 = _c.DictWriter(
            f, delimiter="\t",
            fieldnames=["rank", "gene", "score", "hit", "disease_area", "druggability", "mr_z"],
            extrasaction="ignore",
        )
        w2.writeheader()
        w2.writerows(genes_out)


def log_row(rows, exp_num, result, desc, status):
    top_genes = result.get("top_genes") or list(result.get("all_scores", {}).keys())
    row = {
        "experiment":  exp_num,
        "timestamp":   datetime.now().isoformat(timespec="seconds"),
        "spearman_r":  result["spearman_r"],
        "n_gold":      result["n_gold"],
        "n_orphan":    result["n_orphan"],
        "top_gene":    top_genes[0] if top_genes else "",
        "status":      status,
        "description": desc,
    }
    rows.append(row)
    header = not os.path.exists(RESULTS_TSV) or os.path.getsize(RESULTS_TSV) == 0
    with open(RESULTS_TSV, "a", newline="") as f:
        import csv as _c
        w = _c.DictWriter(f, fieldnames=row.keys(), delimiter="\t")
        if header:
            w.writeheader()
        w.writerow(row)
    return row


def main():
    from evaluate import run_evaluation

    parser = argparse.ArgumentParser()
    parser.add_argument("--max",      type=int, default=len(EXPERIMENTS))
    parser.add_argument("--no-push",  action="store_true")
    args = parser.parse_args()

    import data as d
    print("=" * 60)
    print("RepurposeAgent — Autoresearch Loop")
    print("=" * 60)
    print(f"Gold standard (curated):    {len(d.GOLD_STANDARD)} genes")
    print(f"Scored pool (orphan+gold):  {len(d.ORPHAN_GENES)} genes")
    print(f"Metric: spearman_r (scorer vs cosine-sim-to-gold-centroid ranking)")
    print(f"Experiments: {min(args.max, len(EXPERIMENTS))}")
    print(f"Auto-push: {'off' if args.no_push else 'on → github.com/RitataLU/repurpose-agent'}\n")

    print("Gold standard (curated gene-drug pairs):")
    for gene, info in d.GOLD_STANDARD.items():
        print(f"  {gene:<10} {info['disease']:<30} → {info['drug']}")
    print()

    rows = []
    shutil.copy(SCORER_FILE, SCORER_BAK)

    base = run_evaluation()
    best_spearman = base["spearman_r"]
    best_result   = base
    log_row(rows, 0, base, "baseline: equal weights", "keep")
    save_outputs(rows, best_result)
    git_push(0, best_spearman, args.no_push)
    print(f"[0] Baseline  spearman_r={best_spearman:+.4f}")

    for i, exp in enumerate(EXPERIMENTS[:args.max], 1):
        shutil.copy(SCORER_FILE, SCORER_BAK)
        write_scorer(exp)
        try:
            result = run_evaluation()
            spearman = result["spearman_r"]
            if spearman > best_spearman:
                best_spearman, best_result = spearman, result
                status, marker = "keep", "✓ IMPROVED"
            else:
                shutil.copy(SCORER_BAK, SCORER_FILE)
                status, marker = "discard", "✗ discard"

            log_row(rows, i, result, exp["desc"], status)
            save_outputs(rows, best_result)
            git_push(i, spearman, args.no_push)
            print(
                f"[{i}] {marker}  spearman_r={spearman:+.4f}"
                f"  | {exp['desc'][:55]}"
            )
        except Exception as e:
            shutil.copy(SCORER_BAK, SCORER_FILE)
            empty = {
                "spearman_r": 0.0,
                "n_gold": len(d.GOLD_STANDARD), "n_orphan": len(d.ORPHAN_GENES),
                "top_genes": [], "all_scores": {}, "all_sims": {},
                "weights_snapshot": {},
            }
            log_row(rows, i, empty, f"crash: {str(e)[:60]}", "crash")
            print(f"[{i}] CRASH: {e}")

    print(f"\n{'='*60}")
    top_genes = best_result.get("top_genes") or list(best_result.get("all_scores", {}).keys())
    print(f"Best spearman_r: {best_spearman:+.4f}")
    if top_genes:
        print(f"Top 3 orphan genes: {', '.join(top_genes[:3])}")
    print(f"Dashboard: https://ritataLU.github.io/repurpose-agent/")


if __name__ == "__main__":
    main()
