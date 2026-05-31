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
        "recall_at_k":          round(len(hits) / K, 4),
        "K":                    K,
        "n_hits":               len(hits),
        "n_gold":               len(d.GOLD_GENE_SET),
        "n_orphan":             len(d.ORPHAN_GENES),
        "top_k_genes":          top_K,
        "hit_genes":            hits,
        "all_scores":           dict(ranked),
        "weights_snapshot":     dict(sc.WEIGHTS),
        "interaction_snapshot": list(sc.INTERACTION_TERMS),
        "burden_mul":           sc.BURDEN_MULTIPLIER,
    }
