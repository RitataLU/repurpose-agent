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
