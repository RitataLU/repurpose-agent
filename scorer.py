"""scorer.py — upweight causal features: MR + druggability — auto-written by agent.py"""

WEIGHTS = {
    "mr_z_score": 3.0,
    "druggability_score": 3.0,
    "open_targets_score": 2.0,
    "gwas_pval_log10": 1.5,
    "n_gwas_studies": 1.0,
    "tissue_specificity": 1.0,
    "ppi_degree": 0.5,
    "pubmed_count_5yr": 0.3,
    "eqtl_effect": 1.5
}
BURDEN_MULTIPLIER = False
INTERACTION_TERMS = []

def score_gene(features: dict) -> float:
    score = sum(WEIGHTS.get(k,0)*v for k,v in features.items()
                if k in WEIGHTS and isinstance(v,(int,float)))
    if BURDEN_MULTIPLIER:
        score *= (1 + features.get("burden_daly_m",0)/200)
    for a,b,w in INTERACTION_TERMS:
        score += w * features.get(a,0) * features.get(b,0)
    return score
