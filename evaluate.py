"""
evaluate.py — Fixed evaluation. Never modified by agent.
Metric: Spearman correlation between scorer ranking and
cosine similarity ranking to gold standard centroid.
Range: -1 to 1. Higher = better.
"""
import importlib
import math

def normalize(genes):
    skip = {"disease_area","disease_name","data_source","ensembl_id"}
    sample = next(iter(genes.values()))
    keys = [k for k in sample if k not in skip and isinstance(sample.get(k),(int,float))]
    normed = {g: dict(f) for g,f in genes.items()}
    for k in keys:
        vals = [genes[g][k] for g in genes if isinstance(genes[g].get(k),(int,float))]
        mn,mx = min(vals),max(vals)
        rng = mx-mn if mx!=mn else 1.0
        for g in normed:
            if isinstance(normed[g].get(k),(int,float)):
                normed[g][k] = (normed[g][k]-mn)/rng
    return normed, keys

def cosine_sim(vec_a, vec_b, keys):
    dot = sum(vec_a.get(k,0)*vec_b.get(k,0) for k in keys)
    na  = math.sqrt(sum(vec_a.get(k,0)**2 for k in keys))
    nb  = math.sqrt(sum(vec_b.get(k,0)**2 for k in keys))
    return dot/(na*nb) if na*nb > 0 else 0.0

def spearman_r(rank_a, rank_b):
    n = len(rank_a)
    if n < 3: return 0.0
    d2 = sum((rank_a[i]-rank_b[i])**2 for i in range(n))
    return 1 - (6*d2)/(n*(n**2-1))

def run_evaluation():
    import data as d
    import scorer as sc
    importlib.reload(d)
    importlib.reload(sc)

    FEATURE_KEYS = d.FEATURE_KEYS
    normed_orphan, keys = normalize(d.ORPHAN_GENES)

    # Gold centroid from GOLD_STANDARD feature values
    # Use whatever numeric fields exist in GOLD_STANDARD entries
    # Fall back to using open_targets_score and genetics_score as proxies
    gold_vecs = []
    for gene, info in d.GOLD_STANDARD.items():
        vec = {k: info.get(k, info.get("genetics_score", 0.5)) for k in keys}
        gold_vecs.append(vec)
    centroid = {k: sum(v.get(k,0) for v in gold_vecs)/len(gold_vecs) for k in keys}

    # Similarity ranking (ground truth)
    genes = list(normed_orphan.keys())
    sims  = {g: cosine_sim(normed_orphan[g], centroid, keys) for g in genes}
    sim_ranked = sorted(genes, key=lambda g: sims[g], reverse=True)
    sim_rank   = {g: i for i,g in enumerate(sim_ranked)}

    # Scorer ranking (agent)
    scores = {g: sc.score_gene(normed_orphan[g]) for g in genes}
    sc_ranked  = sorted(genes, key=lambda g: scores[g], reverse=True)
    sc_rank    = {g: i for i,g in enumerate(sc_ranked)}

    # Spearman
    ra = [sim_rank[g] for g in genes]
    rb = [sc_rank[g]  for g in genes]
    sr = spearman_r(ra, rb)

    return {
        "spearman_r":      round(sr, 4),
        "n_orphan":        len(genes),
        "n_gold":          len(d.GOLD_STANDARD),
        "top_genes":       sc_ranked[:20],
        "all_scores":      dict(zip(sc_ranked, [scores[g] for g in sc_ranked])),
        "all_sims":        dict(zip(sim_ranked, [sims[g] for g in sim_ranked])),
        "weights_snapshot":dict(sc.WEIGHTS),
    }
