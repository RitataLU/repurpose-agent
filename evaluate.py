"""
evaluate.py — Fixed evaluation. Never modified by agent.
Metric: Spearman correlation between scorer ranking and
cosine similarity ranking to gold standard centroid.
Range: -1 to 1. Higher = better.
"""
import importlib, math

def normalize(genes):
    skip = {"disease_area","disease_name","data_source","ensembl_id","disease","drug","area","ot_score","genetics_score","gene"}
    sample = next(iter(genes.values()))
    keys = [k for k in sample if k not in skip and isinstance(sample.get(k),(int,float))]
    normed = {g: dict(f) for g,f in genes.items()}
    for k in keys:
        vals = [genes[g][k] for g in genes if isinstance(genes[g].get(k),(int,float))]
        if not vals: continue
        mn,mx = min(vals),max(vals)
        rng = mx-mn if mx!=mn else 1.0
        for g in normed:
            if isinstance(normed[g].get(k),(int,float)):
                normed[g][k] = (normed[g][k]-mn)/rng
    return normed, keys

def cosine(a, b, keys):
    dot = sum(a.get(k,0)*b.get(k,0) for k in keys)
    na = math.sqrt(sum(a.get(k,0)**2 for k in keys))
    nb = math.sqrt(sum(b.get(k,0)**2 for k in keys))
    return dot/(na*nb) if na*nb > 0 else 0.0

def spearman(rank_a, rank_b):
    n = len(rank_a)
    if n < 3: return 0.0
    d2 = sum((rank_a[i]-rank_b[i])**2 for i in range(n))
    return 1 - (6*d2)/(n*(n*n-1))

def run_evaluation():
    import data as d, scorer as sc
    importlib.reload(d); importlib.reload(sc)
    # Normalize orphan + gold using SAME min/max so vectors are comparable
    combined = {**d.ORPHAN_GENES, **d.GOLD_FEATURES}
    normed, keys = normalize(combined)
    gold_vecs = [normed[g] for g in d.GOLD_FEATURES if g in normed]
    centroid = {k: sum(v.get(k,0) for v in gold_vecs)/len(gold_vecs) for k in keys}
    orphan_normed = {g: normed[g] for g in d.ORPHAN_GENES if g in normed}
    genes_list = list(orphan_normed.keys())
    sims = {g: cosine(orphan_normed[g], centroid, keys) for g in genes_list}
    scores = {g: sc.score_gene(orphan_normed[g]) for g in genes_list}
    sim_sorted = sorted(genes_list, key=lambda g: sims[g], reverse=True)
    sc_sorted = sorted(genes_list, key=lambda g: scores[g], reverse=True)
    sim_rank = {g:i for i,g in enumerate(sim_sorted)}
    sc_rank = {g:i for i,g in enumerate(sc_sorted)}
    r = spearman([sim_rank[g] for g in genes_list], [sc_rank[g] for g in genes_list])
    return {
        "spearman_r": round(r, 4),
        "n_orphan": len(genes_list),
        "n_gold": len(d.GOLD_FEATURES),
        "top_genes": sc_sorted[:20],
        "all_scores": {g: scores[g] for g in sc_sorted},
        "all_sims": {g: sims[g] for g in sim_sorted},
        "weights_snapshot": dict(sc.WEIGHTS),
        "interaction_snapshot": list(sc.INTERACTION_TERMS),
        "burden_mul": sc.BURDEN_MULTIPLIER,
    }
