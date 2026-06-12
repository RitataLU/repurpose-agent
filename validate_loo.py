"""
validate_loo.py — Leave-one-out cross-validation on the gold standard.
Usage: python3 validate_loo.py
Output: docs/loo.json
"""
import json, copy, sys, os

sys.path.insert(0, os.path.dirname(__file__))

from evaluate import normalize, cosine, spearman

RESULTS_JSON = "docs/results.json"
LOO_JSON     = "docs/loo.json"

def load_winning_weights():
    """Return weights from the best-spearman_r keep entry, or read scorer.py directly."""
    with open(RESULTS_JSON) as f:
        results = json.load(f)
    candidates = [r for r in results if r.get("status") == "keep" and r.get("spearman_r") is not None]
    best = max(candidates, key=lambda r: r.get("spearman_r", 0)) if candidates else (results[0] if results else {})
    weights = best.get("weights_snapshot", {})
    interaction = best.get("interaction_snapshot", [])
    burden_mul = best.get("burden_mul", False)

    # If results.json doesn't carry weights (older format), import scorer.py directly
    if not weights:
        import scorer as sc
        weights = dict(sc.WEIGHTS)
        interaction = list(sc.INTERACTION_TERMS)
        burden_mul = sc.BURDEN_MULTIPLIER
        print(f"weights_snapshot absent in results.json — reading scorer.py directly")
    else:
        print(f"Using weights from experiment {best.get('experiment','?')} (spearman_r={best.get('spearman_r','?')})")
    return weights, interaction, burden_mul

def score_gene_with(features, weights, interaction, burden_mul):
    score = sum(weights.get(k, 0) * v for k, v in features.items()
                if k in weights and isinstance(v, (int, float)))
    if burden_mul:
        score *= (1 + features.get("burden_daly_m", 0) / 200)
    for a, b, w in interaction:
        score += w * features.get(a, 0) * features.get(b, 0)
    return score

def main():
    import data as d

    gold_features = dict(d.GOLD_FEATURES)
    orphan_genes  = dict(d.ORPHAN_GENES)
    gold_genes    = list(gold_features.keys())
    n_total       = len(gold_genes)

    weights, interaction, burden_mul = load_winning_weights()

    print(f"\nLeave-one-out CV on {n_total} gold standard genes\n")

    per_gene = []
    recovered = 0

    for held_out in gold_genes:
        # Build training gold (remaining genes)
        remaining_gold = {g: f for g, f in gold_features.items() if g != held_out}

        # Candidate pool = orphans + held-out gene
        candidates = copy.deepcopy(orphan_genes)
        candidates[held_out] = copy.deepcopy(gold_features[held_out])

        # Normalize over combined space (remaining gold + candidates)
        combined = {**remaining_gold, **candidates}
        normed, keys = normalize(combined)

        # Centroid from remaining 12 gold genes
        gold_vecs = [normed[g] for g in remaining_gold if g in normed]
        if not gold_vecs:
            print(f"  {held_out}: no remaining gold — skip")
            continue
        centroid = {k: sum(v.get(k, 0) for v in gold_vecs) / len(gold_vecs) for k in keys}

        # Score all candidates
        cand_normed = {g: normed[g] for g in candidates if g in normed}
        scores = {g: score_gene_with(cand_normed[g], weights, interaction, burden_mul)
                  for g in cand_normed}

        ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)
        rank_map = {g: i + 1 for i, (g, _) in enumerate(ranked)}

        rank = rank_map.get(held_out, len(ranked))
        rec  = rank <= 20
        if rec:
            recovered += 1

        per_gene.append({"gene": held_out, "rank": rank, "recovered": rec})
        status = "✓" if rec else "✗"
        print(f"  {status} {held_out:<10}  rank={rank:>3}  {'TOP-20' if rec else 'missed'}")

    rate = round(recovered / n_total, 3) if n_total else 0.0
    result = {
        "n_total": n_total,
        "recovered_top_20": recovered,
        "recovery_rate": rate,
        "per_gene": per_gene,
    }

    with open(LOO_JSON, "w") as f:
        json.dump(result, f, indent=2)

    print(f"\nLOO recovery: {recovered}/{n_total} = {rate:.1%}")
    print(f"Wrote {LOO_JSON}")

if __name__ == "__main__":
    main()
