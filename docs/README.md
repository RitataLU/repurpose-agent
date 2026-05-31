# RepurposeAgent — Results

An autonomous agent that discovers and ranks orphan GWAS genes by drug repurposing potential,
using only free public APIs and an iterative scoring loop.

## Pages
- [Live Dashboard](index.html) — experiment log, gene rankings, gold standard table
- [Results Visualization](results.html) — full experiment story with all charts and narrative
- [Project Showcase](showcase.html) — editorial overview of the problem, method, and findings

## What we found

An agent analysed **169 GWAS-implicated genes** with no approved drug across 6 diseases
(Alzheimer's, Parkinson's, coronary artery disease, type 2 diabetes, schizophrenia, IBD).

Using **Spearman rank correlation** as the metric — measuring how well the scorer's ranking
matches each orphan gene's cosine similarity to a centroid of 13 curated gold-standard
gene-drug pairs — the agent ran 10 scoring experiments:

- **Baseline Spearman r:** 0.4659 (equal-weight formula)
- **Best Spearman r:** 0.5237 (experiment 10: OT score boost + GWAS replication signal)
- **Improvement:** +0.0578 over baseline
- **Top orphan predictions:** KCNQ1, PLG, FTO

The best formula up-weights Open Targets association score (4.0×), GWAS study count (3.0×),
and MR z-score (3.0×) — confirming that the strongest signal comes from causal genetics and
established functional evidence, not from literature popularity.

## How the metric works

1. **Normalize** all orphan + gold gene features together (min-max, same scale)
2. **Gold centroid** = mean of 13 normalized gold gene feature vectors
3. **Cosine similarity** = how close each orphan gene is to the centroid
4. **Spearman r** = rank correlation between scorer ranking and cosine-similarity ranking

A higher Spearman r means the scorer ranks genes the same way that biological similarity to
proven drug targets would suggest.

## Data sources
- [GWAS Catalog](https://www.ebi.ac.uk/gwas/) — genome-wide significant associations (p < 5×10⁻⁸)
- [Open Targets Platform](https://platform.opentargets.org/) — association scores, druggability
- [PubMed E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25499/) — publication counts

## GitHub
**Repository:** https://github.com/RitataLU/repurpose-agent
