# Shared E3 substrate-retrieval evaluation

Run from the repository root with NumPy and Matplotlib installed.

```bash
python3 evaluation/evaluate.py
```

The default run uses the existing normalized Smith–Waterman matrix and the
2,994 evaluation-ready pairs. To compare other representations, pass each
similarity matrix and its row/column index. The following paths are examples.

```bash
python3 evaluation/evaluate.py \
  --matrix ESM-2 esm2.npy esm2_index.csv \
  --matrix "Domain Jaccard" domains.npy domains_index.csv \
  --paralogs paralog_pairs.csv \
  --out evaluation/results_comparison
```

Each NumPy matrix must be square, with higher scores indicating greater
similarity. Its CSV index must have a unique `UniProt accession` per row, in the
order used on both axes. Off-diagonal scores must be finite on the shared
cohort. Omit absent accessions from the index and matrix. Ties use the
original roster order. The optional paralog CSV takes `accession_a` and
`accession_b` columns, with each pair excluded in both directions for all
representations and benchmarks in the annotated-paralog scenario.

For the region comparison, pass the recognition-window and catalytic-window
similarity matrices in one separate invocation and select a separate output
directory. The evaluator will compare both with whole-protein Smith–Waterman
on E3s covered by all three matrices. Include the matching whole-protein
matrix as an additional input if the window scores use another representation.

## Outputs

- `retrieval.png` plots precision and hit rate at k = 5, 10, and 15 for unconstrained retrieval. Dots
  and intervals show the popularity-controlled permutation mean and central
  95% of permutation results.
- `summary.csv` reports precision, hit rate, mean substrate-set Jaccard, and
  same-class retrieval rate. It includes the popularity ranking, three
  permutation nulls, and paired query-bootstrap differences against
  Smith–Waterman. Each resampling calculation uses 1,000 draws with seed 0.
- `coverage.csv` reports labeled accessions available for each representation
  and the intersection used for comparison.
- `per_query.csv` records metrics and retrieved accessions per query.
- `by_class.csv` reports mean observed metrics and query counts by class.

Each control uses a shared candidate mask. Queries with fewer than k
eligible candidates are omitted from that scenario's averages. The
alignment-score masks are sensitivity analyses. The annotated-paralog
control requires a team-supplied annotation file.

All metrics use documented substrate sets. The popularity ranking uses the
full documented substrate counts available for this retrospective comparison.
