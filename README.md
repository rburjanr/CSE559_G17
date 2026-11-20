# Project 2: E3 Ubiquitin-Ligase Substrate Recognition: What Do Different Representations Capture?
## CSE 559 - G17

Our main question is to determine which E3 representation (sequence alignment, ESM-2 embeddings, or whole-protein domain Jaccard) best identifies E3s that share substrates closest together, and whether this reflects substrate recognition or just E3 class and family. For each human E3 with a known substrate, we rank all other E3s by similarity and check how often the closest ones share a substrate, compared with a shuffled baseline. As an extension, we repeat this using only the catalytic and recognition regions of each protein.

## Repository layout
| Folder | Contents | Details |
|---|---|---|
| `data_prep/` | E3-substrate pair table construction, cross-check against the mentor workbook, evaluation-ready pairs | `data_prep/README.md` |
| `alignment/` | Smith-Waterman all-vs-all similarity, substrate-sharing evaluation, shuffled baselines, controls, figures | `alignment/README.md` |
| `evaluation/` | Shared representation evaluator, retrieval figure and result tables | `evaluation/README.md` |

## Status
| Component | Status |
|---|---|
| Deduplicated human E3-substrate pair table (3,163 pairs; 2,994 evaluation-ready, 374 E3s) | Done, cross-checked against the mentor workbook (all 3,110 mentor pairs recovered) |
| Sequence alignment (Smith-Waterman, BLOSUM62) similarity and evaluation | Done |
| Shared representation evaluator | Implemented and verified on alignment at k = 5, 10, and 15 |
| ESM-2 embeddings (mean-pooled, cosine similarity) | Not yet in the repository |
| Whole-protein InterPro domain Jaccard | Not yet in the repository |
| Extensions (additional embedding model, clustering, recognition/catalytic windows) | Not yet in the repository |

## Preliminary results (alignment only)

All cutoffs use the same 374 E3s with documented substrates.

| Cutoff | Alignment precision | Alignment hit rate | Popularity-controlled shuffle precision | Popularity ranking precision |
|---|---:|---:|---:|---:|
| 5 | 14.2% | 37.4% | 6.7% | 29.1% |
| 10 | 11.0% | 47.6% | 6.2% | 26.8% |
| 15 | 9.7% | 51.3% | 6.0% | 23.6% |

Alignment exceeds the popularity-controlled shuffle at every cutoff. Ranking
by documented substrate count gives higher precision at every cutoff.
At k = 5, mean substrate-set Jaccard is 1.82% for alignment and 0.48% for
popularity ranking.

Run `python3 evaluation/evaluate.py` to regenerate the [figure](evaluation/results/retrieval.png)
and [summary table](evaluation/results/summary.csv). See the [evaluation guide](evaluation/README.md)
for accession-indexed matrix inputs and the optional shared paralog filter.
Cross-representation and window results depend on the corresponding matrices.
The [milestone report](milestone-3-report.txt) includes the original alignment controls and the updated retrieval results.

## Team and contributions
| Member | Contribution |
|---|---|
| Rahul Burjanroppa | Repository and project description; sequence-alignment arm and evaluation |
| Zichen Zhao | E3-substrate pair table construction and validation (`data_prep/`) |
| Zhen Liu | Shared representation evaluator, comparative baselines and multi-cutoff figure (`evaluation/`) |


## Data
Mentor-provided input files are expected in `alignment/data/` with their original filenames (see `alignment/README.md` for the list and which script uses each).
