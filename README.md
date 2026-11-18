# Project 2: E3 Ubiquitin-Ligase Substrate Recognition: What Do Different Representations Capture?
## CSE 559 - G17

Our main question is to determine which E3 representation (sequence alignment, ESM-2 embeddings, or whole-protein domain Jaccard) best identifies E3s that share substrates closest together, and whether this reflects substrate recognition or just E3 class and family. For each human E3 with a known substrate, we rank all other E3s by similarity and check how often the closest ones share a substrate, compared with a shuffled baseline. As an extension, we repeat this using only the catalytic and recognition regions of each protein.

## Repository layout
| Folder | Contents | Details |
|---|---|---|
| `data_prep/` | E3-substrate pair table construction, cross-check against the mentor workbook, evaluation-ready pairs | `data_prep/README.md` |
| `alignment/` | Smith-Waterman all-vs-all similarity, substrate-sharing evaluation, shuffled baselines, controls, figures | `alignment/README.md` |

## Status
| Component | Status |
|---|---|
| Deduplicated human E3-substrate pair table (3,163 pairs; 2,994 evaluation-ready, 374 E3s) | Done, cross-checked against the mentor workbook (all 3,110 mentor pairs recovered) |
| Sequence alignment (Smith-Waterman, BLOSUM62) similarity and evaluation | Done |
| ESM-2 embeddings (mean-pooled, cosine similarity) | Not yet in the repository |
| Whole-protein InterPro domain Jaccard | Not yet in the repository |
| Extensions (additional embedding model, clustering, recognition/catalytic windows) | Not yet in the repository |

## Preliminary result (alignment only)
- For each of 374 E3s with at least one documented substrate, the top-5 alignment neighbors share a substrate 14.2% of the time, compared with 5.1% for random E3 pairs and 6.7% under a null that controls for how well-studied each E3 is.
- The signal remains about twice chance after removing close paralogs and after setting aside RBX1, RNF7 and ANAPC11.
- A baseline that ignores sequence and always returns the most-studied E3s scores higher (0.29), so the "shares at least one substrate" metric must be read together with popularity-controlled nulls.
- Full tables, figures and the per-class breakdown are in `alignment/results/`.

## Team and contributions
| Member | Contribution |
|---|---|
| Rahul Burjanroppa | Repository and project description; sequence-alignment arm and evaluation |
| Zichen Zhao | E3-substrate pair table construction and validation (`data_prep/`) |


## Data
Mentor-provided input files are expected in `alignment/data/` with their original filenames (see `alignment/README.md` for the list and which script uses each).