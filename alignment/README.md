# E3 substrate recognition: sequence-alignment arm (Rahul Burjanroppa)

Smith-Waterman similarity between 678 human E3-ome proteins, evaluated by whether
nearest-neighbor E3s share documented substrates. This directory holds the alignment
and evaluation code; the ESM-2 and domain-Jaccard arms are separate pieces of the team project.

## Data (not included; put the mentor's files in `data/` with their ORIGINAL names, no renaming)
| File in `data/` | Used by |
|---|---|
| `E3_substrate_pairs_and_e3_roster.xlsx` (sheets "Known pairs", "E3 roster" with sequences) | all scripts |
| `UbiBrowser_known_E3_substrate_records_raw.txt` | `build_pairs.py` |
| `UbiNet_E3_substrate_archive_2024-09-20.txt` | `build_pairs.py` |
| `UbiBrowser_chain_annotated_records.txt`, `E3ome_candidate_gene_curation.xlsx` | not yet used |

## Dependencies
Python 3.11+; `pip install pandas numpy openpyxl biopython parasail matplotlib`

## Run (about 1 minute total)
```
python src/build_pairs.py          # raw UbiBrowser+UbiNet -> deduplicated human pairs; cross-check vs mentor workbook
python src/align_all.py            # all-vs-all Smith-Waterman (14 s on 2 cores) -> results/sw_*.npy, e3_index.csv
python src/evaluate_sw.py          # neighbors vs substrate sharing; shuffles; paralog / heavily-studied controls
python src/evaluate_degree_null.py # popularity-controlled null (shuffle within substrate-count bins)
python src/extra_analyses.py       # popularity/length baselines, per-class and per-degree tables, example neighbors, auto-ubiquitination check
python src/make_figures.py         # results/fig1_*.png, fig2_*.png
```
`src/evalcore.py` holds shared helpers. All randomness uses fixed seed 0, 1,000 permutations.

## Settings (fixed before any evaluation)
- Local alignment (Smith-Waterman), BLOSUM62, gap open 12 / extend 1 (parasail convention: gap of length k costs 12+(k-1); equals BLAST 11+k).
- Primary normalization: score / sqrt(self_i * self_j). Sensitivity: score / min(self), raw score.
- parasail was validated against Biopython PairwiseAligner: identical scores on 200/200 random pairs.
- Evaluation set: E3s with >=1 documented substrate (374 of 677 unique accessions). Missing pairs are unknown, not negatives.
- TRIM49D1/TRIM49D2 share accession C9J1S8 and an identical sequence; collapsed to one E3 for evaluation. LOC122513141 has no sequence and is excluded.

## Outputs in `results/`
`pairs_dedup_human.csv`, `discrepancy_*.csv`, `sw_norm_sqrt.npy` (+ `sw_norm_min.npy`, `sw_raw.npy`; rows follow `e3_index.csv`),
`e3_index.csv`, `sw_main_neighbors.npy`, `sw_evaluation.csv`, `sw_degree_controlled_null.csv`, `sw_by_degree_bin.csv`,
`sw_per_query.csv`, `sw_extra_analyses.txt`, `fig1_sw_vs_baselines.png`, `fig2_sensitivity.png`.