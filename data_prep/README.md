# Data preparation and validation (Zichen Zhao) 

## Recommended evaluation input
`evaluation_ready_pairs_2994.csv` contains the 2,994 mentor-confirmed Known pairs whose E3 is marked `Yes` in the 679-E3 roster. This is the conservative evaluation-ready table for downstream nearest-neighbor scoring.

## Cross-check
- Raw UbiBrowser + UbiNet processing produced 3,163 unique human E3-substrate accession pairs.
- All 3,110 pairs in the mentor `Known pairs` sheet were recovered (100% coverage).
- 53 additional raw-derived pairs were not in the mentor Known-pairs sheet.
- Of those 53: 29 have an E3 accession that is not a single accession in the mentor roster (mostly complexes or missing accession), 4 have a roster E3 but missing substrate accession, and 20 have a valid roster E3 and substrate accession but are absent from the mentor Known-pairs sheet.
- To avoid silently changing the mentor-curated evaluation set, the 20 plausible extras are retained for review rather than included in the primary evaluation input.

### Files
- `evaluation_ready_pairs_2994.csv`: recommended downstream evaluation table.
- `mentor_known_pairs_3110.csv`: complete mentor Known-pairs sheet in standardized columns.
- `extra_53_pairs_review.csv`: discrepancy table with reason categories.
- `crosscheck_summary.csv`: counts for reporting.

## Follow-up check (added from the alignment arm)
Matching the 20 plausible extras against UbiNet's organism fields: 16 have non-human substrates (viral, bacterial or yeast, e.g. HPV E7, HBV P, HIV tat), 1 is human (WWP1-EGFR, P00533), and 3 could not be verified from the provided files. This supports keeping them out of the primary evaluation set.

## Running the script
`build_e3_substrate_table.py` currently reads from hard-coded `/mnt/data/...` paths with download-style filenames; edit the `RAW`, `UBINET` and `CUR` paths at the top to point at the mentor files before running.