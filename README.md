## E3-substrate data validation 

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

### Suggested Milestone 3 wording
We integrated the raw UbiBrowser known E3-substrate records with the UbiNet archive, restricted the data to human records, standardized accession fields, and deduplicated E3-substrate accession pairs. This yielded 3,163 unique pairs. Cross-checking against the mentor-provided workbook recovered all 3,110 curated known pairs. 
We identified 53 additional raw-derived pairs and retained them for discrepancy review rather than silently adding them to the primary evaluation set. For downstream nearest-neighbor evaluation, we use the 2,994 mentor-confirmed pairs whose E3 is explicitly marked as present in the 679-member E3-ome roster.
