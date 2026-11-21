"""Compare accession-indexed E3 similarity matrices on documented substrate sharing."""
import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "alignment/src"))
from evalcore import coarse_class, share_matrix, topk_neighbors

ALIGNMENT = ROOT / "alignment/results"
ACCESSION = "UniProt accession"
CUTOFFS = (5, 10, 15)
DRAWS = 1000
METRICS = ("precision", "hit", "substrate_jaccard", "same_class")
POPULARITY = "Popularity ranking"


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def load_matrix(path, index_path, allow_duplicates=False):
    accessions = [row[ACCESSION].strip() for row in read_csv(index_path)]
    matrix = np.load(path, allow_pickle=False).astype(float)
    if matrix.shape != (len(accessions), len(accessions)):
        raise ValueError(f"{path}: expected a square matrix matching {index_path}")
    if not all(accessions) or (not allow_duplicates and len(set(accessions)) != len(accessions)):
        raise ValueError(f"{index_path}: expected unique, nonempty accessions")
    indices = {}
    for position, accession in enumerate(accessions):
        indices.setdefault(accession, position)
    return matrix, indices


def permutation_draws(groups, rng):
    blocks = [np.flatnonzero(groups == group) for group in np.unique(groups)]
    draws = np.tile(np.arange(len(groups)), (DRAWS, 1))
    for draw in draws:
        for block in blocks:
            draw[block] = rng.permutation(block)
    return draws


def query_scores(neighbors, shared, jaccard, same_class):
    row = np.arange(len(neighbors))[:, None]
    overlaps = shared[row, neighbors]
    return np.column_stack((overlaps.mean(1), overlaps.any(1),
                            jaccard[row, neighbors].mean(1), same_class[row, neighbors].mean(1)))


def evaluate(cohort, matrices, substrate_sets, classes, eligible, scenario, k):
    metrics = [f"{metric}_at_{k}" for metric in METRICS]
    active = eligible.sum(axis=1) >= k
    if not active.any():
        return [], []
    shared = share_matrix(substrate_sets)
    jaccard = np.array([[len(a & b) / len(a | b) for b in substrate_sets] for a in substrate_sets])
    same_class = classes[:, None] == classes
    degrees = np.array([len(group) for group in substrate_sets])
    rankings = {**matrices, POPULARITY: np.broadcast_to(degrees, eligible.shape)}
    null_draws = {
        "simple": permutation_draws(np.zeros(len(cohort)), np.random.default_rng(0)),
        "class": permutation_draws(classes, np.random.default_rng(0)),
        "popularity": permutation_draws(np.digitize(degrees, [2, 3, 6, 11, 31]),
                                        np.random.default_rng(0)),
    }
    bootstrap = np.random.default_rng(0).integers(active.sum(), size=(DRAWS, active.sum()))
    queries = np.flatnonzero(active)
    summary, per_query = [], []
    reference = None
    for method, matrix in rankings.items():
        neighbors = topk_neighbors(np.where(eligible, matrix, -np.inf), k)
        observed = query_scores(neighbors, shared, jaccard, same_class)[active]
        if reference is None:
            reference = observed
        delta = observed - reference
        lo, hi = np.percentile(delta[bootstrap].mean(axis=1), [2.5, 97.5], axis=0)
        row = {"scenario": scenario, "k": k, "method": method, "queries": int(active.sum()),
               "candidate_pool": len(cohort)}
        for i, metric in enumerate(metrics):
            row[metric] = observed[:, i].mean()
            row[f"{metric}_delta"] = delta[:, i].mean()
            row[f"{metric}_delta_lo"] = lo[i]
            row[f"{metric}_delta_hi"] = hi[i]
        if method != POPULARITY:
            for name, draws in null_draws.items():
                neighbor_labels = draws[:, neighbors[active]]
                query_labels = draws[:, queries, None]
                overlap = shared[query_labels, neighbor_labels]
                null_values = (overlap.mean(axis=(1, 2)), overlap.any(axis=2).mean(axis=1),
                               jaccard[query_labels, neighbor_labels].mean(axis=(1, 2)))
                for metric, values in zip(metrics[:3], null_values):
                    prefix = f"{metric}_{name}_null"
                    row[prefix] = values.mean()
                    row[f"{prefix}_lo"], row[f"{prefix}_hi"] = np.percentile(values, [2.5, 97.5])
                    row[f"{prefix}_p"] = (1 + (values >= row[metric]).sum()) / (DRAWS + 1)
        summary.append(row)
        for query, scores in zip(queries, observed):
            per_query.append({
                "scenario": scenario, "k": k, "method": method, "accession": cohort[query],
                "class": classes[query], "substrate_count": degrees[query],
                "eligible_neighbors": int(eligible[query].sum()),
                "neighbors": " ".join(cohort[i] for i in neighbors[query]),
                **dict(zip(metrics, scores)),
            })
    return summary, per_query


def write_csv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def summarize_by_class(per_query):
    groups = defaultdict(list)
    for row in per_query:
        groups[row["scenario"], row["k"], row["method"], row["class"]].append(row)
    return [
        {"scenario": scenario, "k": k, "method": method, "class": cls, "queries": len(rows),
         **{metric: np.mean([row[metric] for row in rows]) for metric in (f"{name}_at_{k}" for name in METRICS)}}
        for (scenario, k, method, cls), rows in groups.items()
    ]


def plot(summary, path):
    cutoffs = sorted({row["k"] for row in summary if row["scenario"] == "All candidates"})
    fig, axes = plt.subplots(len(cutoffs), 2, figsize=(10, 2.3 * len(cutoffs) + 1),
                             squeeze=False)
    for k, pair in zip(cutoffs, axes):
        rows = [row for row in summary if row["scenario"] == "All candidates" and row["k"] == k]
        y = np.arange(len(rows))
        for ax, base, title in zip(pair, METRICS, ("Precision", "Hit rate")):
            metric = f"{base}_at_{k}"
            values = np.array([row[metric] for row in rows])
            ax.barh(y, values, height=0.52, color=["#16717c"] * (len(rows) - 1) + ["#9ca7b3"])
            prefix = f"{metric}_popularity_null"
            for i, row in enumerate(rows[:-1]):
                ax.hlines(i, row[f"{prefix}_lo"], row[f"{prefix}_hi"], color="#453015")
                ax.plot(row[prefix], i, "o", color="#ffbe58", markeredgecolor="#453015",
                        label="Popularity-controlled shuffle (95% interval)" if i == 0 else None)
            for i, value in enumerate(values):
                ax.text(value + .015, i, f"{value:.1%}", va="center", fontsize=10)
            ax.set(yticks=y, yticklabels=[row["method"] for row in rows],
                   xlim=(0, 1.12), title=f"{title}@{k}")
            ax.set_xticks(np.linspace(0, 1, 6))
            ax.xaxis.set_major_formatter(PercentFormatter(1))
            ax.invert_yaxis()
            ax.grid(axis="x", alpha=.2)
            ax.set_axisbelow(True)
            for side in ("top", "right", "left"):
                ax.spines[side].set_visible(False)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", frameon=False)
    fig.suptitle(f"Retrieval of documented shared substrates | {rows[0]['queries']} E3 queries",
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0.13, 1, 0.90))
    fig.savefig(path, dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", nargs=3, action="append", default=[],
                        metavar=("NAME", "SCORES_NPY", "INDEX_CSV"))
    parser.add_argument("--paralogs", type=Path, help="CSV with accession_a and accession_b")
    parser.add_argument("--out", type=Path, default=ROOT / "evaluation/results")
    args = parser.parse_args()

    roster = {}
    for row in read_csv(ALIGNMENT / "e3_index.csv"):
        roster.setdefault(row[ACCESSION], coarse_class(row))
    substrates = defaultdict(set)
    for row in read_csv(ROOT / "data_prep/evaluation_ready_pairs_2994.csv"):
        substrates[row["e3_accession"]].add(row["substrate_accession"])
    labeled = roster.keys() & substrates.keys()
    matrices = {"Smith-Waterman": load_matrix(ALIGNMENT / "sw_norm_sqrt.npy",
                                              ALIGNMENT / "e3_index.csv", allow_duplicates=True)}
    for name, path, index in args.matrix:
        if name in matrices or name == POPULARITY:
            parser.error(f"Duplicate or reserved method name: {name}")
        matrices[name] = load_matrix(path, index)

    # The roster order is the stable tie-breaker for every method.
    cohort = [a for a in roster if a in labeled and all(a in index for _, index in matrices.values())]
    if len(cohort) <= min(CUTOFFS):
        parser.error(f"At least {min(CUTOFFS) + 1} common labeled E3s are required")
    coverage = []
    for name, (matrix, index) in matrices.items():
        coverage.append({"method": name, "labeled_e3s": len(labeled),
                         "covered_labeled_e3s": len(labeled & index.keys()),
                         "shared_cohort_e3s": len(cohort)})
        order = [index[a] for a in cohort]
        selected = matrix[np.ix_(order, order)]
        if not np.isfinite(selected[~np.eye(len(cohort), dtype=bool)]).all():
            parser.error(f"{name}: nonfinite off-diagonal scores on the shared cohort")
        matrices[name] = selected

    classes = np.array([roster[a] for a in cohort])
    eligible = ~np.eye(len(cohort), dtype=bool)
    same_class = classes[:, None] == classes
    scenarios = {"All candidates": eligible, "Within class": eligible & same_class}
    if args.paralogs:
        positions = {accession: i for i, accession in enumerate(cohort)}
        filtered = eligible.copy()
        for row in read_csv(args.paralogs):
            a, b = row["accession_a"].strip(), row["accession_b"].strip()
            if a in positions and b in positions:
                i, j = positions[a], positions[b]
                filtered[i, j] = filtered[j, i] = False
        scenarios["Exclude annotated paralogs"] = filtered
    for cutoff in (0.5, 0.3, 0.2):
        scenarios[f"Alignment score below {cutoff}"] = eligible & (matrices["Smith-Waterman"] < cutoff)

    summary, per_query = [], []
    for scenario, mask in scenarios.items():
        for k in CUTOFFS:
            results, queries = evaluate(cohort, matrices, [substrates[a] for a in cohort],
                                        classes, mask, scenario, k)
            summary.extend(results)
            per_query.extend(queries)
    args.out.mkdir(parents=True, exist_ok=True)
    write_csv(args.out / "summary.csv", summary)
    write_csv(args.out / "per_query.csv", per_query)
    write_csv(args.out / "coverage.csv", coverage)
    write_csv(args.out / "by_class.csv", summarize_by_class(per_query))
    plot(summary, args.out / "retrieval.png")


if __name__ == "__main__":
    main()
