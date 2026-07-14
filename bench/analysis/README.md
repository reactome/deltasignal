# Failure-analysis toolkit

Scripts for breaking down DeltaSignal benchmark misses by root cause —
"where did the signal go wrong?" — so improvements can target the most
common failure mode instead of moving thresholds.

Everything operates on the **per-case dump** that the benchmarks emit when
`DS_DUMP_CASES=/path/to/dump.tsv` is set. The dump is one row per
benchmark case:

    pathway  gene  direction  key_output  predicted  expected  method  pred_ui

`direction` is `0` (KO, perturbed to UI 0) or `2` (OE, perturbed to UI 80).
`predicted`/`expected` are class labels `0` (DOWN) / `1` (NORM) / `2` (UP).
`pred_ui` is the model's continuous output on the 0–100 UI scale.

## Setup

One-time: build the gene → Reactome stable_id cache from a live Neo4j.
Re-run if new genes appear in the dumps.

    python bench/analysis/build_gene_cache.py

By default reads `/tmp/ds_best_cases.tsv` and writes `/tmp/gene_to_stids.json`.

## Scripts

### `classify_failures.py`

Categorize each wrong case by mechanism:

    python bench/analysis/classify_failures.py /tmp/ds_best_cases.tsv \
        /tmp/ds_classified.tsv

Categories:

| Category | Meaning |
|---|---|
| `GENE_NOT_IN_PATHWAY` | Gene has no UUID in the network. |
| `READOUT_NOT_IN_PATHWAY` | Readout entity has no UUID and no proxy. |
| `NO_PATH_GENE_TO_READOUT` | Both endpoints exist, but no directed path. |
| `PATH_EXISTS_DIRECTION_FLIPPED` | Path exists, prediction has opposite sign. |
| `PATH_EXISTS_FALSE_POSITIVE` | Expected NORM, model predicted change. |
| `PATH_EXISTS_SIGNAL_LOST` | Expected non-NORM, model predicted NORM. |

For Mit_G1 and S_Phase cases the joint network (Mit_G1 ∪ S_Phase) is used,
matching what `benchmark_selective_joint.py` evaluates against.

### `compare_configs.py`

Diff two dumps to see what flipped between configurations:

    python bench/analysis/compare_configs.py \
        /tmp/ds_best_cases.tsv /tmp/ds_or_mean.tsv \
        --label-a or_max --label-b or_mean

Reports overall accuracy delta, per-pathway delta, and the lists of cases
that gained, lost, or moved sideways.

### `trace_paths.py`

Print the shortest gene → readout path through the logic network for one
case (or for every case of a category in a classified dump):

    python bench/analysis/trace_paths.py \
        Transcriptional_Regulation_by_TP53 AKT1 4655344

    python bench/analysis/trace_paths.py \
        --from-dump /tmp/ds_classified.tsv PATH_EXISTS_SIGNAL_LOST

Each edge shows polarity, cluster type (and/or/single), edge category, and
the branching factor at the target — useful for spotting "lost in a 10-way
OR" cases.

### `report.py`

Publication-style summary table. Accepts one or more dumps and prints
overall accuracy, confusion matrix, per-pathway accuracy, failure
breakdown (when the gene cache exists), and pred_ui percentiles per class.
With multiple inputs, also prints pairwise per-pathway deltas.

    python bench/analysis/report.py /tmp/ds_or_max.tsv /tmp/ds_or_mean.tsv \
        --labels or_max,or_mean

### `case_dump_overview.py`

Quick at-a-glance summary of a single dump — no network access required.
Reports confusion matrix, per-pathway accuracy, readouts with ≥4 wrong
cases (structural-issue candidates), and per-class pred_ui percentiles.

    python bench/analysis/case_dump_overview.py /tmp/ds_best_cases.tsv

### `threshold_sweep.py`

Sweep classification cutoffs (DOWN < d, UP ≥ u) against a dump's pred_ui
column and find the single GLOBAL pair that maximizes accuracy. A wide
tie band in the top-10 list signals bimodal output distribution — useful
diagnostic for whether thresholds or the propagator is the bottleneck.

    python bench/analysis/threshold_sweep.py /tmp/ds_best_cases.tsv

### `check_silo_bug.py`

Detects siloed stable_ids in a pathway's logic network — cases where
positional decomposition has emitted multiple UUIDs for the same
biological entity without connecting them, so upstream-reaction-output
UUIDs and downstream-reaction-input UUIDs are different nodes and
signal can't flow.

    python bench/analysis/check_silo_bug.py             # all 9 pathways
    python bench/analysis/check_silo_bug.py --effective # slower BFS check
    python bench/analysis/check_silo_bug.py PATHWAY_DIR

This is a diagnostic for an UPSTREAM (`logic-network-generator`) bug.
Re-run after a generator change; the silo% should drop. See
`project_uuid_silo_bug` memory for the full write-up. Known limitation:
"internally fragmented all-bridge" stids (e.g. TP53 Tetramer) aren't
caught — use `trace_paths.py` for those.

### `run_config_sweep.sh`

End-to-end driver: restart the API server under each named config and
run the 9-pathway benchmark, leaving a `/tmp/ds_<label>.tsv` per config
for `compare_configs.py` / `report.py` to consume.

    bench/analysis/run_config_sweep.sh

### `path_summary.py`

Tabulate path topology (hops, cluster mix, OR branching) per case:

    python bench/analysis/path_summary.py /tmp/ds_classified.tsv \
        PATH_EXISTS_SIGNAL_LOST

Use to distinguish long-path attenuation from paralog-OR masking — the
former is intrinsic to the network, the latter is sensitive to the
`DS_OR_MODE` (max vs. mean) configuration.

## Typical workflow

1. Run a benchmark with `DS_DUMP_CASES=/tmp/ds_<config>.tsv`.
2. `classify_failures.py` to get category counts and a classified TSV.
3. `path_summary.py` on `PATH_EXISTS_SIGNAL_LOST` to see whether OR
   branching is dominating (→ try `DS_OR_MODE=mean`) or path length is
   (→ network is structurally lossy, not a config bug).
4. `trace_paths.py` on the worst offenders (largest readout clusters,
   highest OR branching) to read the actual cascade.
5. Try a candidate config, dump again, `compare_configs.py` to see which
   cases flipped — verify that gains aren't paid for by equal losses.

## Catalog location

All scripts read pathways from `$PATHWAY_CATALOG` (default
`~/gitroot/logic-network-generator/output`). Set this env var if the
catalog lives elsewhere.
