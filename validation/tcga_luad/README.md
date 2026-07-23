# TCGA LUAD DeltaSignal Validation

This folder contains a reproducible validation pipeline for testing DeltaSignal on TCGA lung adenocarcinoma (LUAD) RNA-seq and survival metadata.

The pipeline connects:

1. Logic Network Generator pathway outputs (`logic_network.csv`, `stid_to_uuid_mapping.csv`)
2. DeltaSignal parsing and steady-state solving
3. TCGA LUAD RNA expression observations
4. Survival analysis and shuffled-label controls
5. Publication-style figures and a plain-English result interpretation report

## Example Validation Panel

selected examples:

- `R-HSA-69620`: Cell Cycle Checkpoints
- `R-HSA-453279`: Mitotic G1-G1/S phases
- `R-HSA-1257604`: PIP3 activates AKT signaling

These were selected after the mRNA observability audit and LNG mapping fix.

## Inputs

- TCGA LUAD expression TSV, genes by samples
- TCGA LUAD clinical matrix with survival fields
- LNG output directory containing one folder per pathway
- Observability audit tables:
  - `pathway_suitability.tsv`
  - `root_observability_audit.tsv`

### 1. Clone The Logic Network Generator

```bash
git clone git@github.com:reactome/logic-network-generator.git
```

Use Logic Network Generator to produce one output folder per Reactome pathway. Each pathway output folder must contain:

```text
logic_network.csv
stid_to_uuid_mapping.csv
```

Pass the parent output directory to this pipeline with `--lng-output-root`.

### 2. Prepare TCGA LUAD Data

Prepare two tab-separated files:

Expression matrix:

```text
gene    TCGA-...-01    TCGA-...-01    ...
TP53    12.3           8.7            ...
...
```

Clinical matrix:

```text
sampleID        days_to_death    days_to_last_followup    vital_status
TCGA-...-01     540              NA                       DECEASED
TCGA-...-01     NA               1200                     LIVING
...
```

The runner currently expects TCGA tumor sample IDs ending in `-01`.

### 3. Prepare Observability Audit Tables

The scoring runner uses a prior observability audit to decide which root nodes can be backed by mRNA measurements.

Required `pathway_suitability.tsv` columns:

```text
pathway_id
pathway_name
suitability
```

Required `root_observability_audit.tsv` columns:

```text
pathway_id
root_uuid
mapping_status
gene_symbol
```

Rows with `mapping_status == "mapped"` and a gene present in the expression matrix become DeltaSignal root observations.

## Run DeltaSignal Scoring

From the DeltaSignal repo root:

```bash
python validation/tcga_luad/scripts/run_validation.py \
  --lng-output-root /path/to/logic-network-generator/output \
  --observability-summary /path/to/pathway_suitability.tsv \
  --root-audit /path/to/root_observability_audit.tsv \
  --expression-tsv /path/to/TCGA.LUAD.expression.tsv \
  --clinical-tsv /path/to/TCGA.LUAD.clinicalMatrix.tsv \
  --sample-count 0 \
  --shuffle-count 1000 \
  --solve-mode batch \
  --max-iters 2000 \
  --output-dir validation/tcga_luad/outputs/tcga_luad_three_pathways
```

`--sample-count 0` means score all eligible samples.

The batch solver starts one Julia process per pathway and solves every sample for that pathway in one process. This avoids the large per-sample Julia startup overhead seen in early readiness runs.

Each result JSON records the DeltaSignal commit and effective `DS_*` solver
configuration. `readiness_run_summary.json` also records DeltaSignal and LNG
Git state plus SHA-256 hashes for the expression, clinical, audit, LNG, and
parsed-network inputs.

To compare two runs on the same samples and pathways:

```bash
python validation/tcga_luad/scripts/compare_validation_runs.py \
  --baseline-dir /path/to/baseline-run \
  --candidate-dir /path/to/candidate-run \
  --output-dir /path/to/comparison-output
```

The comparison reports score correlations, absolute shifts, boundary
saturation, top-quartile overlap, and rank-quartile agreement. A new full-cohort
run should proceed only after a small current-stack smoke run passes these
distribution and convergence checks.

## Generate Figures and Report

Use a Python environment with `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`, and `lifelines`.

```bash
python -m pip install -r validation/tcga_luad/requirements.txt
```

```bash
python validation/tcga_luad/scripts/generate_figures_and_report.py \
  --run-dir validation/tcga_luad/outputs/tcga_luad_three_pathways \
  --output-dir validation/tcga_luad/example_results/tcga_luad_three_pathways \
  --shuffle-count 1000 \
  --copy-summary-tables
```

Generated artifacts include:

- Kaplan-Meier survival curves for each pathway
- Box/violin plots of pathway activity by survival event
- Histograms/density plots showing the median split
- Samples x pathways activity heatmap
- Pair plot/scatter matrix of the three pathway scores
- Shuffled-label null plots
- Univariate Cox model forest plot when `lifelines` is installed
- `reports/RESULT_INTERPRETATION.md`

## Historical Demonstration Result

The following result was produced with an older LNG network catalog and older
DeltaSignal solver defaults. It demonstrates that the workflow runs end to end;
it must not be presented as a current-stack benchmark.

In that example run, Cell Cycle Checkpoints (the Reactome pathway `R-HSA-69620`) was the strongest signal:

this pathway is about mechanisms that control whether cells are allowed to keep dividing: DNA damage checkpoints, G1/S transition control, G2/M control, mitotic checkpoint behavior, etc. In cancer, high activity in this kind of pathway often means the tumor is highly proliferative or under replication/checkpoint stress. That can be associated with more aggressive disease.

For this LUAD run:
We scored 502 TCGA lung adenocarcinoma tumor samples.
For each sample, DeltaSignal predicted activity across the Cell Cycle Checkpoints logic network.
We summarized each sample by terminal_mean_activity, meaning: average predicted activity of terminal/output nodes in that pathway network.
Then we split patients into two equal groups:251 patients with lower Cell Cycle Checkpoints activity, 251 patients with higher Cell Cycle Checkpoints activity.

The survival result says:
Low activity group: 73 deaths, KM median survival 1725 days
High activity group: 109 deaths, KM median survival 1229 days
(patients whose tumors had higher predicted Cell Cycle Checkpoints pathway activity died more often and had shorter median survival in this cohort.)
The log-rank p-value: p = 8.92e-05
The shuffle control: 0/1000 shuffled controls were as strong or stronger
(I randomly scrambled the pathway activity labels across patients 1000 times. None of those random fake splits produced a log-rank p-value as small as the real Cell Cycle Checkpoints split.)

This is an association in one cohort, not evidence that DeltaSignal predicts
causal pathway responses accurately. The result must be regenerated after any
LNG graph or DeltaSignal propagation change before it is compared with a newer
run.


## Important Caveats

- This focused run tested only three pathways.
- Median-split survival is a useful sanity check but not the final clinical model.
- Reactome-wide expansion needs multiple-testing correction across pathways and metrics.
- The current expression-to-observation transform is cohort percentile activity.
- Terminal outputs are graph-derived terminal nodes, not manually curated pathway readouts.
- Batch scoring currently skips influence scores because the validation summaries use node activities and solver diagnostics only.
- The historical example predates the current LNG and DeltaSignal defaults.
- Survival association is not a solver-accuracy benchmark; use the fixed
  perturbation benchmarks under `bench/` to assess DeltaSignal itself.
