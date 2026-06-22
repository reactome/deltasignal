# TCGA LUAD DeltaSignal Validation

This folder contains a reproducible validation pipeline for testing DeltaSignal on TCGA lung adenocarcinoma (LUAD) RNA-seq and survival metadata.

The pipeline connects:

1. Logic Network Generator pathway outputs (`logic_network.csv`, `stid_to_uuid_mapping.csv`)
2. DeltaSignal parsing and steady-state solving
3. TCGA LUAD RNA expression observations
4. Survival analysis and shuffled-label controls
5. Publication-style figures and a plain-English result interpretation report

## Focused Validation Panel

The current focused panel is:

- `R-HSA-69620`: Cell Cycle Checkpoints
- `R-HSA-453279`: Mitotic G1-G1/S phases
- `R-HSA-1257604`: PIP3 activates AKT signaling

These were selected after the mRNA observability audit and LNG mapping fix.

## Inputs

Required inputs are intentionally explicit:

- TCGA LUAD expression TSV, genes by samples
- TCGA LUAD clinical matrix with survival fields
- LNG output directory containing one folder per pathway
- Observability audit tables:
  - `pathway_suitability.tsv`
  - `root_observability_audit.tsv`

The local defaults assume the surrounding GSoC workspace layout:

```bash
/Users/chryseis/gsoc/
  deltasignal/
  logic-network-generator/output/
  deltasignal-pipeline/data/
  real-deltasignal-pipeline/outputs/observability_audit_panel_rhsa_after_lng_mapping_fix/
```

For another machine, pass the input paths explicitly.

## Run DeltaSignal Scoring

From the DeltaSignal repo root:

```bash
python validation/tcga_luad/scripts/run_validation.py \
  --sample-count 0 \
  --shuffle-count 1000 \
  --solve-mode batch \
  --max-iters 2000 \
  --output-dir validation/tcga_luad/outputs/tcga_luad_three_pathways
```

`--sample-count 0` means score all eligible samples.

The batch solver starts one Julia process per pathway and solves every sample for that pathway in one process. This avoids the large per-sample Julia startup overhead seen in early readiness runs.

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

## Main Current Result

In the full eligible TCGA LUAD focused run, Cell Cycle Checkpoints (`R-HSA-69620`) was the strongest signal:

- 502 samples scored
- 0 DeltaSignal solve failures
- 502/502 solves converged for the pathway
- terminal mean log-rank p = `8.92e-05`
- 0/1000 shuffled controls were as strong or stronger

Interpret this as a promising association, not causal proof. The result means that, under the current mRNA-to-DeltaSignal pipeline, predicted Cell Cycle Checkpoints activity separates TCGA LUAD patients into groups with meaningfully different survival.

## Important Caveats

- This focused run tested only three pathways.
- Median-split survival is a useful sanity check but not the final clinical model.
- Reactome-wide expansion needs multiple-testing correction across pathways and metrics.
- The current expression-to-observation transform is cohort percentile activity.
- Terminal outputs are graph-derived terminal nodes, not manually curated pathway readouts.
- Batch scoring currently skips influence scores because the validation summaries use node activities and solver diagnostics only.
