# DeltaSignal Benchmarks

DeltaSignal changes must be evaluated on two separate scoreboards:

1. **Solver-only:** keep LNG network files and perturbation cases fixed while
   changing only the DeltaSignal commit or `DS_*` configuration.
2. **End-to-end:** regenerate networks with a specified LNG commit, then test
   the complete LNG-to-DeltaSignal stack.

TCGA survival analysis is an external-use validation. It can show clinical
association, but it does not by itself measure whether the pathway solver made
more accurate causal predictions.

## Solver Invariants

The synthetic benchmark separates propagation math from feedback-loop solving:

```bash
julia --project=. bench/benchmark_solver_invariants.jl \
  --output-dir /tmp/deltasignal-solver-invariants
```

It checks causal polarity, complex assembly, disconnected-component locality,
edge-order invariance, observation pinning, boundedness, and convergence. These
are necessary model properties, not substitutes for empirical validation.

## MP-BioPath Perturbation Cases

Install the workbook reader:

```bash
python -m pip install -r bench/requirements.txt
```

Run the published perturbation cases against a frozen LNG catalog:

```bash
python bench/benchmark_mpbiopath_cases.py \
  --supplementary-workbook /path/to/PredictiveAccuracyOfBiologicalPathways_SupplementaryTables.xlsx \
  --id-map /path/to/db_id_to_name_mapping.txt \
  --catalog /path/to/logic-network-generator/output \
  --ground-truth experimental \
  --output-dir /path/to/benchmark-output
```

The harness starts an isolated DeltaSignal API process, hashes all inputs and
network files, caches one solve per perturbation, and writes:

- `benchmark_cases.tsv`: one auditable row per perturbation/readout case;
- `benchmark_summary.json`: coverage, macro-F1, balanced accuracy, per-class
  metrics, convergence, and paired comparisons;
- `benchmark_manifest.json`: Git state, configuration overrides, thresholds,
  input hashes, and output locations;
- `deltasignal_server.log`: the complete API execution log.

Unmapped cases are reported as unscored with an exact reason. They are never
silently converted to unchanged predictions. Use a development split to choose
thresholds and configuration, then report the final result once on held-out
empirical pathways.
