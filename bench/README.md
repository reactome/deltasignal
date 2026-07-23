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

## Preliminary Frozen-Network Result

The first comparison used the three pathways already present in an older LNG
catalog. Of 344 experimental cases, 204 were scorable after the ID and network
mapping audit. On those same 204 cases:

- legacy DeltaSignal was correct for 88 cases;
- current DeltaSignal was correct for 91 cases;
- MP-BioPath was correct for 145 cases;
- curator predictions were correct for 156 cases.

This is a small improvement for current DeltaSignal, but it remains well behind
MP-BioPath on this fixed subset. The change is not consistent across pathways:
PIP3/AKT improves, Mitotic G1/G1-S regresses, and Cell Cycle Checkpoints is
unchanged. Enabling SCC solving removes the current solver's non-convergence
reports but does not change any discrete classifications here, so the three-case
gain comes from propagation semantics rather than loop solving.

These are diagnostic results on older frozen LNG graphs, not the final
current-stack benchmark. Regenerated current-LNG networks must be reported on
the separate end-to-end scoreboard described above.

## Preliminary Current-Stack Result

A second run regenerated all three networks from Reactome Release 96 with LNG
`09b1597`, then tested DeltaSignal `ff39e3b`. Of 344 eligible experimental
cases, 223 were scorable:

- current DeltaSignal with SCC solving: 164/223 (73.5%);
- MP-BioPath on the same cases: 162/223 (72.6%);
- curator predictions on the same cases: 175/223 (78.5%).

The paired bootstrap difference against MP-BioPath was +0.9 percentage points
with a 95% interval from -3.6 to +4.9 points. This supports approximate parity
on this subset, not a statistically established win. Mean rather than maximum
aggregation across duplicate output UUIDs gave 165/223, so the overall result
was not sensitive to that choice.

The pathway results were uneven:

- PIP3 activates AKT signaling: 82/84;
- Mitotic G1/G1-S: 39/86;
- Cell Cycle Checkpoints: 43/53.

On these same regenerated networks, the legacy propagation configuration
scored 142/223 and reported 49 non-converged cases. Current propagation with
the flat solver scored 166/223 but reported 142 non-converged cases. SCC solving
removed all convergence failures and changed the total from 166 to 164.

For the 204 cases scorable in both the old and regenerated catalogs, current
DeltaSignal improved from 91/204 to 145/204. Nineteen additional cases became
scorable in the regenerated catalog and all 19 were classified correctly.
This large shift is driven mainly by PIP3/AKT and must be replicated on held-out
pathways before it is treated as general predictive improvement.

The regenerated PIP3 graph was also produced twice with a fixed Python hash
seed. UUID-labeled CSV hashes differed because LNG still creates random UUID4
identifiers, but a 20-round attributed graph-refinement check produced the same
semantic graph hash for both runs. Deterministic UUID5 identifiers would make
byte-level provenance and run comparison much cleaner.
