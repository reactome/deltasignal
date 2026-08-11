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

The default panel contains all 847 cases with experimental outcomes across the
ten MP-BioPath pathways. The three pathways examined during development are
fixed in `bench/mpbiopath_ten_pathways.tsv`; the other seven form the held-out
test set. Do not change that split after looking at held-out results.

For release-aware failure classification, first export the benchmark-relevant
identifiers from the exact Reactome Neo4j release used to generate the catalog:

```bash
python bench/export_reactome_id_audit.py \
  --supplementary-workbook /path/to/PredictiveAccuracyOfBiologicalPathways_SupplementaryTables.xlsx \
  --id-map /path/to/db_id_to_name_mapping.txt \
  --neo4j-http http://127.0.0.1:7474/db/graph.db/tx/commit \
  --output /path/to/reactome_id_audit.tsv
```

Pass the resulting file with `--reactome-id-audit`. Key-output reaction proxies
are disabled by default because they change the biological readout being
scored. `--allow-output-proxies` enables only proxies explicitly exported by
LNG and records their role in every case row.

Run the frozen-catalog solver factorial with:

```bash
python bench/run_mpbiopath_factorial.py \
  --catalog release97_diagram_on=/path/to/catalog-diagram-on \
  --catalog release97_diagram_off=/path/to/catalog-diagram-off \
  --supplementary-workbook /path/to/PredictiveAccuracyOfBiologicalPathways_SupplementaryTables.xlsx \
  --id-map /path/to/db_id_to_name_mapping.txt \
  --reactome-id-audit /path/to/reactome_id_audit.tsv \
  --output-dir /path/to/factorial-output
```

This crosses each frozen graph catalog with legacy/current propagation and
flat/SCC solving. Stoichiometric weighting is intentionally absent: current
production defaults do not use it, and Adam's perturbation tests found no
accuracy improvement from enabling it.

After completing the current-SCC diagram-on and diagram-off runs, generate the
auditable scorecard and presentation-ready figures with:

```bash
python bench/generate_evaluation_report.py \
  --diagram-on-dir /path/to/v97-diagram-on-current-scc \
  --diagram-off-dir /path/to/v97-diagram-off-current-scc \
  --tcga-readiness /path/to/readiness_key_findings.tsv \
  --output-dir /path/to/evaluation-report
```

The report distinguishes case-level convergence flags from unique perturbation
solves. One failed solve can be repeated across several key-output cases, so
reporting only the case count can exaggerate the number of numerical failures.
It also reports accuracy on converged cases separately; non-converged outputs
remain diagnostic results and must not be presented as equally reliable model
predictions.

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

Use `--allow-output-proxies` to prefer an exact key-output entity and fall back
to LNG's explicit producing/consuming reaction mapping only when that entity is
absent. Proxy fallback must be reported separately from exact-entity scoring.
For an internal validity control, combine `--derive-output-proxies` with
`--prefer-output-proxies`; the harness then derives adjacent reaction readouts
for exported entities so entity/proxy classification agreement can be measured
on overlap cases. This audit mode is not the production endpoint.

The summary also reports no-change, development-class-frequency, signed
reachability, and shortest-signed-path baselines. Class frequency is fitted
only on the development pathways. Structural baselines are cycle-safe and
return unchanged when positive and negative paths conflict.

## How To Read The Two Scoreboards

The frozen-network and current-stack results answer different questions:

| Scoreboard | What is held fixed? | What does it measure? |
| --- | --- | --- |
| Frozen network | Older LNG graph files and the perturbation cases | The effect of changing DeltaSignal propagation or solver behavior |
| Current stack | Only the published perturbation cases | The combined effect of current LNG, Reactome data, mappings, and DeltaSignal |

A DeltaSignal prediction depends on both layers:

```text
prediction = DeltaSignal propagation(LNG-generated network, perturbation)
```

Changing LNG can alter nodes, edges, signs, logic gates, positional instances,
and gene/output mappings. A current-stack gain therefore cannot automatically
be attributed to DeltaSignal's propagation or SCC solver.

The MP-BioPath and curator predictions are read from Supplementary Table S1;
they are not rerun through LNG. The summary reports them at two scopes:

- `*_on_scored_cases`: paired accuracy on exactly the cases DeltaSignal could
  map and score;
- `*_on_all_eligible_cases`: accuracy across the complete eligible panel.

Paired accuracy is useful for comparing classifications without treating a
missing DeltaSignal mapping as a wrong prediction. It must always be shown
beside DeltaSignal coverage and coverage-adjusted accuracy. Otherwise a method
that scores only a favorable subset can appear better than a method that covers
the entire panel.

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

This is a conditional comparison on the 223 cases DeltaSignal could score. It
does **not** show overall DeltaSignal superiority. Across all 344 eligible
cases:

- MP-BioPath was correct for 259/344 (75.3%);
- curator predictions were correct for 282/344 (82.0%);
- DeltaSignal was correct for 164/344 eligible cases (47.7%
  coverage-adjusted accuracy);
- DeltaSignal left 121/344 cases unscored because the key output database ID
  was absent from the generated network.

The paired bootstrap difference against MP-BioPath was only +0.9 percentage
points, with a 95% interval from -3.6 to +4.9 points. This supports approximate
parity on the mapped subset, not a statistically established win. Mean rather
than maximum aggregation across duplicate output UUIDs gave 165/223, so the
paired result was not sensitive to that choice.

The pathway results were uneven:

- PIP3 activates AKT signaling: 82/84;
- Mitotic G1/G1-S: 39/86;
- Cell Cycle Checkpoints: 43/53.

On these same regenerated networks, the legacy propagation configuration
scored 142/223 and reported 49 non-converged cases. Current propagation with
the flat solver scored 166/223 but reported 142 non-converged cases. SCC solving
removed all convergence failures and changed the total from 166 to 164.

For the 204 cases scorable in both the old and regenerated catalogs, current
DeltaSignal improved from 91/204 to 145/204. MP-BioPath was also correct for
145/204 of these same cases. The regenerated catalog made 19 additional cases
scorable; DeltaSignal classified all 19 correctly, while MP-BioPath classified
17 correctly. Those additional cases produce the 164-versus-162 headline.

The complete reconciliation is:

| Case set | Current DS, old LNG | Current DS, new LNG | MP-BioPath |
| --- | ---: | ---: | ---: |
| 204 cases scorable in both catalogs | 91/204 | 145/204 | 145/204 |
| 19 cases newly scorable with current LNG | unscored | 19/19 | 17/19 |
| Current-stack paired total | not applicable | 164/223 | 162/223 |

Among the 204 shared cases, 55 DeltaSignal classifications changed from wrong
to correct and one changed from correct to wrong. The pathway-level changes
were:

| Pathway | Current DS, old LNG | Current DS, new LNG | Net change |
| --- | ---: | ---: | ---: |
| PIP3 activates AKT signaling | 30/74 | 72/74 | +42 |
| Mitotic G1/G1-S | 35/85 | 38/85 | +3 |
| Cell Cycle Checkpoints | 26/45 | 35/45 | +9 |

### Attribution And Implications

| Component | Evidence | Supported interpretation |
| --- | --- | --- |
| DeltaSignal propagation changes | 88/204 legacy to 91/204 current on the same old graphs | Small, pathway-dependent internal improvement |
| SCC solver on old graphs | Removed convergence failures without changing classifications | Better numerical reporting, no accuracy gain on that subset |
| Current LNG/Reactome networks | 91/204 to 145/204 under current DeltaSignal semantics | Primary source of the large shared-case gain |
| Additional mapping coverage | 19 new cases: DeltaSignal 19/19, MP-BioPath 17/19 | Explains the two-case paired headline advantage |
| Overall coverage | DeltaSignal scored 223/344; MP-BioPath predicted all 344 | DeltaSignal is not yet the better complete system |

On the regenerated graphs, legacy propagation scored 142/223 with 49
non-converged cases, while current propagation with SCC scored 164/223 with no
convergence failures. This suggests an interaction between propagation
semantics and the newer topology, but the legacy accuracy is not a clean
production comparator when 49 solves do not converge.

The defensible conclusion is:

1. the integrated current LNG plus DeltaSignal stack improved substantially on
   these three pathways;
2. the large gain is attributable mainly to LNG topology and mapping, especially
   PIP3/AKT, rather than SCC solving or DeltaSignal internals alone;
3. DeltaSignal reaches paired accuracy parity with MP-BioPath where it can
   produce a prediction, but still has substantially worse end-to-end coverage;
4. the PIP3-heavy gain and perfect 19/19 newly mapped result require replication
   on pathways not used during recent LNG/DeltaSignal development.

The regenerated PIP3 graph was also produced twice with a fixed Python hash
seed. UUID-labeled CSV hashes differed because LNG still creates random UUID4
identifiers, but a 20-round attributed graph-refinement check produced the same
semantic graph hash for both runs. Deterministic UUID5 identifiers would make
byte-level provenance and run comparison much cleaner.
