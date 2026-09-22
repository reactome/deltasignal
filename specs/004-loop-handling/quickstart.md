# Quickstart: reproducing feature 004's figures

Prerequisites are the same as
`specs/002-upregulation-propagation/quickstart.md` (Julia on `PATH`, the
MP-BioPath workbook and id map in `~/codes_and_results`, `PYTHONHASHSEED=0`,
a Release97 Neo4j only if you re-derive the identifier check).

```bash
export S=/tmp/claude-1001/.../scratchpad
export W=~/codes_and_results
export PATH=$HOME/julia-1.10.10/bin:$PATH
export MPB=~/gitroot/mp-biopath-pathways
```

## 1. The acyclicity asymmetry (R1, FR-010)

```bash
python3 bench/analysis/cycle_structure.py --catalog $S/cat
python3 bench/analysis/cycle_structure.py --mpbiopath $MPB/pathways
```

Expected: LNG 4,561 cycle-resident nodes across 34 components, largest 836
(TP53), 643 (Cell Cycle Checkpoints), 465 (ERBB2), 443 (RAF). MP-BioPath:
zero self-loops, four of nine fully acyclic, largest component 11 nodes.

## 2. Component classification (R3)

```bash
python3 bench/analysis/cycle_structure.py --catalog $S/cat --classify --ratio-threshold 15
```

Expected: six components classified `recycling_artifact`, holding 3,339
nodes, with nodes-per-reaction 20.9–160.8; the rest 3.0–9.3. Re-run with
`--ratio-threshold 10` and `--ratio-threshold 20` — the partition must not
change, and the report states which threshold produced it.

## 3. The control — DeltaSignal on MP-BioPath's networks (US1)

```bash
python3 bench/analysis/mpbiopath_network_adapter.py \
  --pathways $MPB/pathways --out $S/mpb_as_ds
python3 bench/benchmark_mpbiopath_cases.py \
  --supplementary-workbook $W/PredictiveAccuracyOfBiologicalPathways_SupplementaryTables.xlsx \
  --id-map $W/db_id_to_name_mapping.txt \
  --catalog $S/mpb_as_ds --output-dir $S/k_control --port 8321
```

**Verify before trusting it**: the adapter must reproduce MP-BioPath's own
published predictions when its scoring rule is applied, on a spot-check of
cases. If it does not, the control measures the adapter, not the propagator.

Read the result against `$S/i_defaults` (DeltaSignal on its own networks,
365/564) and the `mpbiopath_prediction` column (407/564) restricted to the
common case set.

## 4. Intervention arms (US3, US4)

```bash
python3 bench/analysis/loop_interventions.py --catalog $S/cat \
  --arm drop_diagram_bridge --out $S/cat_nodb
```

Then benchmark exactly as in feature 002's quickstart, one arm at a time, on
a fresh port, killing the Julia server between arms. Arms to run, in order:
`drop_diagram_bridge`, `break_recycling`, both, and `dagify` (upper bound
only).

## 5. Score an arm

Use the same scoring as feature 002 plus the two additions this feature
requires:

- **Coverage delta first.** A case scoreable in the baseline and not in the
  arm is a coverage loss, not an accuracy gain. Read this before the score.
- **Cyclic/acyclic split.** Report DeltaSignal and MP-BioPath separately on
  cases whose readout is cycle-resident and on those whose is not.

## A/B protocol

Unchanged and mandatory: same catalog directory for both arms, per-pathway
net breakdown, and the count of changed predictions where **both** arms
converged. 179 of 564 cases do not converge and their values depend on
sweep order, so a handful of changed cases concentrated in TP53 with no
both-converged support is the known artifact, not a result.
