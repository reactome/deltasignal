# Quickstart: reproducing every figure in this feature

Every number in `spec.md` and `research.md` comes from the commands below.
All arms in a comparison **must** run against one catalog build — see the
A/B protocol note at the end for why.

## Prerequisites

- Neo4j holding Reactome **Release97**, reachable at `NEO4J_URL` (dev
  credentials are in the logic-network-generator `.env`).
- Julia 1.10 on `PATH`: `export PATH=$HOME/julia-1.10.10/bin:$PATH`
- The MP-BioPath comparison inputs:
  - `~/codes_and_results/PredictiveAccuracyOfBiologicalPathways_SupplementaryTables.xlsx`
  - `~/codes_and_results/db_id_to_name_mapping.txt`
- `PYTHONHASHSEED=0` exported **before** the interpreter starts, or catalog
  builds are not reproducible.

Set a scratch root once; the rest of this document uses it:

```bash
export S=/tmp/claude-1001/.../scratchpad     # session scratchpad
export W=~/codes_and_results
export PATH=$HOME/julia-1.10.10/bin:$PATH
export PYTHONHASHSEED=0
```

## 1. Build the shared catalog (once)

Ten pathways, logic-network-generator `main` at 4ff0408:

`create-pathways.py` takes **one** stable id per run (`--pathway-id`) or a
TSV list (`--pathway-list`, columns `id` and `pathway_name`). The list form
is the one to use here:

```bash
cd ~/gitroot/logic-network-generator
printf 'id\tpathway_name\n' > $S/pathways.tsv
for p in R-HSA-1227986 R-HSA-1257604 R-HSA-195721 R-HSA-3700989 R-HSA-453279 \
         R-HSA-5673001 R-HSA-5693567 R-HSA-68875 R-HSA-69242 R-HSA-69620; do
  printf '%s\t%s\n' "$p" "$p" >> $S/pathways.tsv
done
python3 bin/create-pathways.py --pathway-list $S/pathways.tsv --output-dir $S/cat
ls $S/cat | wc -l    # expect 10
```

The directory names the benchmark expects embed the pathway name, so put the
real Reactome names in the second column rather than repeating the id if you
are building the catalog from scratch.

`create-pathways.py` re-execs itself with `PYTHONHASHSEED=0` if it is not
already set, so the export above is belt-and-braces rather than required.

**A stale-cache trap:** regenerating into a directory that already holds
pre-fingerprint caches is a silent no-op — `_cache_is_reusable` adopts them
and reproduces the old networks while reporting success. `rm -rf
<dir>/*/cache` before any regeneration you intend to be real.

## 2. Pin the AND behaviour (T002–T004)

```bash
cd ~/gitroot/deltasignal
julia --project=. test/test_and_curves.jl
```

All five testsets pass on the current defaults. To confirm the test can
actually fail — the point of writing it before the change:

```bash
DS_AND_MODE=hill_log DS_HILL_SAT_EPS=0.001 DS_ASSEMBLY_LIMITING=1 \
  julia --project=. test/test_and_curves.jl
# expect: the three curve testsets still pass, and
#         "defaults implement the design intent" fails 3 of 3
```

## 3. Experimental benchmark, both arms (T009, T013)

New defaults — no env overrides, the code defaults *are* the arm:

```bash
python3 bench/benchmark_mpbiopath_cases.py \
  --supplementary-workbook $W/PredictiveAccuracyOfBiologicalPathways_SupplementaryTables.xlsx \
  --id-map $W/db_id_to_name_mapping.txt \
  --catalog $S/cat --output-dir $S/i_defaults --port 8301
```

Old defaults, for the A/B:

```bash
DS_AND_MODE=hill_log DS_HILL_SAT_EPS=0.001 DS_ASSEMBLY_LIMITING=1 \
python3 bench/benchmark_mpbiopath_cases.py \
  --supplementary-workbook $W/PredictiveAccuracyOfBiologicalPathways_SupplementaryTables.xlsx \
  --id-map $W/db_id_to_name_mapping.txt \
  --catalog $S/cat --output-dir $S/a_clampON --port 8302
```

Expected: 365/564 (macro-F1 0.5781) and 334/564 (0.5601).

## 4. Curator benchmark, both arms (T011)

Identical to step 3 with `--ground-truth curator` and fresh ports and output
directories. Expected: **2628/3914** (macro-F1 0.6408) new, **2479/3914**
(0.5798) old.

## 5. Score any pair of arms

`benchmark_summary.json` reports accuracy; macro-F1, per-class recall,
per-pathway attribution and the both-arms-converged count come from
`benchmark_cases.tsv`. The columns that matter are `expected`,
`prediction` (`0`=DOWN, `1`=NO_CHANGE, `2`=UP; empty = unscored),
`converged`, `pathway_name`, and the baseline columns
`mpbiopath_prediction` and `shortest_signed_path_prediction`.

Restrict every baseline to the rows DeltaSignal actually scores — MP-BioPath
has a prediction on 847 rows but only 564 are comparable, and quoting 625
instead of 407 against DeltaSignal's 365 compares different case sets.

## Mandatory A/B protocol

uuid4 node ids are minted fresh per catalog build, `Dict` order follows them,
and that sets Gauss-Seidel sweep order inside an SCC — so **a non-converged
solve returns different values for a structurally identical network**. With
179 of 564 experimental cases non-converged under the new defaults this is
not a hypothetical. Before believing any result:

1. Both arms must read the *same* catalog directory.
2. Report per-pathway net change. A gain in a pathway the change did not
   touch is the tell.
3. Count changed predictions where **both** arms converged. For this
   feature: 59 of 86 (experimental) and 282 of 462 (curator).
