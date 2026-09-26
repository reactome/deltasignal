# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Quick Start (Docker-First)

This repo is the DeltaSignal **engine + HTTP API**. The UI lives in the
WebsiteAngular workspace and talks to this API over HTTP — see `docs/API.md`
for the contract.

```bash
# Development: start the Julia API (binds 0.0.0.0:8080, published on the host)
docker compose -f docker-compose.dev.yml up julia-api
# Then point the Angular dev server's proxy at http://localhost:8080.

# Production: single API service (TLS / reverse proxy / UI hosting are handled
# by the WebsiteAngular deployment, not here)
docker compose -f docker-compose.prod.yml up

# Health check:
curl http://localhost:8080/api/health
```

## Catalogs: what is being served, and what produced it

Generated networks live on **real disk** under `~/deltasignal-catalogs/`,
never in `/tmp`. On this machine `/tmp` is a tmpfs; a reboot on 2026-09-24
wiped every catalog behind specs/018–021, which still cite them by name. Nothing
recorded what the dev API was serving either, and it quietly served a July
catalog for two months while every benchmark ran against newer ones.

```bash
scripts/catalog.sh build    # generate from generator HEAD, verify, move `current`
scripts/catalog.sh status   # what `current` is, and whether the generator has moved on
scripts/catalog.sh list     # every build and which one is current
scripts/catalog.sh bench    # benchmark `current` on both axes, filed under the build
curl localhost:8080/api/health | jq .catalog    # what the RUNNING API is serving
```

- Each build is `builds/<date>_<generator-sha>/` with a `BUILD.json` recording
  the generator commit, environment, the exact pathway list given, and whether
  every pathway built. `current` only moves after a build verifies complete, so
  a partial build can never become what the API serves.
- The pathway list is versioned in `bench/catalog_pathways.tsv`.
- `docker-compose.dev.yml` mounts `~/deltasignal-catalogs/current`. **The link
  is resolved every time the container STARTS**, not when it is created — so a
  plain restart, a crash-restart under `restart: unless-stopped`, or a reboot
  picks up a newly moved `current`. After a build, recreate the API to switch
  deliberately: `docker compose -f docker-compose.dev.yml up -d --force-recreate julia-api`.
- The mount uses `create_host_path: false`, so a missing catalog **fails the
  container start** instead of docker silently creating an empty directory.
- `bench` resolves the build once, refuses unless `/api/health` reports it, and
  re-checks after each axis in case the container restarted mid-run. It also
  refuses if `src/` or `bench/` is uncommitted, and records nothing if a metric
  cannot be parsed. Results go to `builds/<id>/results/<solver-sha>/`, so every
  number is traceable to a catalog build and a solver commit. **Cite that pair
  in a spec, not a scratch path.**
- `BUILD.json` records the Reactome release, release date and content checksum
  from the database's own `DBInfo` node, because the same generator commit
  against a different release builds a different catalog.

## Benchmark protocol: what gets pinned

**Protocol of record (Adam, 2026-09-25): perturb only ROOT inputs that are, or
contain, the gene**, as the MP-BioPath publication does. The perturbation then
propagates from root input to terminal readout; if the ends are right, the
middle is taken to be right. `DS_PIN_SCOPE` (benchmark-side) selects it:

- `root_cycle` (**default**, specs/024): `root`, plus, for a gene with no
  true root, the form a catalytic cycle regenerates (fed only by its own
  downstream or by dissociation/depletion edges). This is how MP-BioPath's
  hand-cut loops left enzymes as roots, reproduced while keeping the loops.
- `root`: nodes with no incoming edge that are or contain the gene. A gene
  with no root form is not perturbed, and its cases are invalid.
- `entry`: the resolved nodes no other resolved node reaches. It also pins a
  gene's first mid-pathway occurrence. Measured, not adopted.
- `all`: every node whose `member_leaves` include the gene. This was the
  protocol from 2026-07-14 (`0bd4565`) to 2026-09-25: 9,378 pins, 89% of them
  mid-pathway complexes, which were *set* rather than computed. It is kept only
  to reproduce numbers from that period, which are **not comparable** to root
  numbers (specs/023).

Every benchmark log prints `Protocol:` and `Pinned:` lines. Check them instead
of inferring the protocol from flags.

**Run arms with `scripts/run_arm.sh NAME [--server DS_X=v] [--bench DS_Y=v]`,
not ad-hoc scripts.** It pins a detached worktree at HEAD, serves it with the
dev compose environment plus the named overrides on the resolved build, and
refuses to run in three cases:
- the container does not show an override;
- the pinned code never reads it;
- `src/` or `bench/` is dirty.

It writes `ARM.json` next to the results.

## Development Commands

### Running Tests
```bash
# Run all tests in Docker (recommended)
docker-compose -f docker-compose.dev.yml --profile test up test-runner

# Run specific tests
docker-compose -f docker-compose.dev.yml run --rm julia-api julia --project=/app /app/test/test_basic.jl
docker-compose -f docker-compose.dev.yml run --rm julia-api julia --project=/app /app/test/test_steady_state.jl

# Local Julia (if available)
julia test/test_basic.jl
```

### Building and Dependencies
```bash
# Build development environment (recommended)
docker-compose -f docker-compose.dev.yml build

# Build production image
./scripts/build.sh

# Start development with auto-rebuild
docker-compose -f docker-compose.dev.yml up --build

# Local Julia setup (if needed)
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### CLI Usage
```bash
# Available commands: parse, solve, export
# Use --help with any command for options

# 1. Parse logic network from TSV files
docker compose -f docker-compose.dev.yml run --rm julia-api sh -c \
  'mkdir -p /app/data/output && julia --project=/app /app/cli/deltasignal.jl parse \
  --logic /app/examples/sample_logic_network.tsv \
  --uuid-map /app/examples/sample_uuid_mapping.tsv \
  --set-map /app/examples/sample_set_mappings.tsv \
  --output /app/data/output/network.json \
  --validate'

# 2. Solve steady-state network
docker compose -f docker-compose.dev.yml run --rm julia-api julia --project=/app /app/cli/deltasignal.jl solve \
  --network /app/data/output/network.json \
  --observations /app/examples/sample_observations.csv \
  --output /app/data/output/results.json
# --mu / --gamma exist but are read ONLY by DS_SCC_METHOD=minimize. Under the
# default they are inert and reported as `nothing`. Do not copy `--gamma 0.1`:
# it is measured harmful under the minimiser (a 100x perturbation reads 51x).

# 3. Export results in various formats
docker compose -f docker-compose.dev.yml run --rm julia-api julia --project=/app /app/cli/deltasignal.jl export \
  --results /app/data/output/results.json \
  --format csv \
  --output /app/data/output/export.csv

# Export formats: pathway-browser, cytoscape, csv
# Aggregation methods: stoichiometry_weighted, mean, max, min, confidence_weighted, geometric_mean

# Complete workflow example (parse → solve → export):
docker compose -f docker-compose.dev.yml run --rm julia-api sh -c '
  mkdir -p /app/data/output &&
  julia --project=/app /app/cli/deltasignal.jl parse \
    --logic /app/examples/sample_logic_network.tsv \
    --uuid-map /app/examples/sample_uuid_mapping.tsv \
    --set-map /app/examples/sample_set_mappings.tsv \
    --output /app/data/output/network.json &&
  julia --project=/app /app/cli/deltasignal.jl solve \
    --network /app/data/output/network.json \
    --observations /app/examples/sample_observations.csv \
    --output /app/data/output/results.json &&
  julia --project=/app /app/cli/deltasignal.jl export \
    --results /app/data/output/results.json \
    --format csv \
    --output /app/data/output/export.csv'

# Interactive Julia shell (for development)
docker compose -f docker-compose.dev.yml --profile cli run --rm cli
```

## Architecture Overview

### Core Mathematical Pipeline
The system implements a biologically-realistic reaction model with the following mathematical transformations:

1. **Sensitivity Transform** (`src/core/sensitivity.jl`): Adaptive input response based on activity levels
2. **Multi-Input Aggregation** (`src/core/aggregators.jl`): Combines activator inputs using geometric mean or other methods
3. **Hill Functions** (`src/core/hill_functions.jl`): Implements sigmoidal activation and inhibition
4. **Reaction Model** (`src/core/reaction_model.jl`): Combines all transformations into complete reaction dynamics

### Module Organization

**Core Components** (`src/core/`):
- `reaction_model.jl`: The complete reaction model — `compute_reaction_output_vec` (the live, config-driven propagator) and the AND/OR/inhibition/assembly aggregation modes.
- `aggregators.jl`, `hill_functions.jl`: aggregation primitives used by the reaction model.
- `sensitivity.jl`: `apply_sensitivity_transform`, applied per-activator in the reaction model.

**Solvers** (`src/solvers/`):
- `steady_state.jl`: the live solver — SCC-condensation feed-forward with observation pinning (`solve_steady_state` → `solve_steady_state_penalty` → `solve_scc_ordered!`), plus `compute_influence_scores` for explainability.

> **Quarantined (`attic/`):** the "biological realism" module cluster
> (biological_realism, compartmentalization, experimental_constraints,
> temporal_dynamics, stochastic_effects, feedback_enhancements), the
> `enhanced_steady_state` / `time_dynamic` solvers, and `parameter_learning`
> were moved to `attic/` — a call-graph audit confirmed none are reachable from
> the live parse/solve path. See `attic/README.md` to revive any of them.

**I/O** (`src/io/`):
- `tsv_parser.jl`: Parses logic networks from TSV format with UUID mapping
- `reactome_mapper.jl`: Maps between expanded networks and Reactome pathways

### Data Flow
1. TSV logic network → Parser → ReactionNetwork JSON
2. ReactionNetwork + Observations → Solver → Node activities + influence scores
3. Results can be aggregated back to pathway-level view using Reactome mappings

## Key Data Structures

### Input Formats
- **Logic Network TSV**: `Parent_UUID | Child_UUID | AND/OR | Pos/Neg | Stoichiometry`
- **UUID Mapping TSV**: `Network_UUID | Reactome_ID | Entity_Type | Set_ID`
- **Set Mappings TSV**: `Set_ID | Original_Name | Member_UUIDs`
- **Observations CSV**: `node_uuid,activity,confidence`. The node column may
  also be `node`, and the activity column may also be `value`
  (`examples/sample_observations.csv` uses `node,condition,value,confidence`).
  Activity is on the **0-100 scale** (0 = none, 1 = normal baseline, 100 = 100x)
  and is range-checked; a value outside 0-100 is an error, not a clamp.
  A node listed twice with **different** values is an error, not a
  last-write-wins merge — two contradictory measurements are a conflict, not a
  set. An exactly repeated row is fine.
  Confidence is 0-1, but it is a **binary gate, not a weight**: anything above
  `OBS_CONFIDENCE_TOL` (1e-6) pins exactly as hard as 1.0, and anything at or
  below is discarded. A soft weighted pin is a design direction, not current
  behaviour. The gate is pinned by `test/test_observation_pinning.jl`; it used
  to be inoperative for root nodes, which is the case a perturbed gene almost
  always is. Any extra column — including `condition` — is ignored, so
  a file holding several conditions is applied as one simultaneous
  perturbation set.

### Internal Scales
- **User Interface (I/O files)**: 0-100 scale for node activities
  - **0** = no activity
  - **1** = normal baseline activity (1× normal)
  - **100** = one hundred times normal activity (100× normal)
  - All values constrained to [0, 100] range
- **Internal computation**: 0-1 normalized scale (divide by 100)
  - 0 → 0.0 (no activity)
  - 1 → 0.01 (normal baseline)
  - 100 → 1.0 (100× normal)
- Conversion: `internal = ui_value / 100.0` and `ui_value = internal * 100.0`
- All conversions handled automatically by solvers and CLI

## Important Implementation Notes

### No pathway-specific parameter tuning
Earlier versions of this file described per-context parameter tuning
(signalling vs metabolic vs transcriptional, bistable cell-cycle switches).
**That code is in `attic/` and is not reachable from the live parse/solve
path** — a call-graph audit confirmed it. Every solve uses the same
`ReactionEvalConfig`. Per-edge-type parameters are a design direction, not a
current behaviour.

### Solver Configuration
The reaction model is driven by `DS_*` environment variables, resolved once per
solve into a `ReactionEvalConfig` (see `resolve_reaction_eval_config` in
`reaction_model.jl`). The **code defaults are the validated winning config** —
env vars only override for benchmark sweeps:
- `DS_INHIBITION_MODE=divide`, `DS_AND_MODE=hill_sat`, `DS_OR_MODE=mean`,
  `DS_ASSEMBLY_LIMITING=0` (**since 2026-09-25, specs/023**: under root pinning
  the min() blocked every single-subunit overexpression; off is held-out +401),
  `DS_SELF_INHIBITOR_WEIGHT=0.1` (specs/022+023, below),
  `DS_INHIBITOR_EPS=1e-12`, `DS_HILL_SAT_EPS=1e-9`,
  `DS_DEPLETION_H_MIN=0.1` (= 1/`DS_DEPLETION_H_MAX`, so depletion may suppress
  at most as hard as it may de-repress — it was previously floored at ZERO and
  could suppress without limit; +28 held-out at p<0.0001, and it costs PIP3
  93.75%→90.63%. `specs/011` has the record),
  `DS_HILL_LOG_ZMAX=10.0`. Behaviour is pinned by `test/test_and_curves.jl`.
  Re-measured 2026-09-16 on 23,022 wide-curator cases, one variable at a time:
  `specs/009-solver-defaults/research.md`. That supersedes the `hill_sat` /
  `assembly_limiting=0` choice in `specs/002`, which was made on 564
  experimental cases and cost 286 on the curator set.
  **`docker-compose.dev.yml` must be kept in step with these.** It drifted for
  six days — production runs the code defaults, every benchmark runs the dev
  container, and the tests pin the code defaults, so the divergence made the
  suite fail in the project's own container while measuring a config production
  never ran.
- `DS_AND_MODE=hill_sat` implements the stated AND intent exactly: AND
  multiplies fold-changes capped at 100, so 0.5*0.5 = 0.25, 0.1*0.1 = 0.01,
  `0 x anything` = 0, 10*10 = 100 and 100*100 = 100, all to 0.0%. It costs
  **16 held-out cases (0.085pp, macro-F1 -0.0014)** against `hill_log` on
  23,908 conditioned cases, and buys depth-invariant magnitudes: `hill_log`
  compresses so hard along a cascade that a 100x source reads 74x at one hop
  and **19x at ten**, so two readouts with identical biology and different
  path lengths get different predicted folds. `specs/010` has the full
  justification and the numbers.
  The earlier record that `hill_sat` "scores 90 cases worse" is **void** — it
  was measured against a `hill_sat` that could not represent a knockout and
  inverted its own saturation on wide reactions. Why the remaining 16-case
  downward compression helps is still unexplained.
- SCC-condensation solve is on by default (`DS_SCC_SOLVE=1`); the legacy flat
  iteration and the `"fixed_point"` `SteadyStateParams.method` are not used by
  the CLI or API (both use the penalty/SCC path).
- `DS_COFACTOR_MODE=inert` — metabolic cofactors (ATP, ADP, NAD+, H2O, Pi …)
  are pinned at baseline, so they remain AND inputs but cannot carry a
  perturbation. The networks still contain them: the generator represents what
  curators recorded, and this is the propagator declining to route through it.
  Set `propagate` for the previous behaviour. The list and the measured A/B are
  in `src/core/cofactors.jl`.
- `DS_SCC_METHOD=fixed_point` (default) | `minimize` | `pool` | `pool_parity` | `pool_all`:
  the last three treat a cyclic component as a conserved pool (specs/017);
  measured, refuted as a blanket rule by a traced case (product of entries
  destroys isoform redundancy), not adopted.
- `DS_SCC_BREAK_ROLES` (default empty): derived-edge recycling closures
  (assembly, depletion; catalyst) read at their entry value and excluded from
  component detection — specs/018. Measured: the generator-side fix absorbed
  it (the `assembly` arm is now ≈0), so it stays off as the record.
- `DS_COMPOSITION_MODE=assembly` and `DS_DEPLETION_OWN_PRODUCT=full` are the
  byte-identical defaults for two measured-but-not-adopted alternatives
  (`limit`: a `composition` edge can lower a container but never raise it;
  `limit_novel`: the same, skipping the 54% of composition edges that repeat a
  producing reaction; `suppress_only`: a depleter built from the target cannot
  de-repress it).
  Both A/B'd on the shared catalog, neither adopted; `specs/016` has the arms
  and the traced mechanisms.
- `DS_SELF_INHIBITOR_WEIGHT=0.1` (default; `off` restores the old product): an
  inhibitor that *contains* its own reaction's input (a sequestering complex
  such as WIF1:WNT, or RUNX1 mRNA:miR-675 RISC) counts that input twice under
  `divide`: one such inhibitor cancels it, two invert it. The part of the
  inhibitor's change that the shared input explains (their log-fold overlap)
  is kept only at weight `w`; the rest keeps full strength. A fully tracking
  inhibitor reads `x^(1-w)`; a partly tracking one is damped, not deleted. It
  applies only where the shared input reaches the inhibitor in the network,
  and it can only weaken an inhibitor. It needs the bundle's
  `containment.csv`; the solve reports `self_inhibitor_rule` (on / off /
  inert) and how many inhibitor slots were flagged. Measured numbers are in
  specs/023. specs/022 has the rule; specs/023 has the adoption.
- Export aggregation default: `stoichiometry_weighted`.

### Where design decisions live

This repo uses spec-kit. **`specs/NNN-name/` is the record of why**; CLAUDE.md
carries only current state. Before changing solver behaviour, read the spec
for the last feature that touched it rather than re-deriving from the code.

- `specs/002-upregulation-propagation/` — the AND/assembly defaults above.
  `research.md` holds the measured A/B including the negative results.
- `specs/010-and-multiplication-fidelity/` — the AND default above, the
  measured deviation from the specification, and the depth-invariance argument.
- `specs/011-depletion-suppression-bound/` — the depletion floor above, found
  by root-causing every gap case and tracing one; includes the floor sweep and
  the PIP3 regression it costs.
- `specs/012-self-contained-inhibitors/` — TWO NEGATIVE results on structural
  double-counts. An inhibitor containing its own substrate cancels the signal
  exactly (150 reactions); dropping those edges is −61 held-out, p<0.0001 —
  **but that arm's flag set was hash-seed dependent** (fixed 2026-09-25), so
  the −61 is not reproducible as recorded and must be re-run before citing. An
  entity that is both catalyst and substrate squares its own fold-change
  (7,013 reactions, 15.8%); deduplicating is −15 held-out, p=0.0015. Both
  mechanisms are real and both are load-bearing — on these networks,
  *bounding* a runaway operator has paid off and *deleting* a wrong term has not.
- `specs/003-solver-objective/` — **negative on accuracy.** Levenberg-Marquardt
  made the minimiser fast but it loses at every gamma (held-out −71 / −87 /
  −148). The mechanism findings stand (with gamma = 0 the all-zero state is a
  global minimum); the accuracy claim does not. `mu` and `gamma` are read
  **only** by `DS_SCC_METHOD=minimize`; under the default fixed-point method
  they are inert and are reported as `nothing` with `mu_gamma_read = false`,
  deliberately, so they cannot be read as having shaped the answer. One design
  term is still untested (FR4, soft observations), and `minimize` stays
  reachable because the parameter-learning direction needs it.
- `specs/006-bounded-derepression/` — the de-repression ceiling in the list
  above, and the epsilon attribution behind it.
- `specs/007-cofactor-conduction/` — the record for `DS_COFACTOR_MODE=inert`.
- `specs/013-solver-label-invariance/` — relabelling a verified-isomorphic
  network **moves predictions**, because sweep order inside a strongly
  connected component comes from Julia `Dict` hash order. Quantified, not
  fixed. It is the blocking prerequisite for parameter learning, which
  specs/014 and specs/003 state rather than 013 itself.
- `specs/014-loop-elasticity/` — the positive-cycle knife-edge and the
  sigmoid-epsilon arms.
- `specs/015-dissociation-sinks/` — released-subunit readout handles traced
  end to end; the largest class of severed curator routes, and not a lever.
- `specs/016-curator-oracle/` — read the pathway from Neo4j and diff it
  against the generated network. Also the composition-edge arms: **an
  Interferon α/β fix and a loss everywhere else.** The apparent +40 held-out
  is **−60 outside IFN α/β** (016 E4), and re-measured once the welds were
  gone it is **−185** (018), because composition edges are the derived class
  that closes cycles. Refuted, not merely concentrated — an earlier revision
  of this file said "declined on concentration rather than on harm", which was
  wrong.
- `specs/017-loop-pool/` — treating a cyclic component as a conserved pool.
  A traced case refuted it as a blanket rule; the `pool_all` arm refuted it on
  both axes. Nothing adopted.
- `specs/018-derived-edge-loops/` — **the largest accuracy finding.** Our own
  boundary expansion welded cycles Reactome does not contain (1,994 of 2,077
  cycle-carrying assembly edges). Fixing it is held-out **+173**, p < 1e-4,
  94% concentrated in DSB Repair. Records one residual defect, measured at 4
  edges in 2 pathways and therefore not a lever.
- `specs/019-sink-bridges/` — three interventions in the sink-bridge family,
  all null or negative. The Interferon α/β result (+100 fixed, 0 broken) is
  recorded as live and unclaimed.
- `specs/020-variant-node-sharing/` — variant sharing adopted: nodes −34.6%,
  edges −21.5%, no pathway gains cyclic nodes, held-out net zero on both axes.
  Also the flag-expiry policy: a flag is removed once its question is answered.
- `specs/021-empirical-holdout-axis/` — the phospho-site validation design,
  the target list, and what a magnitude claim can and cannot be. Open.
- `specs/022-self-contained-inhibition/` — inhibitors built from their own
  reaction's input; the `DS_SELF_INHIBITOR_WEIGHT` rule and its A/B.
- `specs/030-boundary-hierarchy/` — the generator decomposes a root complex
  one `hasComponent` level at a time, joining components the pathway already
  has (LNG #97, default on). Held-out **+197** (243 / 46, p = 1.6e-33), +175 of
  it IFN α/β. Records two traced losses: TP53 (uuid order in the loop) and 25
  STAT1 cases (root-copy preference).
- `specs/031-severed-route-anatomy/` — what the 252 remaining severed routes
  are, classified by script. 108 need composition hops into complexes no
  reaction forms (a curation-review list). 96 entity → reaction breaks are copy
  multiplicity; the targeted join is untried, near a family that failed three
  times. About 41 cases for an any-copy pool. Analysis only.
- `specs/009-solver-defaults/` — the one-variable-at-a-time re-measurement
  behind the `DS_*` defaults above (cited in that section too).
- `.specify/memory/constitution.md` — project principles the specs are
  checked against.

**Not listed here: 001, 004, 005 and 008.** They are OR semantics, loop
handling, node identity and cycle handling — superseded as decisions by 013,
014, 017 and 018, which are listed. The index covers the specs behind current
behaviour, not the whole directory; `ls specs/` is the complete list.

Per-feature numbers belong in that feature's `research.md`, not here.

### Testing Strategy
**Most of the test suite cannot fail.** Thirteen files contain assertions; the
other seven execute code and print output.

| file | assertions | note |
|---|---|---|
| `test/test_config_validation.jl` | 196 | |
| `test/test_loop_elasticity.jl` | 173 | + 1 `@test_broken` |
| `test/test_propagator_invariants.jl` | 121 | |
| `test/test_loop_pool.jl` | 81 | |
| `test/test_solver_determinism.jl` | 80 | |
| `test/test_and_curves.jl` | 71 | |
| `test/test_scc_break_roles.jl` | 56 | |
| `test/test_self_inhibition.jl` | 80 | specs/022, added 2026-09-25 |
| `test/test_cli_observations.jl` | 39 | |
| `test/test_api_errors.jl` | 44 | |
| `test/test_cycle_handling.jl` | 34 | + 2 `@test_broken` |
| `test/test_observation_pinning.jl` | 23 | |
| `test/test_worked_example.jl` | 9 | |

Counts are the `Pass` column of each file's outer `Test Summary`, not `@test`
occurrences — several testsets generate assertions in loops. **Enforced on
every CI run:** `.github/workflows/test.yml` runs every suite in
`test/asserting_suites.txt`, and `scripts/check_doc_counts.py` fails the build
if this table disagrees with what they print (it caught the 2026-09-25
changes). Locally, the loop in the test-runner service re-reads them.

Earlier revisions of this table were wrong in several places at once: they said
"three files" while listing twelve, split the table with a stray blank line,
and gave `test_and_curves.jl` as 17 (actually 71), `test_cli_observations.jl`
as 20 (39) and `test_cycle_handling.jl` as 39 + 2 (34 + 2). Re-read them from
CI rather than trusting a figure in this file.

**A failing testset used to hide every later one.** A top-level `@testset`
throws when it finishes with a failure, which aborts the file. In the dev
container, whose `DS_*` overrides then contradicted the code defaults (since
fixed), that meant
`test_config_validation.jl` ran 10 assertions and silently skipped **141** —
the cofactor tests, the silo-bridge tests and the observation-membership test
all looked green because they never executed. The files now nest their testsets
inside one outer testset so failures accumulate and everything runs.

`test_basic.jl`, `test_steady_state.jl`, `test_hill_function.jl`,
`test_feedback_loops.jl`, `test_inhibition_focused.jl`,
`test_pathway_propagation.jl` and `test_biological_fixes.jl` execute code and
print output; they pass whatever the code does. Treat a green run of those as
"it did not throw", nothing more.

**The docker runner used to be a false green and is now fixed.** It ran
`test_basic.jl` and `test_steady_state.jl` — zero assertions between them —
then chained a `test_full_pipeline.jl` that does not exist, and printed "All
tests passed". It now runs the suites in `test/asserting_suites.txt`, the same
list CI reads, and mounts `Project.toml` so the package can actually load.

**The dev container could not load DeltaSignal at all until 2026-09-22.** The
compose file mounted `src`, `cli`, `test`, `bench` and `examples`, but **not
`Project.toml`** — so the dependency set is frozen at image-build time while
the source is live. `/app/Project.toml` in the running image is dated May 25
and predates the `SparseArrays` dependency, so every test file dies at load
with `Package DeltaSignal does not have SparseArrays in its dependencies`.
CI was unaffected: it instantiates from the repo. **Both `julia-api` and
`test-runner` now mount `Project.toml`**, which is sufficient for this failure:
it was a `[deps]`-table check and `SparseArrays` is a stdlib, so no
`Manifest.toml` entry is needed. Verified 2026-09-25 — `julia-api` had been
crash-looping on exactly this error, and recreating it with the mount brought
`/api/health` back in about twelve seconds. A genuinely new **external**
dependency still needs `docker compose -f docker-compose.dev.yml build`.
`Manifest.toml` is gitignored, so mounting it is not an option.

When adding behaviour, add assertions to one of the files that assert,
or start a new one — do not extend a file from the zero-assertion list and assume it is
covering anything.

## Project Dependencies
- Julia 1.10+ (LTS required)
- Key packages: Optim.jl, JSON3.jl, DataFrames.jl, CSV.jl, HTTP.jl, ArgParse.jl
- Docker for containerized deployment
- UI: separate — Angular project in the WebsiteAngular workspace, consuming this repo's HTTP API (`docs/API.md`)