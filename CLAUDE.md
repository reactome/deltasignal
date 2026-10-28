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

**`DS_PIN_ROOT_IGNORE=depletion` (default since 2026-10-15, specs/061
adoption)**: an in-edge of type depletion does not stop a node being a root,
so a form fed only by depletion is pinned as a root (as in MP-BioPath's
networks, specs/058). `DS_PIN_ROOT_IGNORE=none` restores the specs/024 rule,
for reproducing numbers from before 2026-10-15.

**Rule decisions set RAF aside** (specs/054 D1, Adam 2026-10-13): its
readouts are a bistable switch, not signal. Gates on any case axis are read
with `bench/analysis/holdout_report.py --compare B --set-aside
RAF_MAP_kinase_cascade`, and both figures are reported. Headlines keep RAF.

**`DS_KO_AGG=mean_entity` (default since 2026-10-26, specs/105, adopted at the freeze, specs/111)**: a readout
that resolves to several nodes weighs each curated entity once (a proxy reaction copy by its output entity, a set
by its member), not each variant copy. `DS_KO_AGG=mean` reproduces numbers from before 2026-10-26.

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
  `DS_HILL_LOG_ZMAX=10.0`. Behaviour is pinned by `test/test_and_curves.jl`,
  including the band just below the 100× cap (specs/074).
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
- `DS_SCC_ORDER=flow` (**default since 2026-10-03**, specs/013 amendments 2–3)
  | `structure` | `label`: the Gauss-Seidel sweep order inside a cyclic component.
  `flow` sweeps downstream, by breadth-first distance from the component's
  entries, with ties broken by a stable-id structural colour. It is
  label-free, so rebuilding a catalog (which re-mints every uuid) no longer
  moves predictions. Values agree to `pred_ui` precision, not bit for bit:
  tied positions keep label order, and RAF/IFN differ by up to 6e-16. The order
  depends on Julia's `hash`, so re-check invariance when the Julia version
  changes. Under `label`, a rebuild of identical content moved 158
  TP53 predictions. Held-out and convergence are unchanged. The solve reports
  `scc_order` and `scc_order_ties`. `DS_SCC_SWEEP=jacobi` is also
  label-free, but costs held-out −94 and converges less (amendment 1).
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
- `DS_DEPLETION_OWN_PRODUCT=suppress_only_entity` (specs/060 amendment 1; not
  a default): as `suppress_only`, and a depleter whose entity CONTAINS its
  target's entity (containment table) also counts as own product, so a
  depleter built from another copy of the target (SMAD2:SMURF2 -| recycled
  SMAD2) cannot de-repress it. The solve reports `own_entity_depletions`.
- `DS_DEPLETION_OWN_PRODUCT=suppress_only_source` (**default since 2026-10-15**,
  specs/061, adopted with `DS_RECYCLED_TWIN=on` and the bench's
  `DS_PIN_ROOT_IGNORE=depletion`; held-out +32 (33 / 1), experimental +7 (7 / 0)
  against the previous defaults; `full` restores them):
  as `suppress_only`, and an entity-own depleter is capped at 1 only while a
  SOURCE copy (a copy of the target's entity it is built from, found walking
  back along material edges, not catalyst or regulator ones (amendment 2);
  never the target itself, amendment 1)
  reads below baseline. **Amendments 3–5 (2026-10-15): no cap is applied
  during an iteration.** Pass 1 solves with every own-product depletion free
  (as `full`); an entity-own edge whose depleter AND a source copy fell in
  pass 1 is removed; pass 2 solves without it, again uncapped. Node-copy own
  products are NOT capped under this mode (as `full`): read live, the
  `min(factor, 1)` cap and the source switch made TP53 and TGF-beta oscillate,
  and the scored value was the iteration cap's choice (specs/061 R5). The
  solve reports `own_source_depletions` (edges removed) and
  `own_product_removed`. A SMAD2 knockout
  cannot de-repress SMAD2 through SMAD2:SMURF2; an MDM2 knockout still
  de-represses TP53 through MDM2:TP53. The solve reports
  `own_source_depletions`.
- `DS_COMPOSITION_MODE=assembly` (default) and `DS_DEPLETION_OWN_PRODUCT=full`
  (the default until 2026-10-15) are the baselines for two measured-but-not-adopted alternatives
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
- `DS_SELF_ACTIVATOR_WEIGHT=0.1` (**default since 2026-10-20**, specs/083,
  adopted on faithfulness: no class change on either axis, 322 curator values
  move >5%) | `off` (the old product) | w in [0, 1]: the mirror of specs/022 for a positive regulator or
  catalyst that contains its reaction's own `input` and is computed from it (a
  TF:gene complex regulating that gene's transcription). Its overlap with the
  input is kept at w/n, so the input is counted once rather than once per such
  activator. Reports `self_activator_rule` and `self_activators`.
- `DS_DRUG_MODE=inert` (**default since 2026-10-13**, specs/055; curator
  held-out +8, 8 fixed / 0 broken, experimental +2 with RAF set aside) |
  `propagate` (specs/032): drug-derived nodes (a Reactome `Drug`, a complex
  containing one, a set that is all drugs), listed by the generator's
  `drugs.csv`, held at baseline like cofactors, for a cell without the drug.
  An observation still wins. The solve reports `drug_rule` and `drugs_held`;
  a bundle without the table reports `inert: no drug table` and the bench
  refuses it, so re-run such a build with `--server DS_DRUG_MODE=propagate`.
  RAF moves from its switch's high state to its low one (specs/054).
- `DS_PATHOGEN_MODE=inert` (**default since 2026-10-10**, specs/048 amendment
  2) | `propagate`: nodes listed in the generator's `pathogens.csv` are held
  at baseline, as drugs are, for a cell without the pathogen. An observation
  still wins.
  - Under the generator's `LNG_PATHOGEN_PROTEIN` (default on), the list holds
    only entities that carry a pathogen PROTEIN (an EWAS on a
    `ReferenceGeneProduct`): a protein leaf, a complex with such a component,
    or a set whose members all are, e.g. SARS-CoV-2 N:M:PDPK1.
  - Viral RNAs are not listed. They are the ligand DDX58/IFIH1 senses, and
    holding them broke 12 curator cases, net −11 (arm `path`).
  - Builds from before 2026-10-10 have no `pathogens.csv`, so the bench refuses
    them under the default; re-run one with `--server DS_PATHOGEN_MODE=propagate`
    (the specs/047 MP-BioPath-network run included).
  - Against `propagate`: experimental +7 (all PDPK1 knockout), held-out +1,
    nothing broken.
  - The solve reports `pathogen_rule` and `pathogens_held`. A bundle without
    the table reports `inert: no pathogen table`, and the bench refuses it.
- `DS_UBIQUITIN_MODE=propagate` (default) | `inert` (specs/098, **not adopted**: held-out +6 USP1, experimental −4 CCNE1; free ubiquitin is the only
  link of the p27 degradation cycle): free
  ubiquitin nodes (the three ubiquitin stIds exactly, not conjugates) held at baseline like
  drugs, so the free pool cannot close loops over a pathway's ubiquitinations. Reports
  `ubiquitin_rule` and `ubiquitin_held`.
- `DS_SET_POOL_MODE=product` (default) | `extreme` | `geomean` | `mean` | `max`
  (specs/033): how a set POOL node (the generator's `LNG_SET_POOL`: one node
  per set-valued catalyst or regulator, fed by `set_member` edges) combines its
  members. It has no effect on networks without `set_member` edges.
- `DS_CONSERVED_MODE=off` (default) | `inert` (specs/035): a loop input released
  UNCHANGED by its own loop (every producer consumes a complex containing it:
  RAS GAPs, scaffold partners, F-actin) is held at baseline like a cofactor.
  Signal carriers transformed inside the loop are not touched. The solve reports
  `conserved_rule` and `conserved_held`.
- `DS_CYCLE_MODE=balance` (**default since 2026-09-28**, specs/039 amendment 7;
  held-out +77, experimental 0 against `off`) | `off`, with `DS_CYCLE_PHI=0.1`
  (explained in docs/MODEL.md §4): each interconversion pool in the
  generator's `pools.csv` is solved as a conserved Markov chain (π = πP)
  instead of being iterated. A modification curated in steps (bind, modify,
  release) is one pool: its enzyme complexes are intermediates with no
  baseline share that read the flux through them (amendment 2).
  `DS_CYCLE_CARRIERS=on` (default; `off` for the isolation arm) makes the
  enzyme's free form read its producers other than the release steps of
  EVERY pool that lists it (specs/049 amendment 7: the RAS GAPs, carriers of
  four pools, were written from the other pools' release flux and RAF never
  converged; fixing it is RAF tuning +21, experimental +14, nothing else).
  `DS_CYCLE_REEVAL_CARRIERS=carrier` (default since 2026-10-18, specs/075;
  `full` restores the old drives) applies the same rule where another pool's
  drives re-evaluate the carrier. The solve reports `cycle_rule`, `cycle_pools_solved`, `cycle_pools_multistep`,
  `cycle_carriers` and `cycle_carriers_held`; under any `DS_SCC_METHOD` other than
  `fixed_point`, or `DS_SCC_SOLVE=0`, pools are not solved and `cycle_rule`
  says `balance: inert under ...`. It needs the generator's pool
  tables; a bundle without them reports `balance: no pool table` (loops then
  iterate as under `off`), and the bench refuses such a build. φ₀ refers to
  the pool's ACTIVE form, which the generator marks (`LNG_POOL_ACTIVE_VIA`).
  A pool may mark several base (resting) states (specs/081); a `recycle` edge
  (generator `LNG_RECYCLE_POOLS`, default on since specs/081) joins a modification cycle into one
  component and is never an input. The solve reports `recycle_edges`
  (`src/core/recycle.jl`); the bench never counts it against a root.
  `DS_POOL_BASE_WEIGHT=split` (**default since 2026-10-23**, specs/091; tuning +35,
  held-out 0, experimental +1 against `full`) | `full` (specs/081, the earlier
  default): a pool's k resting copies weigh 1/k each, so the next form rests at
  φ₀ whatever k is (β-catenin's ceiling was 181×). Reports `pool_base_weight` and
  `cycle_bases_split`.
  `DS_POOL_SUPPLY_AND=off` (default) | `on` (specs/090, **not adopted**: tuning −6, RBL2
  overexpression reads UP through DREAM ⊣ RBL1): a pool
  state's AND-joined producers (a recycle base assembled from the protein and its
  partners) are one supply term, their capped product. Reports `pool_supply_and`
  and `cycle_supply_and_groups`.
  `DS_POOL_DRIVE_TWINS=off` (default) | `hold` (specs/104, measured, no class change on either axis;
  not adopted at the freeze): a pool's drives hold at baseline any non-pool copy of the pool's own entities that
  a step reads (WNT: root copies of the PP2A step's substrate read 80x under APC overexpression and
  emptied the intermediate the ubiquitination readout leaves from, 0.014 -> 1.11). Reports
  `pool_drive_twins_rule` and `pool_drive_twins`.
  A release pool (generator `LNG_RELEASE_POOLS`, default on since specs/086) lists
  the released copy as its carrier, reading the free state through a
  `released_twin` edge; paths may be up to 10 steps (`CYCLE_MAX_STEPS`).
- `DS_SELF_FED_MODE=off` (default; `entry` measured and **refuted**: TP53 −190,
  PIP3 −126, experimental −66) | `entry` (specs/040 rule A): in an
  iterated cyclic component, an input fed by nothing except what it feeds
  (every signal path from a component entry passes through the input's own
  consumers; set pools transparent) reads its component-entry value, which
  is baseline. Entries are per solve: pins off baseline and targets with an
  outside activator off baseline. The solve reports `self_fed_rule`,
  `self_fed_nodes` and `self_fed_edges_held`. `catalyst_loop` (specs/109, measured, gate passes, not a
  default) runs the same test only in positive components closed solely through catalyst or
  regulator slots (autocatalysis, a TF switching on its own gene) and holds only those edges;
  it also reports `self_fed_catalyst_components`. Also fixed there: the final
  residual now reads held edges (specs/018 closures too) at the value the
  iteration used, not the live state, which had reported residual ≈ 1 on
  converged solves.
- `DS_SELF_INHIBITOR_COPIES=on` (**default since 2026-10-23**, specs/095; tuning +14,
  experimental +9, held-out 0, nothing broken against `off`) | `off`: the
  specs/022 pair test also accepts an inhibitor that another copy of the input's
  entity reaches (PTEN's mRNA:miRNA complexes, built from the transcribed mRNA copy,
  repressing translation of the root copy). Reports `self_inhibitor_copies` and
  `self_inhibitor_copy_pairs`.
- `DS_SELF_INHIBITOR_LEAVES=1` (**default since 2026-09-28**; curator held-out
  +48, experimental +20 at p 0.0012 against `0`) | `0` (specs/040 rule B): the
  specs/022 pair test also matches an inhibitor that shares a non-cofactor
  leaf with the input when a carrier of that leaf reaches both (built from
  the same species). Same weight, weaken-only. Reports
  `self_inhibitor_leaves` and `self_inhibitor_leaf_pairs`.
- `DS_CARRIER_LOOP_MODE=closed_twin` (**default since 2026-10-12**, specs/049
  amendment 6; held-out +122, experimental −2 at p 0.73 against `off`) | `off` |
  `closed` | `closed_supply` | `all` (`closed`: held-out +81, experimental −1;
  `all`: held-out −55). A carrier loop is a protein bound
  and released unchanged (E → E:X → E), which the generator's `find_pools`
  rejects as a pool and lists in `carrier_loops.csv`. Under `closed`, the free
  form of a loop that nothing outside it produces loses its activating
  in-edges from the loop and so reads baseline: total carrier is conserved,
  and bound forms hold a negligible share. A silo bridge (off by default) or a
  boundary edge into it would still feed it. `closed_supply` (amendment 4)
  also drops depletion in-edges into closed free forms; `closed_twin`
  (amendment 5) adds to that an in-edge from each root copy of the same
  entity, so a closed carrier follows its own gene. Under `all`, every listed free form reads only its
  producers outside the loop. Bound forms are not touched. The solve reports
  `carrier_loop_rule`, `carrier_loop_forms`, `carrier_loop_edges_dropped` and
  `carrier_loop_twin_edges`.
  A bundle without the table reports `<mode>: no carrier table` and solves as
  `off`, and the bench refuses it.
  `DS_CARRIER_COMPONENT_TWIN=on` (**default since 2026-10-21**, specs/084
  amendment 1, post hoc; tuning +6, experimental +4, held-out 0, nothing broken)
  | `off`: under `closed_twin`, a closed free form with no starting copy of its
  own entity and no activating in-edge left after the drop (a complex the
  pathway only regenerates, e.g. the ATM-carrying DSB complex in Cell Cycle
  Checkpoints) gains AND in-edges (`carrier_component`) from the starting
  copies of the maximal entities it contains. Reports
  `carrier_component_rule` and `carrier_component_edges`.
  `DS_CARRIER_LOOP_MOVES=drop` (default) | `keep` (specs/093, **not adopted**: held-out −10; CBL in MET fixed, HGF/MET knockout rail UP): a release
  the generator marks `moves` in `carrier_loops.csv` (it carries the free form across
  compartments: MET recycled to the membrane) is left out of the loop table, so no
  carrier-loop rule drops it. Reports `carrier_loop_moves`.
  `DS_CARRIER_LOOP_OPEN=fed_single_ub` (**default since 2026-10-23**, specs/088 amendment 1;
  tuning +2, experimental +1, held-out 0 against `fed_single`: ubiquitin is not a displacer,
  so a ligase released after ubiquitinating its target, MDM2:MDM4, stops reading its own
  release) | `fed_single` (default 2026-10-19 to 2026-10-23, specs/077
  amendment 1, post hoc; held-out +65 (84 / 19), experimental −5 (every loss
  traced: release artefacts or threshold shifts), faithful build +272) | `off`
  (the earlier solver) | `fed` (fails its gate: intrinsic apoptosis −26).
  `fed_single` is `fed` that keeps a release with a second non-cofactor input
  (a displacer, BAD freeing tBID). The rule: an OPEN loop's free form that keeps a producer outside
  its loops stops reading its own release (RBBP8 released by MMEJ in DSB
  Repair on the faithful build); a form fed only by loop steps (the IFN α/β
  receptor scaffold, which `all` cut) keeps them. Reports
  `carrier_loop_open_rule` and `carrier_loop_open_edges`. The table comes from LNG #111, so builds
  before it (`20261010-1415_c0e1b2a` and earlier) are re-run with
  `--server DS_CARRIER_LOOP_MODE=off`. Canonical build `20261026-1349_8b6527e` (the frozen build, specs/111: generator tag `freeze-2026-10-26`, identical predictions); before it `20261024-2211_7f4b28a` (since 2026-10-24, specs/102: Release 97 plus the summation supplement, LNG #127; Release 97 alone is `20261024-1634_f2c6a71` (specs/101, LNG #126); `20261024-0801_7a754fc` (specs/100, adds `rnas.csv`), `20261022-2011_233665c` (specs/089), `20261021-2218_f99c4ab` (specs/086) and `20261020-2038_23ef3e1` (specs/081) before).
- `DS_RECYCLED_TWIN=on` (**default since 2026-10-15**, specs/061 adoption) |
  `off` (specs/059): a RECYCLED copy of an
  entity (every producer, except depletion and dissociation, downstream of
  itself) that no root copy of the same entity reaches along non-depletion
  edges gains an OR in-edge (`recycled_supply`) from each root copy (in-edges,
  if any, all depletion; never a carrier free form, specs/061 amendment 2). The one-way form of the identity twin that
  `closed_twin` handles for carrier loops (SMAD2 in TGF-beta). The solve
  reports `recycled_twin_rule` and `recycled_twin_edges`.
- `DS_DEPLETION_COMBINE=product` (default) | `share` (specs/103, **not adopted**: held-out 0, tuning
  +13, experimental 0, PIP3 −13: VAV1 through PTEN's four ligases): several depletion edges into one
  node read as consumers of one pool, b / mean(x), instead of the product of their factors. MDM4
  knockout then de-represses TP53 (9 curator cases). Reports `depletion_combine` and
  `depletion_shared_nodes`.
- `DS_METABOLITE_POOLS=off` (default) | `on` (specs/107 rule 1, measured: no class change on either axis; fails B4 by two
  pairs, PIP2 under PI3K feedback; adoption Adam's): small
  molecules converted into one another by one-in/one-out reactions (LACT ⇄ PYR, R5P ⇄ RU5P) are one
  pool solved by the cycle balance, every form resting with an equal share. Reports
  `metabolite_pools_rule` and `metabolite_pools`.
- `DS_CARRIER_MULTISTEP=off` (default) | `strict` | `wide` | `compartment` (specs/107 rules 2–3, measured, gate fails: `strict` held-out +6 is
  free ubiquitin through USP1, the specs/098 gain, and B4 −2; `wide` held-out −50, SMAD2/3 phosphorylation read as a
  carrier and free ubiquitin in Mitotic G1): a
  form carried round a path of up to 10 reactions and released unchanged (BAP31, microtubules) joins the
  carrier-loop table, so the carrier rules apply. `strict`: protein carriers, no set nodes; `wide`: also
  small-molecule carriers (PE through LC3:PE), set nodes, and copies of one reference entity across
  compartments or modifications, from the generator's `leaf_refs.csv` (LNG #128; without it that part
  reports `wide: no reference table`). `compartment` (amendment 1, measured, guard fails: held-out −10, all Class I MHC knockouts DOWN → NORMAL once
  the BAP31 loop stops railing; rails 32% → 20%, behaviour checks pass): `wide`, joining copies
  only when they also share modified residues (the table's `modifications` column, LNG #129; without it
  `compartment: no modification column` and the bench refuses), and never free ubiquitin as a carrier.
  Reports `carrier_multistep_rule`, `carrier_multistep_forms`, `carrier_multistep_rows`.
- `DS_RNA_TWIN=on` (**default since 2026-10-24**, specs/100; tuning +28 (28 / 0), experimental
  +4 (4 / 0), held-out 0: every TP53 case in PIP3) | `off`: a root copy of an RNA
  species (the generator's `rnas.csv`, LNG #125) reads its produced copies (`rna_supply`, OR, one-way,
  solver-only), so TP53 raising PTEN transcription reaches translation. Reports `rna_twin_rule`
  and `rna_twin_edges`; a bundle without the table reports `on: no RNA table` and the bench
  refuses it, so re-run an older build with `--server DS_RNA_TWIN=off`.
- `DS_SEQUESTRATION=off` (default) | `deadend_sole` | `deadend` | `sole_direct` (specs/063,
  all **refuted**: held-out −1,180 / −1,686 / −1,121): a binding whose complex nothing in the pathway uses
  (only dissociation sinks) gains a depletion edge from the complex to its
  partner, so binding takes the partner out of circulation (follistatin,
  USP18, IL13RA2). `deadend_sole` only where one input is used nowhere else
  (the sequesterer) and exactly one is a partner; `sole_direct` draws the edge
  from the sequesterer instead of the complex. Reads `nodes.csv`; reports
  `sequestration_rule` and `sequestration_edges`, and the bench refuses
  `<mode>: no node table`.
- `DS_ZERO_SUPPLY_MODE=off` (default; `hold` measured and **not adopted**,
  specs/052) | `hold`: under a pin to 0, a
  node built from a leaf (protein or small molecule) that no unpinned source
  can supply is held at 0
  (a knocked-out receptor's complexes cannot be regenerated by a cycle).
  Reads the generator's `nodes.csv` (`node_kind`, `member_leaves`); only
  entity kinds are held, never pools or reactions; only material-moving
  edges carry a protein. The solve reports `zero_supply_rule`,
  `zero_supply_leaves` and `zero_supply_held`; a bundle without `nodes.csv`
  reports `hold: no node table` and the bench refuses it.
- `DS_NONCONV_FALLBACK=lm` (**default since 2026-10-16**, specs/066; held-out
  +1 (2 / 1), experimental 0, converged solves 1,654 → 1,682 curator and
  222 → 234 experimental) | `off`: a cyclic component still above tolerance
  at the iteration cap is re-solved by Levenberg-Marquardt on |F(x) − x|²
  (no baseline prior; `DS_NONCONV_LM_ITERS`, default 400), started from the
  last iterate, and the result is kept only if its residual is lower.
  Converged components are never touched; components holding cycle-balance
  pool nodes are handled by `DS_NONCONV_POOLS` below. Why: Goodwin rings (a complex depleting each of its
  own free subunits across an AND, from the generator's ubiquitin pass, the
  specs/057 class) have a fixed point that repels any damped relaxation, so
  the iteration cycled and the cap chose the class. The solve reports
  `nonconv_fallback`, `nonconv_components`, `nonconv_improved`,
  `nonconv_rescued` and `nonconv_skipped_pools`. DSB Repair's 454-node
  component improves but still ends above tolerance (RNF8/RNF168). That is
  no longer the cap's choice: across caps 500–4000 none of the 638 cases
  still flagged non-converged changes class (specs/066 amendment 2).
  `DS_NONCONV_POOLS=joint` (**default since 2026-10-18**, specs/073; canonical
  identical, faithful build TP53 non-converged 663 → 102) | `skip`: a
  component holding pool nodes is re-solved too, with each pool-written
  node's row G(x) − x (G the pool update). Reports `nonconv_pools`,
  `nonconv_pool_components`, `nonconv_pool_improved`, `nonconv_pool_rescued`.
  `DS_NONCONV_RESTARTS=4` (**default since 2026-10-19**, specs/076; canonical
  no class change, converged curator 23,117 → 23,284; faithful build every
  TP53 case converges; `0` restores the previous solver) with
  `DS_NONCONV_RESTART_SWEEPS=250`: a component the fallback leaves above
  tolerance is iterated that many more sweeps from its last iterate and the
  fallback retried, up to N times; the best state is kept.
  Reports `nonconv_restarts_max`, `nonconv_restarts`, `nonconv_restart_rescued`.
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
  the target list, and what a magnitude claim can and cannot be. The decryptM
  drug axis (`bench/decryptm_bench.py`): the frozen Release97 rules **lose to
  reachability** on its held-out regulated cases (45 vs 70 of 190; by-drug 95%
  interval [−0.194, −0.053]); 116 of the 190 cases are unreachable inside the
  pathway. Under the measured kinobeads target map (specs/051 R4): 86 vs 115
  of 278. decryptM 1 has been FIT since that first scoring.
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
- `specs/032-inert-drugs/` — `DS_DRUG_MODE`: drug-derived nodes held at
  baseline. Not adopted then (held-out +30, experimental −5, all RAF);
  adopted in specs/055 once RAF was set aside.
- `specs/033-set-pool-nodes/`, `specs/038-pool-max/` — one OR node per
  set-valued catalyst or regulator (`LNG_SET_POOL`), feeding the reaction once.
  Adoptable on faithfulness under `DS_SET_POOL_MODE=product`: held-out +6 / −2,
  experimental +3 / +8 (the +8 concentrated in 4 perturbations) on two
  measurements. `max` (any member suffices) is refuted (held-out −120,
  experimental −57); `mean` does not gain.
- `specs/034-regulator-wiring/` — inventory of how catalysts and regulators are
  wired; the RAF collapse traced by iteration; MP-BioPath's network and its
  1,800 hand edits compared with ours. UUIDs already resolve 86% of their
  loop artefacts. Our loops are mostly curated interconversion cycles, and the
  219 inhibitions they added are mostly the reversal enzyme of such cycles
  (187 absent from ours). The design input for what comes next.
- `specs/035-conserved-participants/` — `DS_CONSERVED_MODE`: loop inputs
  released unchanged held at baseline. Not adopted: it held PIP3 itself, and
  experimental was −122.
- `specs/036-capped-bundles-as-pools/`, `specs/037-combined-raf-fixes/` — the
  variant cap bundles a set's alternatives into one all-required node
  (`LNG_CAP_POOLS` fixes the structure; default on in the generator since specs/045). Combined with de-duplication
  it was held-out −162; de-duplication is the damage. Not adopted.
- `specs/040-self-feedback/` — a node's own downstream multiplied back into
  it (RAF traced: recycled set-pool members, bundle leaves fed only by their
  step's downstream, a leaf-sharing inhibitor). Two blind derivations; the
  rules `DS_SELF_FED_MODE` and `DS_SELF_INHIBITOR_LEAVES` above are
  pre-registered there. **Rule B adopted** (held-out +48, experimental +20,
  worst pathway −3); rule A refuted (TP53 −190, PIP3 −126).
- `specs/046-variant-complex-nodes/` — a complex holding sets becomes one node
  per member combination, and each reaction one copy per consistent choice
  across all participants, catalysts and regulators included (generator
  `LNG_VARIANT_NODES`, default on since LNG #109, canonical build
  `20261009-2307_132d4b8`; `0` restores the old networks). Arm vn7: held-out **+130** (p 1e-12), tuning −18, experimental
  −10 (ns), coverage identical, cut hand-offs 2,504 → 0. Two failed arms
  (vn4, vn5) are recorded with every mechanism traced. The solver reads the
  new `variant_pool` / `variant_split` edges as OR inputs (mean);
  `test_set_pool.jl` pins that.
- `specs/047-experimental-gap/` — the experimental gap to MP-BioPath, case by
  case. Our propagator on MP-BioPath's own networks scores 610/849 and
  matches or beats it in 8 of 10 pathways: **the gap is network
  construction** (RAF, Mitotic G1, Cell Cycle Checkpoints, TP53), traced to
  five mechanism classes. `bench/analysis/experimental_gap.py` reproduces it.
- `specs/048-experimental-levers/` — two levers from 047.
  - **Pathogen-derived nodes held at baseline: adopted**, as the protein-only
    rule of amendment 2. The first rule also held viral RNA and broke 12
    cases in DDX58/IFIH1 (net −11).
  - **Homodimer binding fix (`LNG_BIND_STOICH`): adopted** on faithfulness;
    it moves no prediction.
  - Canonical build `20261010-1415_c0e1b2a` (LNG #110, DS #101): experimental
    591/849 (+7), held-out +1.
  - Pools that share intermediates (the RAF knockout→UP class) continue in
    specs/049.
- `specs/049-shared-intermediate-pools/` — pools that share intermediates
  (generator `LNG_POOL_SHARED`, solver `update_pools!`) and carrier loops
  (`DS_CARRIER_LOOP_MODE`). **`closed_twin` adopted as the default**
  (amendment 6); every earlier arm was not.
  - `closed` gains curator held-out +81, but experimental is −1, losing the
    cases whose signal ran through a carrier drained to zero.
  - `all` is refuted (IFN α/β −168).
  - The shared pools' RAF knockout fix was itself a 14-3-3 drain. With the
    drain removed, a drug-complex de-repression (specs/032) reads RAF UP.
  - The iteration cap decides no canonical prediction.
  - Amendment 3 (P8): `closed` + `DS_DRUG_MODE=inert` on K (`kcd`) passes the
    registered gate (held-out +89, experimental +1 at p 1), but its RAF moves
    are a loop collapse; **not adopted** (Adam, an override of the rule). On S (`scd`) the
    ARAF/BRAF knockouts read DOWN through the real path, and RAS
    overexpression is lost (experimental −6).
  - Amendments 4–5 (P11–P12): `closed_supply` (also drop depletion into
    closed free forms) fails the old gates. `closed_twin` (also a closed free
    form reads its gene's root copy) gives held-out +122, experimental −2
    (p 0.73), and passes the faithfulness gate registered in amendment 5;
    **adopted as the default** (amendment 6, Adam). Canonical build
    `20261012-1401_01206ca` (LNG #111): curator held-out in release 89.29%,
    experimental 589/849, reproducing `ktwin` on every prediction. Open: an enzyme knockout should raise its substrate (SMAD7,
    PTPN11).
- `specs/051-dataset-roles/` — every dataset is FIT, SELECT or TEST
  (`bench/datasets.tsv`; `bench/dataset_roles.py` `require_role` refuses
  fitting on a non-FIT set; every fitting script must call it, as
  `bench/analysis/train_116.py`, `train_117.py` and `lr_118.py` do). **The curator held-out split is SELECT**, not a test set: it has chosen
  rules since 2026-09. decryptM 1 is FIT after its one scoring.
  The degree-preserving rewired-network null (`bench/rewire_network.py`)
  loses to the real wiring on every axis, experimental included (R3).
- `specs/052-zero-supply/` — `DS_ZERO_SUPPLY_MODE`: a protein that is not
  there is in no complex. From the decryptM EGFR trace (the curated CBL
  cycle rails at 100× under an EGFR knockout). **Not adopted**: curator
  held-out −20 (p 1e-5), experimental +4 n.s. Traced: the HDR loss is root
  pinning zeroing partner proteins, the PTK6 loss a pin cutting a cycle
  pool's recycling path (the specs/039 supply rule treats every pin as a
  reach barrier; not exact-zero numerics), the RAF gain NF1 (enzyme
  knockout raising its substrate).
- `specs/053-nf1-gap-twin/` — why an NF1 knockout cannot reach RAS (the GAP
  is an identity twin fed only by its own release), and the generator join
  that merged such copies (`LNG_CLOSED_TWIN_JOIN`, LNG PR #112, unmerged).
  **Not adopted:** curator held-out −493 (DSB −488); a merge is two-way and
  re-welds cycles. NF1 still fails: one of ten GAPs, and the RAF inversion.
- `specs/054-raf-rail/` — RAF readouts are the high state of a bistable
  switch (threshold about 3.5% either way on HRAS), so every RAS/RAF
  perturbation rails p-MAPK at about 50×. The loop is curated (input-and-
  catalyst complexes; ERK negative feedback lives in another pathway).
  **Decision (Adam, 2026-10-13): RAF is set aside from rule decisions**;
  gates are read with `holdout_report.py --set-aside RAF_MAP_kinase_cascade`
  and both figures reported. Re-read: rule B's experimental +20 is +5 (p 0.06)
  without RAF; inert drugs' −5 is +2.
- `specs/055-inert-drugs-canonical/` — `DS_DRUG_MODE=inert` re-measured on
  the canonical build under specs/054 D1: held-out +8 (8 / 0, all KIT
  overexpression through KIT:TKI complexes), experimental +2 with RAF set
  aside (CDK4 knockout), nothing broken outside RAF. **Adopted as the
  default** (Adam, 2026-10-13), mainly because benchmark cells have no drug.
- `specs/056-reactome-context/` — each pathway solved with the pathways that
  feed it (`bench/context_build.py`, a separate catalog build; a second
  protocol, never the headline). **Not supported:** feed-in curator held-out
  −164, experimental −14 (RAF set aside); with feedback −577 / −30. The
  cell-cycle cases gain, but importing a neighbour's free forms (EGF at the
  10× de-repression ceiling, KDM1A drained to 0) breaks more. Compare such
  arms with `holdout_report.py --pair-on-valid`.
- `specs/057-ub-substrate/` — a ubiquitination depleted every partner of
  a complex modified in place (CBL monoubiquitinating MET depleted free GRB2,
  inverting every MET/HGF case). Generator `LNG_UB_SUBSTRATE_CHANGED`
  (default off, LNG PR #113): 301 of 727 depletion edges. **Not adopted:**
  curator held-out −39 (p 0.0011); with the "Ub-X" complex test
  (`LNG_UB_SUBSTRATE_UBX`, amendment 1) −15 (49/64, p 0.19), experimental 0.
  MET +30 as predicted; CBL-family regulators (−38) were right only through
  the removed edges.
- `specs/058-depletion-roots/` — bench `DS_PIN_ROOT_IGNORE=depletion`: a
  form fed only by depletion counts as a root (as in MP-BioPath's networks).
  Alone not adopted: held-out +9 (33/24, p 0.29), experimental +1. All losses
  are SMAD2/3 (TGF-beta), which score only through the specs/024 fallback.
  **Default since 2026-10-15** as part of the specs/061 adoption.
- `specs/059-recycled-twin/` — `DS_RECYCLED_TWIN` (default on since
  2026-10-15, with specs/061): a recycled copy of an entity reads its root
  copies (one-way, like `closed_twin`). Alone: held-out 0; with specs/058's pins held-out +21
  (p 0.0025) but TGF-beta −12 (SMAD2/3 knockouts read UP: SMAD2:SMURF2, built
  from SMAD2, de-represses SMAD2 when the knockout zeroes it).
- `specs/060-own-product-depletion/` — `DS_DEPLETION_OWN_PRODUCT=suppress_only`
  re-measured on the canonical build: held-out 0 (no case moves), tuning +46,
  experimental +4 (p 0.125), all TP53. Gate fails (needs a held-out gain).
  With specs/058 + 059 (`rtwin58s`): held-out +21, tuning +112, experimental
  +7 (p 0.016); fails only on TGF-beta −12. The entity-level variant
  (`suppress_only_entity`) fixes SMAD2/3 but loses MDM2 KO (experimental −17).
- `specs/061-own-source-gate/` — `DS_DEPLETION_OWN_PRODUCT=suppress_only_source`:
  an entity-own depleter is capped only while a source copy (never the
  target itself, amendment 1) is below baseline. With specs/058 + 059
  (`rtwin58g2`) **the registered gate passes**: held-out +32 (33 / 1),
  tuning +112, experimental +7 (7 / 0), worst pathway −1. Amendment 2 (two
  code-review fixes, `rtwin58g3`) reproduces it exactly. Post hoc; the
  independent korkut2015 check was inconclusive. **Adopted as the default
  2026-10-15** (Adam: "include something that passes more tests"), with
  specs/058's pins and specs/059's twins; Release 98 remains the replication.
  The pre-amendment rule failed (self-sourced cap locked TP53 low under MDM2 KO).
- `specs/062-korkut-test/` — korkut2015 (TEST) scored once against specs/061's
  candidate (`bench/korkut_bench.py`, `run_arm.sh --axis korkut[_predict]`):
  **inconclusive** (net +18, p 0.17; 3 drugs, 2 pathways; both arms below
  always-NOT). korkut2015 is now FIT. Release 98, the other check, is not
  published yet.
- `specs/063-sequestration/` — `DS_SEQUESTRATION`: a dead-end binding
  complex depletes its free partner (348 held-out misses, 274 inhibitor-like).
  No structural or curated signature separates sequestering complexes from
  terminal products (6,295 dead ends). **Refuted in three forms**:
  complex-to-partner edges feed the partner's own fall back (−1,180 /
  −1,686); sequesterer-to-partner edges (`sole_direct`, −1,121) read RUNX1's
  co-activators (SWI/SNF, PRC1, TAL1 complex) as sequesterers. The
  distinction is not encoded in Reactome; a curated flag would be.
- `specs/064-variant-duplicates/` — the generator emitted |members|
  identical copies of a reaction whose output had a bound and a shared slot
  (2,906 duplicate reaction nodes). Fixed in LNG #114; held-out −4 (one
  perturbation, traced), nothing else moves. Canonical build since
  2026-10-16: `20261016-0006_da08370`.
- `specs/065-sig-cache/` — the generator kept set-member signatures across
  pathways of one build (LNG #115 clears them per pathway; #116 five more
  derived caches, amendment 1). Measured both times: all 92
  networks identical in content, no predicted class changes; no rebuild.
- `specs/066-limit-cycles/` — the 23 held-out cases the iteration cap
  decided are Goodwin rings (traced in Fanconi: period-30 limit cycle; no
  damping, linear or log, converges). `DS_NONCONV_FALLBACK=lm` finds the
  fixed point by least squares where relaxation fails. **Adopted as the
  default** (Adam, 2026-10-16): held-out +1, experimental 0.
- `specs/067-specialization-bridges/` — Reactome curates one species under
  two ids: a generic complex of sets and a specific member combination of
  it, ordered only by `precedingEvent` (PDGF: autophosphorylation makes the
  specific p-receptor, every signalling step consumes the generic one).
  Generator `LNG_SPECIALIZATION_BRIDGE=add` joins the two copies with the
  same member leaves (61 pairs, 203 edges, 24 pathways). **Adopted as the
  default** (Adam, 2026-10-17): held-out +59 (101 / 42, p 8.7e-7),
  experimental 0, worst pathway −7; `replace` failed (DNA Damage Bypass
  −14). Open: TLS polymerases coupled through the generic template, cyclic
  growth in Hedgehog and TGF-beta. Its PLCG1 loss is fixed by specs/068.
- `specs/068-specialization-closed-loops/` — a specialization bridge is
  identity, not outside supply, when the generator decides whether a carrier
  loop is closed (`LNG_SPECIALIZATION_CLOSED=twin`, LNG #118). The bridges
  had reopened one loop, PLC-gamma's release of the generic PDGF receptor,
  so PLCG1 set the receptor's level. **Adopted as the default** (Adam,
  2026-10-17): held-out +18 (18 / 0, all PLCG1), experimental 0.
- `specs/071-pool-supply-reach/`, `specs/072-carrier-release-copies/` — two
  pool rules for the faithful build (`DS_POOL_SUPPLY_REACH=material`,
  `DS_CYCLE_CARRIER_RELEASE=reaction`), both no-ops on canonical; kept off by
  default until the faithful build is canonical (Adam, 2026-10-18).
- `specs/075-reeval-carriers/` — `DS_CYCLE_REEVAL_CARRIERS=carrier`: a
  carrier that another pool's drives re-evaluate (specs/039 amendment 5)
  reads by the carrier rule, not its whole reaction (USP7 read 50.5× in the
  TP53 pool's drive under MDM2 overexpression). **Adopted as the default**
  (Adam, 2026-10-18): canonical identical (0 pairs); faithful build TP53
  curator +252, experimental +15, but four TP53 perturbations no longer
  converge. `full` restores the old drives.
- `specs/073-pool-component-lm/` — `DS_NONCONV_POOLS=joint`: a non-converged
  component holding pool nodes is re-solved by LM with the pool update in
  the residual. **Adopted as the default** (Adam, 2026-10-18): canonical
  identical; faithful build TP53 non-converged 663 → 102, and MDM2/MDM4
  overexpression read UP at the fixed point (curators: DOWN).
- `specs/076-lm-restarts/` — `DS_NONCONV_RESTARTS=4`: the non-converged
  fallback retried from later iterates (LM's answer depended on where it
  started). **Adopted as the default** (Adam, 2026-10-19): canonical no class
  change, +167 converged; faithful build TP53 non-converged 306 → 0, curator
  −82 because the fixed point reads AKT1/AKT2/MDM4 knockout NORMAL as
  canonical does (specs/075's gains there were the cap's answer).
- `specs/077-open-carrier-release/` — `DS_CARRIER_LOOP_OPEN=fed_single`: an
  open carrier loop's free form with another producer stops reading its own
  release (RBBP8 released by MMEJ; the SCF(Skp2) ligase, PTK6's ERBB2
  heterodimers), unless a displacer drives the release. **Adopted as the
  default** (Adam, 2026-10-19): canonical held-out +65, experimental −5,
  worst pathway −2; faithful DSB +182. `fed` failed (apoptosis −26). Open: the
  generator's ubiquitin pass drains a ligase's own subunits (CKS1B), the
  specs/057 class.
- `specs/079-error-census/`, `specs/080-no-path-audit/` — analysis. Half of all
  errors have no route in the scored pathway's network. Curators agree with the
  experiments on only 82% of experimental cases, and 134 of 244 experimental errors
  give the curators' answer. Checked against Neo4j itself, 392 of 1,316 curator
  `no_path` errors have a route in the pathway's own reactions that our network loses
  (`db_reach.py`, `breaks.py`).
- `specs/081-recycle-edges/` — p130/p107 return to their E2F4 root complexes, which
  Reactome does not curate as a reaction. Generator `LNG_RECYCLE_POOLS` writes the cycle
  as a pool with several resting states plus a structural `recycle` edge.
  **Adopted** (Adam, 2026-10-20): tuning +70 (70 / 0, all Mitotic G1), experimental
  +8 (10 / 2), held-out 0. The p27 (CDKN1B) inversion is the curated assembly-factor role
  of p21/p27 in the CDK4/6 catalyst set.
- `specs/082-pten-mkrn1-loop/` — analysis: PTEN overexpression in PIP3 reads 83.5×
  UP because the curated AKT → MKRN1 ⊣ PTEN ⊣ PIP3 → AKT loop has gain 1 and rails;
  switch census (28 curator / 12 experimental wrong cases are rails).
- `specs/083-self-inhibitor-sum/` — `DS_SELF_ACTIVATOR_WEIGHT`: a positive regulator
  or catalyst built from its reaction's own input (TF:gene) counts it once (251
  reactions, 36 pathways). **Adopted** (Adam, 2026-10-20) on faithfulness: no class
  change. PTEN stays railed (two-copy mRNA, the MKRN1 loop).
- `specs/084-regenerated-complexes/` — `DS_CARRIER_COMPONENT_TWIN`: a closed complex
  the pathway only regenerates follows its components' starting copies. The first
  arm failed (PDGF −8: it fed copies a specialization bridge already supplies).
  Amendment 1 (only unsupplied forms) is **adopted** (Adam, 2026-10-21): ATM in Cell
  Cycle Checkpoints +6 curator, +4 experimental, nothing broken.
- `specs/085-silo-replication/` — silo bridges re-measured on today's networks:
  held-out −90, 73% of losses NORMAL turned into a change call. The positional silo
  ("two copies") is closed as an accuracy lever.
- `specs/086-released-substrate/` — a protein bound, modified and released by a
  complex, whose released copy carries on (β-catenin and the destruction complex),
  as a free/bound pool whose released copy reads the free state (generator
  `LNG_RELEASE_POOLS`). The first arm captured 1 of 20 variant chains; amendment 1
  (all chains) is **adopted** (Adam, 2026-10-21): WNT tuning +73 (APC, AMER1,
  TCF7L2), experimental +8, held-out 0.
- `specs/088-ub-not-displacer/` — `DS_CARRIER_LOOP_OPEN=fed_single_ub`: ubiquitin
  does not count as a release's displacer, so a ligase released after ubiquitinating its
  target (MDM2:MDM4) stops reading its own release. Gate passes, net 0 (held-out 0, tuning
  −4, experimental 0): MDM2 knockout fixed in Cell Cycle Checkpoints, Mitotic G1 −6 / −2
  through the SCF(Skp2) release exposing the CKS1B drain. Not adopted then (Adam,
  2026-10-22). Amendment 1, re-measured once specs/089 removed the drain: tuning +2,
  experimental +1 (both MDM2 knockout cases, ATM overexpression), held-out 0; **adopted**
  (Adam, 2026-10-23).
- `specs/089-ub-ligase-unit/` — the generator's ubiquitin pass no longer drains a
  ligase's own subunits (`LNG_UB_LIGASE_UNIT`, LNG #122): a leaf of the catalyst's curated
  `activeUnit`, or of a complex released unchanged, that leaves the reaction unmodified is not
  the substrate (CKS1B under SCF(Skp2)); the substrate's partners keep their drains, so CBL
  does not move (specs/078 lost 39 there). 73 of 727 depletion edges removed. **Adopted**
  (Adam, 2026-10-22) despite a failed gate: tuning +22 (p 0.0001), experimental +3, held-out
  −2 (p 0.63; SMAD7 right only through a removed drain, SHH where the arm reads GLI3R UP, a
  threshold pair). Not combined with specs/088.
- `specs/090-recycle-supply/` — `DS_POOL_SUPPLY_AND`: a recycle pool's AND-assembled
  resting complexes as one supply term, so an RBL1/RBL2 knockout empties its pool (review
  2026-10-22, finding 2). **Not adopted** (Adam, 2026-10-23): tuning −6, held-out 0,
  experimental 0; RBL2 overexpression reads UP because the p107 pool's RBL1 copy is
  DREAM-repressed transcription (R-HSA-8964525).
- `specs/091-base-weights/` — `DS_POOL_BASE_WEIGHT=split`: a pool's resting copies share
  one base weight, so magnitudes no longer depend on how many variant copies the generator
  emits (review 2026-10-22, finding 3). **Adopted** (Adam, 2026-10-23): tuning +35 (40 / 5,
  p < 1e-4, Mitotic G1 +37, all p107/p130 kinase/phosphatase perturbations), held-out 0,
  experimental +1. APC/AMER1 at cytosolic β-catenin were right only through the old 181×
  ceiling.
- `specs/092-solved-reach/` — analysis: the solver routes 2 of 1,100 curator `no_path`
  errors, and its carrier-loop drops cut 80 curator / 17 experimental wrong cases' routes
  (RB1, CBL in MET, SMAD7, PTPN11/SPRY2).
- `specs/093-recycling-release/` — `DS_CARRIER_LOOP_MOVES=keep`: an open release that moves
  the carrier across compartments (MET recycled to the membrane) is kept. **Not adopted**
  (Adam, 2026-10-23): held-out −10 (10 / 20): CBL's 10 MET cases fixed, but HGF and MET
  knockouts read UP because the kept recycling closes a self-sustaining receptor loop.
- `specs/094-experimental-census/` — the experimental errors: 220, of which 142 are
  curator–experiment disagreements and 78 agreed errors (PTEN in PIP3 the largest class).
  specs/014's loop elasticity re-measured and **refuted** (held-out −40, tuning −126,
  experimental −42).
- `specs/095-self-inhibitor-copies/` — `DS_SELF_INHIBITOR_COPIES=on`: the specs/022 pair
  test also accepts an inhibitor reached from another copy of the input's entity (PTEN's
  12 mRNA:miRNA complexes, built from the transcribed mRNA copy, shut off translation of the
  pinned root copy). **Adopted** (Adam, 2026-10-23): tuning +14, experimental +9 (9 / 0),
  held-out 0; all 9 PTEN overexpression experimental cases fixed.
- `specs/096-validation-plan/` — how accuracy will be validated: development continues on
  agreed errors and unhandled loop types with a stop rule; freeze; one external dataset scored
  once; blind curator pathways; justified disagreements recorded, scores kept raw.
- `specs/097-loop-audit/` — structural coverage: the 92 benchmark pathways match a 101-pathway
  sample (4.4% of nodes cyclic in both); unhandled loops cost ~8 points on routes through them;
  catalyst- and material-closed positive loops remain unhandled.
- `specs/098-free-ubiquitin/` — `DS_UBIQUITIN_MODE=inert`: free ubiquitin held at baseline.
  **Not adopted** (Adam, 2026-10-23): held-out +6 (USP1), experimental −4 (CCNE1's curated
  route through the p27 degradation cycle runs through free ubiquitin).
- `specs/100-transcribed-copy/` — `DS_RNA_TWIN=on`: a root copy of an RNA species reads the copies
  transcription produces (PTEN mRNA in PIP3: TP53 raised only the copy feeding the 12 miRNA
  repressors, inverting AKT). **Adopted** (Adam, 2026-10-24): tuning +28 (28 / 0), experimental
  +4 (4 / 0), held-out 0. Needs the generator's `rnas.csv` (LNG #125).
- `specs/101-agreed-gap-traces/` — agreed errors traced. T1: RB1 knockout in Mitotic G1 is not
  curated as a route (the repression lives only in a reaction name); a curation-review item, no
  rule. T2: the generator emitted set-to-set reactions as the full cross product (MDM2 alone made
  MDM2:MDM4; p-STAT6 made the STAT3 dimer); `LNG_BALANCE_CHOICES` drops copies whose protein
  content does not balance when better-balanced copies exist. **Adopted** (Adam, 2026-10-24):
  held-out +30 (34 / 4), tuning −12, experimental −4. The MDM4 knockout it was built for
  inverts: the freed MDM2 homodimer raises the ligase pool.
- `specs/102-summation-regulations/` — a curation supplement (generator `LNG_SUMMATION_REGULATIONS`,
  default on): three regulations Reactome states in a reaction's summation but not in its event
  graph, written as structural rules (RB1:E2F ⊣ E2F1/2/3-stimulated transcription; AKT's
  p-T-p21/p27 ⊣ Cip/Kip binding to cyclin:CDK; the p-Y-Cip/Kip CDK4/6 complex → CDK4/6 catalysis).
  Each regulator is one pool of its variants, added after enumeration. **Adopted** (Adam,
  2026-10-24): tuning +173, experimental +6 (641/849), held-out 0; mostly in-sample. Headlines are
  reported with it and with Release 97 alone. AKT1 inverts through p27's assembly-factor role.
- `specs/103-depletion-share/` — `DS_DEPLETION_COMBINE=share`: consumers of one pool share it
  (MDM4 knockout: the heterodimer alone de-repressed free MDM2 10×). **Not adopted** (gate fails, PIP3
  −13); experimental 0. specs/101 T4 records the agreed errors no structural rule reaches.
- `specs/104-pool-drive-twins/` — `DS_POOL_DRIVE_TWINS=hold`: a pool's drives do not read another
  copy of the pool's own entity. No class change (31 WNT curator values move). Not adopted (Adam, 2026-10-26, at the freeze).
- `specs/105-readout-entity-mean/` — bench `DS_KO_AGG=mean_entity`: a set readout read through proxy
  reaction copies or members weighs each curated entity once, not each variant copy (ERBB2: EGFR was
  9 of 225 copies). Gate passes: experimental +2 (643/849), tuning +6, held-out 0. **Adopted as the bench default**
  (Adam, 2026-10-26, at the freeze, specs/111); `DS_KO_AGG=mean` reproduces earlier numbers.
- `specs/106-loop-coverage/` — analysis: every loop in the network as generated and the rule that
  handles it. Rules touch 93% of benchmark loops and 88% in 101 other pathways (carrier rules most;
  amendment 1 counts cofactor holds); about 40% of loop nodes stay in unhandled remnants in both.
  Unhandled: positive loops closed by material steps (20 sample pathways after amendment 2, which counts a complex that is
  both input and catalyst as material), catalysts (7) or free ubiquitin; self-depleting rings only numerically. No further rules are chosen from experimental traces (Adam, 2026-10-25).
- `specs/107-loop-rules/` — three rules for loop types specs/106 left unhandled, judged by behaviour checks on networks
  without ground truth (convergence, rails, dose monotonicity, relabelling) plus a curator held-out no-harm guard.
  Rule 1 (`DS_METABOLITE_POOLS`): no class change, B4 −2. Rule 2 (`DS_CARRIER_MULTISTEP=strict`): held-out +6, all
  free ubiquitin (specs/098's USP1), B4 −2. Rule 3 (`wide`): rails 32% → 10% and monotone 89% → 97% on the benchmark,
  but held-out −50: the reference match joins SMAD2 with p-SMAD2, so a phosphorylation cycle reads as a closed carrier,
  and free ubiquitin closes in Mitotic G1. Amendment 1, rule 3b (`compartment`: compartment copies only, never free
  ubiquitin): behaviour checks pass, held-out −10, every case a Class I MHC knockout matched only through BAP31's railed
  loop. **None adopted** (Adam, 2026-10-26); defaults off; generator `leaf_refs.csv` is LNG #128, its `modifications` column LNG #129.
- `specs/108-justified-disagreements/` — the definition (Adam, 2026-10-25): a disagreement is justified only if traced
  to one element, matched by an expectation independent of the label, the old answer came through an artefact (a rail:
  the readout at a mild 0.5× / 2× dose reaches ≤ 0.01× or ≥ 99×), and the new answer fits; the same test flags gains as
  not credible; scores stay raw. Audit: 8.0% of correct held-out curator cases and 9.5% of correct experimental cases
  are rail-dependent (34 pathways). Of specs/107 rule 3b's ten MHC losses, one is justified.
- `specs/109-catalyst-loops/` — `DS_SELF_FED_MODE=catalyst_loop`: rule A (specs/040) only in positive loops closed through
  catalyst or regulator slots, holding only those edges; judged by the loop-rule gate (Adam, 2026-10-26: behaviour checks,
  un-railed cases keep direction, held-out losses read by specs/108). Gate passes: no class change on any axis, rails gone
  in small autocatalytic and self-activating loops (PAK1, CEBPA, HNF1B/ONECUT1), but large catalyst loops untouched
  (benchmark rails 61.0% → 60.5%): with several entries no node is self-fed. Amendment 1: a pair that is also an
  input (a complex modifying itself, 35% of catalyst pairs) is not a catalyst slot; same scores. The NHEJ and TGF-beta
  loops are material cycles of that kind. Not adopted (Adam, 2026-10-26, at the freeze: it changes no score).
- `specs/110-self-modification/` — `DS_DEDUP_ACTIVATORS=1` (specs/012) re-measured under the loop-rule gate: a complex
  that is input and catalyst of its own reaction counted once. Rails fall (benchmark 37% → 24%) but monotonicity
  falls too, and held-out −60 (3 / 63, 61 not rail-dependent: SOCS1 in IFN-gamma, BER glycosylases): the doubling
  carries correct moderate calls past the thresholds. Gate fails; not adopted.
- `specs/111-freeze/` — the freeze (Adam, 2026-10-26): specs/096's stop rule (four consecutive arms failing their
  gate, met in specs/107). specs/105 adopted, 104 and 109 not; generator and solver tagged `freeze-2026-10-26`; one
  canonical catalog; backlog for after the external test.
- `specs/112-hpn-dream/` — the external test: the frozen model on HPN-DREAM sub-challenge 1A (mTOR inhibitor descendants,
  Hill et al. 2016), networks of the 13 pathways acting on mTOR built with the frozen generator. Mean AUROC over the 30
  official contexts **0.693** (reachability 0.667, constant 0.5); above reachability in 29 of 30. 20 of 48 readouts are
  in a network. Among the 74 teams it would place 10th by mean AUROC (median 0.532; reachability alone 11th), and
  10th by the challenge's mean-rank rule (bootstrap over contexts 7th–13th; amendments 2–3).
  The dataset is now FIT.
- `specs/113-interaction-edges/` — Adam's idea of adding Reactome's IntAct interactions, tested as reachability on LINCS
  P100, scored once (the frozen model was not: 4 of 96 sites in its networks). Curated reach 0.499, plus one interaction
  hop at each end 0.514 (B − A +0.014, interval excludes 0), interactions alone 0.505. Randomly relabelled interactions
  gain most of it (B at the 86.5th percentile, short of the registered 95th). Not adopted; no rule. P100 is now FIT.
- `specs/114-whole-reactome/` — the frozen model on the 28 non-disease top-level pathways (build
  `20261027-1346_8b6527e_top114`, each network solved on its own), beside the frozen 92-pathway build, on FIT data.
  More cases are scored (korkut 11,062 → 17,836; decryptM 2,334 → 5,096), but the pattern does not change: decryptM
  regulated 187 vs reachability 301 of 725, korkut 0.723 below always-NOT 0.913, P100 reach 0.498.
  `run_arm.sh --axis decryptm --decryptm-cases F --decryptm-phospho F` runs the decryptM bench.
- `specs/074-hill-sat-cancellation/` — the `hill_sat` AND's smooth 100× cap
  cancelled below the cap: with ε = 1e-9, products between 87.5× and 99.8×
  were off by up to 8.9 and not monotone. The form is now chosen by the sign
  of raw − cap. **Adopted** (Adam, 2026-10-18): no class changes on either
  axis or build; converged curator cases 22,930 → 23,117 and experimental
  804 → 824 on canonical.
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
**Most of the test suite cannot fail.** Twenty-four files contain assertions; the
other seven execute code and print output.

| file | assertions | note |
|---|---|---|
| `test/test_config_validation.jl` | 219 | |
| `test/test_loop_elasticity.jl` | 175 | + 1 `@test_broken` |
| `test/test_propagator_invariants.jl` | 133 | specs/103 shared depletion added 2026-10-25 |
| `test/test_loop_pool.jl` | 81 | |
| `test/test_solver_determinism.jl` | 126 | |
| `test/test_and_curves.jl` | 96 | specs/074 cap band added 2026-10-18 |
| `test/test_scc_break_roles.jl` | 56 | |
| `test/test_self_inhibition.jl` | 115 | specs/022, added 2026-09-25; self-contained activators (specs/083) 2026-10-20 |
| `test/test_drugs.jl` | 134 | specs/032, added 2026-09-26; pathogens (specs/048) 2026-10-10; free ubiquitin (specs/098) 2026-10-23 |
| `test/test_set_pool.jl` | 79 | specs/033, added 2026-09-26; `variant_pool` (specs/046) 2026-10-09 |
| `test/test_conserved.jl` | 15 | specs/035, added 2026-09-26 |
| `test/test_cycles.jl` | 573 | specs/039, added 2026-09-27; drive twins (specs/104) 2026-10-25; shared pools (specs/049) 2026-10-11; supply reach (specs/071), carrier release copies (specs/072), joint pool components (specs/073), re-evaluated carriers (specs/075) 2026-10-18; several base states and recycle edges (specs/081) 2026-10-20 |
| `test/test_self_fed.jl` | 133 | specs/040, added 2026-09-28; catalyst loops (specs/109) 2026-10-26 |
| `test/test_carrier_loops.jl` | 244 | specs/049, added 2026-10-11; open loops (specs/077) 2026-10-19; component twins (specs/084) 2026-10-21 |
| `test/test_zero_supply.jl` | 43 | specs/052, added 2026-10-13 |
| `test/test_recycled_twin.jl` | 96 | specs/059–061, added 2026-10-14; RNA twins (specs/100) 2026-10-23 |
| `test/test_sequestration.jl` | 26 | specs/063, added 2026-10-15 |
| `test/test_limit_cycles.jl` | 56 | specs/066, added 2026-10-16; restarts (specs/076) 2026-10-19 |
| `test/test_material_loops.jl` | 54 | specs/107, added 2026-10-25; compartment copies (amendment 1) same day |
| `test/test_cli_observations.jl` | 50 | |
| `test/test_api_errors.jl` | 61 | |
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