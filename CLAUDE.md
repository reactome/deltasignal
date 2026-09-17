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
  --output /app/data/output/results.json \
  --mu 1.0 \
  --gamma 0.1

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
  Confidence is 0-1. Any extra column — including `condition` — is ignored, so
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
- `DS_INHIBITION_MODE=divide`, `DS_AND_MODE=hill_log`, `DS_OR_MODE=mean`,
  `DS_ASSEMBLY_LIMITING=1`, `DS_INHIBITOR_EPS=1e-12`, `DS_HILL_SAT_EPS=1e-5`,
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
- `DS_AND_MODE=hill_log` is an **empirical** default, not a principled one.
  `hill_sat` implements the stated AND intent more faithfully (10x10 = 99.98 vs
  74.06) and still scores 90 cases worse. Why compression helps is unexplained.
- SCC-condensation solve is on by default (`DS_SCC_SOLVE=1`); the legacy flat
  iteration and the `"fixed_point"` `SteadyStateParams.method` are not used by
  the CLI or API (both use the penalty/SCC path).
- `DS_COFACTOR_MODE=inert` — metabolic cofactors (ATP, ADP, NAD+, H2O, Pi …)
  are pinned at baseline, so they remain AND inputs but cannot carry a
  perturbation. The networks still contain them: the generator represents what
  curators recorded, and this is the propagator declining to route through it.
  Set `propagate` for the previous behaviour. The list and the measured A/B are
  in `src/core/cofactors.jl`.
- Export aggregation default: `stoichiometry_weighted`.

### Where design decisions live

This repo uses spec-kit. **`specs/NNN-name/` is the record of why**; CLAUDE.md
carries only current state. Before changing solver behaviour, read the spec
for the last feature that touched it rather than re-deriving from the code.

- `specs/002-upregulation-propagation/` — the AND/assembly defaults above.
  `research.md` holds the measured A/B including the negative results.
- `specs/003-solver-objective/` — the solver runs a damped fixed-point
  iteration, not the specified minimisation; `mu` and `gamma` are reported
  but read by nothing. Open.
- `.specify/memory/constitution.md` — project principles the specs are
  checked against.

Per-feature numbers belong in that feature's `research.md`, not here.

### Testing Strategy
**Most of the test suite cannot fail.** Only three files contain assertions:

| file | `@test`s |
|---|---|
| `test/test_config_validation.jl` | 156 |
| `test/test_and_curves.jl` | 17 |
| `test/test_cycle_handling.jl` | 29 + 2 `@test_broken` |
| `test/test_worked_example.jl` | 9 |

Counts are CI-verified (`.github/workflows/test.yml` runs these four by name
and prints each `Test Summary`), not `@test` occurrences — several
testsets generate assertions in loops. The earlier figures in this table (42 /
22 / 8) were wrong in both directions.

**A failing testset used to hide every later one.** A top-level `@testset`
throws when it finishes with a failure, which aborts the file. In the dev
container, where `DS_*` overrides contradict the code defaults, that meant
`test_config_validation.jl` ran 10 assertions and silently skipped **141** —
the cofactor tests, the silo-bridge tests and the observation-membership test
all looked green because they never executed. The files now nest their testsets
inside one outer testset so failures accumulate and everything runs.
| the other seven | **0** |

`test_basic.jl`, `test_steady_state.jl`, `test_hill_function.jl`,
`test_feedback_loops.jl`, `test_inhibition_focused.jl`,
`test_pathway_propagation.jl` and `test_biological_fixes.jl` execute code and
print output; they pass whatever the code does. Treat a green run of those as
"it did not throw", nothing more.

The documented docker runner is also broken: `docker-compose.dev.yml`'s
`test-runner` chains `test/test_full_pipeline.jl`, which does not exist.

When adding behaviour, add assertions to one of the three real files or start
a new one — do not extend a file from the zero-assertion list and assume it is
covering anything.

## Project Dependencies
- Julia 1.10+ (LTS required)
- Key packages: Optim.jl, JSON3.jl, DataFrames.jl, CSV.jl, HTTP.jl, ArgParse.jl
- Docker for containerized deployment
- UI: separate — Angular project in the WebsiteAngular workspace, consuming this repo's HTTP API (`docs/API.md`)