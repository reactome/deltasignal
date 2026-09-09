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

### Biological Realism Enhancements
The system includes pathway-specific parameter tuning based on biological context:
- **Signaling pathways**: Fast dynamics, sensitive activation
- **Metabolic pathways**: Balanced kinetics, substrate availability
- **Transcriptional networks**: Slow dynamics, strong cooperativity
- **Stress response**: Rapid activation, strong feedback
- **Cell cycle**: Bistable switches, checkpoint dynamics

### Solver Configuration
The reaction model is driven by `DS_*` environment variables, resolved once per
solve into a `ReactionEvalConfig` (see `resolve_reaction_eval_config` in
`reaction_model.jl`). The **code defaults are the validated winning config**
(commit c9805b2) — env vars only override for benchmark sweeps:
- `DS_INHIBITION_MODE=divide`, `DS_AND_MODE=hill_log`, `DS_OR_MODE=mean`,
  `DS_ASSEMBLY_LIMITING=1`, `DS_HILL_LOG_ZMAX=10.0`
- SCC-condensation solve is on by default (`DS_SCC_SOLVE=1`); the legacy flat
  iteration and the `"fixed_point"` `SteadyStateParams.method` are not used by
  the CLI or API (both use the penalty/SCC path).
- Export aggregation default: `stoichiometry_weighted`.

### Testing Strategy
Tests are organized by functionality:
- Basic parsing and I/O validation
- Mathematical correctness of transformations
- Solver convergence and accuracy
- Biological realism validation
- Feedback loop handling
- Signal propagation through pathways

## Project Dependencies
- Julia 1.10+ (LTS required)
- Key packages: Optim.jl, JSON3.jl, DataFrames.jl, CSV.jl, HTTP.jl, ArgParse.jl
- Docker for containerized deployment
- UI: separate — Angular project in the WebsiteAngular workspace, consuming this repo's HTTP API (`docs/API.md`)