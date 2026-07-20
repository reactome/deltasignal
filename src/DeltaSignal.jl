module DeltaSignal

using LinearAlgebra
using Statistics
using Optim
using ForwardDiff
using JSON3
using CSV
using DataFrames
using HTTP
using Random

# Input/Output handling (must be first for data structures)
include("io/tsv_parser.jl")
include("io/reactome_mapper.jl")

# Core mathematical operations
include("core/sensitivity.jl")
include("core/aggregators.jl")
include("core/hill_functions.jl")
include("core/reaction_model.jl")

# Solvers
include("solvers/steady_state.jl")

# API server
include("api/server.jl")

# NOTE: the "biological realism" module cluster (biological_realism,
# compartmentalization, experimental_constraints, temporal_dynamics,
# stochastic_effects, feedback_enhancements), the enhanced_steady_state /
# time_dynamic solvers, and parameter_learning were moved to attic/ (2026-07).
# A call-graph audit confirmed none are reachable from the live parse/solve
# path (CLI + API) — they were only ever exercised by their own tests. See
# attic/README.md. Revive from there if/when the learning or time-dynamic
# modes are actually wired to an entry point.

# Core data structures
export NetworkNode, SetExpansionMapping, ReactionNetwork
export LogicNetworkEdge, ReactionParams, SolverResult

# Parsing functions
export parse_logic_network, parse_uuid_mapping, parse_set_mappings, parse_complete_network, create_reaction_network
export convert_to_reaction_network, aggregate_to_pathway_view, validate_network_pathway_mapping
export compute_influence_scores, SteadyStateParams

# Solver functions
export solve_steady_state

# Server functions
export start_server

end # module DeltaSignal