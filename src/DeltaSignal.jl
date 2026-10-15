module DeltaSignal

using LinearAlgebra
using SparseArrays
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
include("core/cofactors.jl")
include("core/drugs.jl")
include("core/silo_bridges.jl")
include("core/reaction_model.jl")
include("core/conserved.jl")
include("core/cycles.jl")
include("core/self_fed.jl")
include("core/carrier_loops.jl")
include("core/recycled_twin.jl")
include("core/sequestration.jl")
include("core/zero_supply.jl")

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
export NetworkNode, SetExpansionMapping, ReactionNetwork, CyclePool
export LogicNetworkEdge, ReactionParams, SolverResult

# Parsing functions
export parse_logic_network, parse_uuid_mapping, parse_set_mappings, parse_complete_network, create_reaction_network
export parse_containment, network_containment_json, containment_from_json
export parse_drug_list, drug_stids_json, drug_stids_from_json, drug_mode, drug_uuids, pools_json, pools_from_json
export parse_pathogen_list, pathogen_stids_json, pathogen_stids_from_json, pathogen_mode, pathogen_uuids
export parse_carrier_loops, carrier_loops_json, carrier_loops_from_json, carrier_loop_mode, apply_carrier_loops, solved_carrier_network
export recycled_twin_mode, apply_recycled_twins
export sequestration_mode, apply_sequestration
export parse_pools, cycle_mode, cycle_carriers, PoolPath
export conserved_mode, conserved_uuids
export parse_node_leaves, node_leaves_json, node_leaves_from_json, zero_supply_mode, zero_supply_uuids
export convert_to_reaction_network, aggregate_to_pathway_view, validate_network_pathway_mapping
export compute_influence_scores, SteadyStateParams

# Solver functions
export solve_steady_state

# Server functions
export start_server

end # module DeltaSignal