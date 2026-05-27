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
include("core/feedback_enhancements.jl")
include("core/biological_realism.jl")
include("core/compartmentalization.jl")
include("core/experimental_constraints.jl")
include("core/temporal_dynamics.jl")
include("core/stochastic_effects.jl")

# Solvers
include("solvers/steady_state.jl")
include("solvers/enhanced_steady_state.jl")
include("solvers/time_dynamic.jl")

# Parameter learning
include("learning/parameter_learning.jl")

# API server
include("api/server.jl")

# Core data structures
export NetworkNode, SetExpansionMapping, ReactionNetwork
export LogicNetworkEdge, ReactionParams, SolverResult
export TimeDynamicParams, TimeDynamicResult

# Parsing functions
export parse_logic_network, parse_uuid_mapping, parse_set_mappings, parse_complete_network, create_reaction_network
export convert_to_reaction_network, aggregate_to_pathway_view, validate_network_pathway_mapping
export compute_influence_scores, SteadyStateParams

# Solver functions
export solve_steady_state, solve_steady_state_enhanced, rollout_time_dynamic
export learn_parameters

# Biological realism functions
export PathwayType, SIGNALING, METABOLIC, TRANSCRIPTIONAL, STRESS_RESPONSE, CELL_CYCLE, APOPTOSIS
export create_biologically_enhanced_reaction_network, compute_enhanced_reaction_output
export analyze_pathway_parameters, infer_pathway_type

# Compartmentalization functions
export Compartment, CYTOPLASM, NUCLEUS, MITOCHONDRIA, ENDOPLASMIC_RETICULUM, GOLGI, LYSOSOME
export create_compartmentalized_reactions, compute_compartmentalized_reaction_output
export analyze_compartment_effects

# Experimental constraints functions
export ExperimentType, ExperimentalDataset, create_experimental_constraints
export create_experimentally_constrained_network, load_experimental_dataset

# Temporal dynamics functions
export ProcessType, ENZYME_CATALYSIS, PROTEIN_BINDING, TRANSPORT, PROTEIN_SYNTHESIS
export simulate_realistic_time_dynamics, analyze_temporal_characteristics, infer_biological_processes

# Stochastic effects functions  
export NoiseType, TRANSCRIPTIONAL_BURSTING, POISSON_NOISE, THERMAL_NOISE
export simulate_with_stochastic_effects, analyze_cell_variability

# Server functions
export start_server

end # module DeltaSignal