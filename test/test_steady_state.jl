#!/usr/bin/env julia

# Test steady-state solver
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV
using DataFrames
using JSON3

# Add the src directory to the load path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
include("../src/io/tsv_parser.jl")
include("../src/core/sensitivity.jl")
include("../src/core/aggregators.jl")
include("../src/core/hill_functions.jl")
include("../src/core/reaction_model.jl")
include("../src/solvers/steady_state.jl")

function test_steady_state_solver()
    println("🧮 Testing DeltaSignal steady-state solver...")
    
    try
        # Load network from examples
        logic_path = joinpath(@__DIR__, "..", "examples", "sample_logic_network.tsv")
        uuid_path = joinpath(@__DIR__, "..", "examples", "sample_uuid_mapping.tsv")
        set_path = joinpath(@__DIR__, "..", "examples", "sample_set_mappings.tsv")
        
        println("📁 Loading network...")
        network = parse_complete_network(logic_path, uuid_path, set_path)
        
        # Set up observations (convert from sample CSV format)
        observations = Dict{String, Tuple{Float64, Float64}}(
            "parent-001" => (75.0, 0.9),  # activity=75%, confidence=0.9
            "parent-003" => (50.0, 0.8),  # activity=50%, confidence=0.8
            "child-002" => (25.0, 0.7)    # activity=25%, confidence=0.7
        )
        
        println("🎯 Observations:")
        for (uuid, (activity, confidence)) in observations
            println("  $uuid: $(activity)% (confidence: $confidence)")
        end
        
        # Test penalty solver
        println("\n⚖️  Testing penalty solver...")
        params_penalty = SteadyStateParams(1.0, 0.1, 200, 1e-4, "penalty")
        result_penalty = solve_steady_state(network, observations, params_penalty)
        
        println("Penalty method results:")
        println("  Converged: $(result_penalty.converged)")
        println("  Iterations: $(result_penalty.iterations)")
        println("  Solve time: $(round(result_penalty.solve_time, digits=3))s")
        
        # Display results
        println("\n📊 Final node activities (penalty method):")
        for (uuid, activity) in sort(collect(result_penalty.node_activities))
            activity_ui = round(activity * 100, digits=1)  # Convert back to UI scale
            
            # Mark if this was observed
            obs_marker = haskey(observations, uuid) ? " [OBS]" : ""
            println("  $uuid: $(activity_ui)%$obs_marker")
        end
        
        # Test influence scores
        println("\n🎯 Computing influence scores...")
        reactions = convert_to_reaction_network(network)
        influence_scores = compute_influence_scores(result_penalty, reactions)
        
        if !isempty(influence_scores)
            println("Top influential nodes:")
            sorted_influence = sort(collect(influence_scores), by=x->x[2], rev=true)
            for (i, (uuid, score)) in enumerate(sorted_influence[1:min(5, length(sorted_influence))])
                println("  $(i). $uuid: $(round(score, digits=4))")
            end
        end
        
        println("\n✨ Steady-state solver test completed successfully!")
        return true
        
    catch e
        println("❌ Test failed: $e")
        println("Stack trace:")
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
        return false
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_steady_state_solver()
end