#!/usr/bin/env julia

# Focused test on inhibition effectiveness
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

function test_inhibition_focused()
    println("🚫 Testing Inhibition Effectiveness")
    println("=" ^ 40)
    
    try
        # Load the propagation test network
        logic_path = joinpath(@__DIR__, "..", "examples", "propagation_test_network.tsv")
        uuid_path = joinpath(@__DIR__, "..", "examples", "propagation_test_uuid_mapping.tsv")
        
        network = parse_complete_network(logic_path, uuid_path, nothing)
        
        # Focus on inhibition scenarios
        test_cases = [
            ("Inhibitor only", Dict("ROOT6" => 80.0)),
            ("Activator only", Dict("ROOT5" => 80.0)), 
            ("Balanced competition", Dict("ROOT5" => 60.0, "ROOT6" => 60.0)),
            ("Strong inhibitor", Dict("ROOT5" => 40.0, "ROOT6" => 90.0)),
            ("Weak inhibitor", Dict("ROOT5" => 90.0, "ROOT6" => 40.0))
        ]
        
        println("Testing ROOT5(+) and ROOT6(-) → G1 → H1 → TERM5 pathway")
        println("Baseline expectation: ~20%\n")
        
        for (case_name, perturbations) in test_cases
            observations = Dict{String, Tuple{Float64, Float64}}()
            for (root, activity) in perturbations
                observations[root] = (activity, 0.9)
            end
            
            params = SteadyStateParams(1.0, 0.05, 150, 1e-5, "penalty")
            result = solve_steady_state(network, observations, params)
            
            # Extract key activities
            g1_activity = get(result.node_activities, "G1", 0.2) * 100
            h1_activity = get(result.node_activities, "H1", 0.2) * 100
            term5_activity = get(result.node_activities, "TERM5", 0.2) * 100
            
            println("$case_name:")
            for (root, activity) in perturbations
                println("  $root: $(activity)%")
            end
            println("  → G1: $(round(g1_activity, digits=1))%")
            println("  → H1: $(round(h1_activity, digits=1))%")
            println("  → TERM5: $(round(term5_activity, digits=1))%")
            
            # Analysis
            if case_name == "Inhibitor only"
                if term5_activity < 25.0
                    println("  ✅ Inhibition effective")
                else
                    println("  ⚠️  Inhibition weak (expected < 25%)")
                end
            elseif case_name == "Strong inhibitor"
                if term5_activity < 30.0
                    println("  ✅ Strong inhibition effective")
                else
                    println("  ⚠️  Strong inhibition insufficient")
                end
            end
            
            println()
        end
        
        println("🔍 Testing negative feedback loop:")
        println("ROOT10 → Q1 → R1 → S1 → [TERM8, Q1(-)]")
        
        neg_fb_cases = [
            ("Neg FB up", Dict("ROOT10" => 80.0)),
            ("Neg FB down", Dict("ROOT10" => 20.0))
        ]
        
        for (case_name, perturbations) in neg_fb_cases
            observations = Dict{String, Tuple{Float64, Float64}}()
            for (root, activity) in perturbations
                observations[root] = (activity, 0.9)
            end
            
            params = SteadyStateParams(1.0, 0.05, 150, 1e-5, "penalty")
            result = solve_steady_state(network, observations, params)
            
            root10_input = perturbations["ROOT10"]
            q1_activity = get(result.node_activities, "Q1", 0.2) * 100
            s1_activity = get(result.node_activities, "S1", 0.2) * 100  
            term8_activity = get(result.node_activities, "TERM8", 0.2) * 100
            
            println("$case_name (ROOT10: $(root10_input)%):")
            println("  → Q1: $(round(q1_activity, digits=1))%")
            println("  → S1: $(round(s1_activity, digits=1))% (inhibits Q1)")
            println("  → TERM8: $(round(term8_activity, digits=1))%")
            
            # Check if feedback is working
            if abs(term8_activity - 20.0) < abs(root10_input - 20.0) * 0.9
                println("  ✅ Negative feedback dampening detected")
            else
                println("  ⚠️  Negative feedback may be weak")
            end
            
            println()
        end
        
        return true
        
    catch e
        println("❌ Test failed: $e")
        return false
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_inhibition_focused()
end