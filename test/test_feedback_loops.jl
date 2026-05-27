#!/usr/bin/env julia

# Test feedback loop behavior
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

function test_feedback_loops()
    println("🔄 Testing Feedback Loop Behavior")
    println("=" ^ 50)
    
    try
        # Load feedback test network
        logic_path = joinpath(@__DIR__, "..", "examples", "feedback_test_network.tsv")
        uuid_path = joinpath(@__DIR__, "..", "examples", "feedback_test_uuid_mapping.tsv")
        
        println("📁 Loading feedback test network...")
        network = parse_complete_network(logic_path, uuid_path, nothing)
        
        println("Network structure:")
        println("  - $(length(network.nodes)) nodes")
        println("  - $(length(network.edges)) edges")
        
        # Analyze the network structure
        println("\n🔍 Network topology analysis:")
        for edge in network.edges
            effect = edge.is_positive ? "activates" : "inhibits"
            println("  $(edge.parent_uuid) → $effect → $(edge.child_uuid)")
        end
        
        # Identify feedback loops
        println("\n🔄 Feedback loop analysis:")
        println("  Positive feedback: A→B→C→A (should amplify)")
        println("  Negative feedback: D→E→F⊣D (should stabilize/oscillate)")
        println("  Positive feedback: G⇄H (should amplify)")
        println("  Negative feedback: I→J→K⊣I (should stabilize)")
        
        # Test with realistic biological perturbation levels
        # Using internal scale: 1.0 = baseline (1x), 2.0 = 2x normal, 5.0 = 5x normal
        test_cases = [
            ("1.5x normal input", Dict{String, Tuple{Float64, Float64}}(
                "A" => (1.5, 0.9), "D" => (1.5, 0.9), 
                "G" => (1.5, 0.9), "I" => (1.5, 0.9)
            )),
            ("3x normal input", Dict{String, Tuple{Float64, Float64}}(
                "A" => (3.0, 0.9), "D" => (3.0, 0.9), 
                "G" => (3.0, 0.9), "I" => (3.0, 0.9)
            )),
            ("5x normal input", Dict{String, Tuple{Float64, Float64}}(
                "A" => (5.0, 0.9), "D" => (5.0, 0.9), 
                "G" => (5.0, 0.9), "I" => (5.0, 0.9)
            ))
        ]
        
        for (case_name, observations) in test_cases
            println("\n" * "="^30)
            println("TEST CASE: $case_name")
            println("="^30)
            
            # Solve with fixed-point method (better for feedback loops)
            params = SteadyStateParams(1.0, 0.05, 200, 1e-5, "fixed_point")  # Lower gamma for less baseline pull
            result = solve_steady_state(network, observations, params)
            
            println("🧮 Solver results:")
            println("  Converged: $(result.converged)")
            println("  Iterations: $(result.iterations)")
            println("  Final residual: $(round(result.final_residual, digits=6))")
            
            println("\n📊 Node activities:")
            
            # Group by feedback loops
            pos_loop1 = ["A", "B", "C"]  # Positive feedback A→B→C→A
            neg_loop1 = ["D", "E", "F"]  # Negative feedback D→E→F⊣D  
            pos_loop2 = ["G", "H"]       # Positive feedback G⇄H
            neg_loop2 = ["I", "J", "K"]  # Negative feedback I→J→K⊣I
            
            println("  Positive Loop A→B→C→A:")
            for node in pos_loop1
                activity = get(result.node_activities, node, 0.0) * 100
                obs_str = haskey(observations, node) ? " [OBS: $(observations[node][1])%]" : ""
                println("    $node: $(round(activity, digits=1))%$obs_str")
            end
            
            println("  Negative Loop D→E→F⊣D:")
            for node in neg_loop1  
                activity = get(result.node_activities, node, 0.0) * 100
                obs_str = haskey(observations, node) ? " [OBS: $(observations[node][1])%]" : ""
                println("    $node: $(round(activity, digits=1))%$obs_str")
            end
            
            println("  Positive Loop G⇄H:")
            for node in pos_loop2
                activity = get(result.node_activities, node, 0.0) * 100
                obs_str = haskey(observations, node) ? " [OBS: $(observations[node][1])%]" : ""
                println("    $node: $(round(activity, digits=1))%$obs_str")
            end
            
            println("  Negative Loop I→J→K⊣I:")
            for node in neg_loop2
                activity = get(result.node_activities, node, 0.0) * 100
                obs_str = haskey(observations, node) ? " [OBS: $(observations[node][1])%]" : ""
                println("    $node: $(round(activity, digits=1))%$obs_str")
            end
            
            # Analyze feedback behavior
            println("\n🔬 Feedback analysis:")
            
            # Check if positive loops amplify
            a_activity = get(result.node_activities, "A", 0.0) * 100
            b_activity = get(result.node_activities, "B", 0.0) * 100
            c_activity = get(result.node_activities, "C", 0.0) * 100
            pos_loop1_avg = (a_activity + b_activity + c_activity) / 3
            
            g_activity = get(result.node_activities, "G", 0.0) * 100
            h_activity = get(result.node_activities, "H", 0.0) * 100
            pos_loop2_avg = (g_activity + h_activity) / 2
            
            # Check if negative loops stabilize/dampen
            d_activity = get(result.node_activities, "D", 0.0) * 100
            e_activity = get(result.node_activities, "E", 0.0) * 100
            f_activity = get(result.node_activities, "F", 0.0) * 100
            neg_loop1_avg = (d_activity + e_activity + f_activity) / 3
            
            i_activity = get(result.node_activities, "I", 0.0) * 100
            j_activity = get(result.node_activities, "J", 0.0) * 100
            k_activity = get(result.node_activities, "K", 0.0) * 100
            neg_loop2_avg = (i_activity + j_activity + k_activity) / 3
            
            input_level = first(values(observations))[1]  # All observations are same in each test
            
            println("  Input level: $(input_level)x normal")
            println("  Pos Loop 1 avg: $(round(pos_loop1_avg, digits=1))x (should be > $(input_level)x)")
            println("  Pos Loop 2 avg: $(round(pos_loop2_avg, digits=1))x (should be > $(input_level)x)")
            println("  Neg Loop 1 avg: $(round(neg_loop1_avg, digits=1))x (should be ≤ $(input_level)x)")
            println("  Neg Loop 2 avg: $(round(neg_loop2_avg, digits=1))x (should be ≤ $(input_level)x)")
            
            # Check for problematic zero values
            all_activities = [get(result.node_activities, node, 0.0) * 100 for node in keys(network.nodes)]
            zero_count = count(x -> x < 0.1, all_activities)  # Essentially zero
            
            if zero_count > 0
                println("  ⚠️  WARNING: $zero_count nodes near zero - may indicate feedback issues")
            end
        end
        
        println("\n" * "="^50)
        println("🎯 FEEDBACK BEHAVIOR ASSESSMENT")
        println("="^50)
        
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
    test_feedback_loops()
end