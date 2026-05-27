#!/usr/bin/env julia

# Test time-dynamic simulation
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
include("../src/solvers/time_dynamic.jl")

function test_time_dynamic_simulation()
    println("🕐 Testing Time-Dynamic Simulation")
    println("=" ^ 40)
    
    try
        # Load a simple network for testing
        logic_path = joinpath(@__DIR__, "..", "examples", "propagation_test_network.tsv")
        uuid_path = joinpath(@__DIR__, "..", "examples", "propagation_test_uuid_mapping.tsv")
        
        network = parse_complete_network(logic_path, uuid_path, nothing)
        
        println("📁 Network loaded with $(length(network.nodes)) nodes, $(length(network.edges)) edges")
        
        # Test scenarios
        test_cases = [
            ("Pulse response", Dict("ROOT1" => (80.0, 0.9)), 5.0),
            ("Step response", Dict("ROOT10" => (60.0, 0.9)), 10.0),
            ("Multi-input", Dict("ROOT5" => (70.0, 0.9), "ROOT6" => (30.0, 0.9)), 8.0)
        ]
        
        for (case_name, perturbations, t_end) in test_cases
            println("\\n" * "="^50)
            println("TEST CASE: $case_name")
            println("="^50)
            
            # Set up time-dynamic parameters
            params = TimeDynamicParams(
                0.0,         # t_start
                t_end,       # t_end  
                0.1,         # dt
                1e-5,        # tolerance
                1000,        # max_iters
                0.05,        # substrate_consumption_rate
                0.005,       # substrate_replenishment_rate
                0.08,        # product_accumulation_rate
                0.02,        # product_decay_rate
                0.3          # observation_weight
            )
            
            println("🎯 Perturbations:")
            for (node, (activity, _)) in perturbations
                println("  $node: $(activity)%")
            end
            println("⏱️  Time span: $(params.t_start) → $(params.t_end) s (dt=$(params.dt))")
            
            # Run time-dynamic simulation
            result = rollout_time_dynamic(network, perturbations, params)
            
            println("\\n🧮 Simulation results:")
            println("  Converged: $(result.converged)")
            println("  Final time: $(round(result.final_time, digits=2))s")
            println("  Solve time: $(round(result.solve_time, digits=3))s")
            println("  Time points: $(length(result.time_points))")
            
            # Analyze key nodes at different time points
            key_times = [0.0, t_end/4, t_end/2, t_end]
            key_nodes = ["ROOT1", "A1", "TERM1", "ROOT10", "Q1", "TERM8"]
            
            println("\\n📊 Activity time series (selected nodes):")
            println("Time\\t", join(key_nodes[1:min(4, length(key_nodes))], "\\t"))
            
            for t in key_times
                t_idx = argmin(abs.(result.time_points .- t))
                activities = String[]
                
                for node in key_nodes[1:min(4, length(key_nodes))]
                    if haskey(result.node_activities, node)
                        activity = result.node_activities[node][t_idx] * 100
                        push!(activities, "$(round(activity, digits=1))%")
                    else
                        push!(activities, "N/A")
                    end
                end
                
                println("$(round(t, digits=1))\\t", join(activities, "\\t"))
            end
            
            # Check substrate/product dynamics if present
            if !isempty(result.substrate_levels)
                println("\\n⚗️  Substrate dynamics:")
                for (substrate, timeseries) in result.substrate_levels
                    initial = timeseries[1]
                    final = timeseries[end]
                    change = final - initial
                    println("  $substrate: $(round(initial, digits=2)) → $(round(final, digits=2)) (Δ=$(round(change, digits=2)))")
                end
            end
            
            if !isempty(result.product_levels)
                println("\\n🧪 Product dynamics:")
                for (product, timeseries) in result.product_levels
                    initial = timeseries[1]
                    final = timeseries[end]
                    accumulation = final - initial
                    println("  $product: $(round(initial, digits=2)) → $(round(final, digits=2)) (+$(round(accumulation, digits=2)))")
                end
            end
            
            # Validate expected behaviors
            println("\\n✅ Validation:")
            
            # Check that perturbed roots maintain their activity levels
            all_valid = true
            for (root, (target_activity, _)) in perturbations
                if haskey(result.node_activities, root)
                    final_activity = result.node_activities[root][end] * 100
                    if abs(final_activity - target_activity) < 15.0  # 15% tolerance
                        println("  ✅ $root maintained target activity: $(round(final_activity, digits=1))% (target: $(target_activity)%)")
                    else
                        println("  ⚠️  $root deviated from target: $(round(final_activity, digits=1))% (target: $(target_activity)%)")
                        all_valid = false
                    end
                end
            end
            
            # Check that some downstream propagation occurred
            downstream_changed = false
            for (node, timeseries) in result.node_activities
                if !(node in keys(perturbations))  # Not a directly perturbed node
                    initial = timeseries[1] * 100
                    final = timeseries[end] * 100
                    if abs(final - initial) > 5.0  # 5% change threshold
                        downstream_changed = true
                        break
                    end
                end
            end
            
            if downstream_changed
                println("  ✅ Downstream propagation detected")
            else
                println("  ⚠️  Limited downstream propagation")
                all_valid = false
            end
            
            if all_valid
                println("  🎉 Test case passed!")
            else
                println("  ⚠️  Test case has issues")
            end
        end
        
        # Compare with steady-state solver
        println("\\n" * "="^60)
        println("🔄 COMPARISON: Time-Dynamic vs Steady-State")
        println("="^60)
        
        test_obs = Dict("ROOT1" => (70.0, 0.9), "ROOT2" => (30.0, 0.9))
        
        # Steady-state result
        ss_params = SteadyStateParams(1.0, 0.05, 150, 1e-5, "fixed_point")
        ss_result = solve_steady_state(network, test_obs, ss_params)
        
        # Time-dynamic result (final state)
        td_params = TimeDynamicParams(0.0, 10.0, 0.1)
        td_result = rollout_time_dynamic(network, test_obs, td_params)
        
        println("\\nFinal state comparison (selected nodes):")
        println("Node\\tSteady-State\\tTime-Dynamic\\tDifference")
        
        comparison_nodes = ["ROOT1", "ROOT2", "A1", "TERM1"]
        for node in comparison_nodes
            if haskey(ss_result.node_activities, node) && haskey(td_result.node_activities, node)
                ss_val = ss_result.node_activities[node] * 100
                td_val = td_result.node_activities[node][end] * 100
                diff = abs(td_val - ss_val)
                
                println("$node\\t$(round(ss_val, digits=1))%\\t\\t$(round(td_val, digits=1))%\\t\\t$(round(diff, digits=1))%")
            end
        end
        
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
    test_time_dynamic_simulation()
end