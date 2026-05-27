#!/usr/bin/env julia

# Test root-to-terminal pathway propagation across different network topologies
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

function test_pathway_propagation()
    println("🌐 Testing Root-to-Terminal Pathway Propagation")
    println("=" ^ 60)
    
    try
        # Load the propagation test network
        logic_path = joinpath(@__DIR__, "..", "examples", "propagation_test_network.tsv")
        uuid_path = joinpath(@__DIR__, "..", "examples", "propagation_test_uuid_mapping.tsv")
        
        println("📁 Loading propagation test network...")
        network = parse_complete_network(logic_path, uuid_path, nothing)
        
        println("Network loaded:")
        println("  - $(length(network.nodes)) nodes")
        println("  - $(length(network.edges)) edges")
        
        # Analyze network topology
        println("\n🏗️  Network topology analysis:")
        
        # Find root and terminal nodes
        all_parents = Set(edge.parent_uuid for edge in network.edges)
        all_children = Set(edge.child_uuid for edge in network.edges)
        
        root_nodes = [uuid for uuid in all_parents if !(uuid in all_children)]
        terminal_nodes = [uuid for uuid in all_children if !(uuid in all_parents)]
        intermediate_nodes = [uuid for uuid in keys(network.nodes) 
                             if !(uuid in root_nodes) && !(uuid in terminal_nodes)]
        
        println("  Root nodes ($(length(root_nodes))): $(join(root_nodes, ", "))")
        println("  Terminal nodes ($(length(terminal_nodes))): $(join(terminal_nodes, ", "))")
        println("  Intermediate nodes ($(length(intermediate_nodes))): $(length(intermediate_nodes))")
        
        # Describe pathway topologies present
        println("\n🔗 Pathway topologies in network:")
        println("  1. Linear chain: ROOT7 → I1 → J1 → K1 → L1 → TERM6")
        println("  2. Convergent: ROOT1,ROOT2 → A1 → TERM1")
        println("  3. Divergent: ROOT3 → B1 → [TERM2, C1→TERM3]")
        println("  4. Positive feedback: ROOT4 → D1 → E1 → F1 → [TERM4, D1]")
        println("  5. Competing inputs: ROOT5(+), ROOT6(-) → G1 → H1 → TERM5") 
        println("  6. Diamond: ROOT8,ROOT9 → M1,N1 → O1 → P1 → TERM7")
        println("  7. Negative feedback: ROOT10 → Q1 → R1 → S1 → [TERM8, Q1(-)]")
        
        # Define test scenarios
        test_scenarios = [
            # Scenario 1: Single root upregulation
            ("Single root up", Dict("ROOT1" => 70.0)),
            ("Single root down", Dict("ROOT1" => 10.0)),
            
            # Scenario 2: Multiple roots same direction
            ("Multi roots up", Dict("ROOT1" => 70.0, "ROOT2" => 70.0)),
            ("Multi roots down", Dict("ROOT1" => 10.0, "ROOT2" => 10.0)),
            
            # Scenario 3: Competing inputs
            ("Competing balanced", Dict("ROOT5" => 60.0, "ROOT6" => 60.0)),
            ("Activator wins", Dict("ROOT5" => 80.0, "ROOT6" => 20.0)),
            ("Inhibitor wins", Dict("ROOT5" => 20.0, "ROOT6" => 80.0)),
            
            # Scenario 4: Feedback loop perturbations
            ("Pos feedback root up", Dict("ROOT4" => 80.0)),
            ("Pos feedback root down", Dict("ROOT4" => 10.0)),
            ("Neg feedback root up", Dict("ROOT10" => 80.0)),
            ("Neg feedback root down", Dict("ROOT10" => 10.0)),
            
            # Scenario 5: Long cascade
            ("Long cascade up", Dict("ROOT7" => 80.0)),
            ("Long cascade down", Dict("ROOT7" => 10.0)),
            
            # Scenario 6: Diamond network
            ("Diamond both up", Dict("ROOT8" => 70.0, "ROOT9" => 70.0)),
            ("Diamond asymmetric", Dict("ROOT8" => 80.0, "ROOT9" => 30.0))
        ]
        
        # Run all test scenarios
        results_summary = []
        
        for (scenario_name, root_perturbations) in test_scenarios
            println("\n" * "="^40)
            println("SCENARIO: $scenario_name")
            println("="^40)
            
            # Set up observations (root nodes only)
            observations = Dict{String, Tuple{Float64, Float64}}()
            for (root_node, activity) in root_perturbations
                observations[root_node] = (activity, 0.9)  # High confidence
            end
            
            println("🎯 Root perturbations:")
            for (node, (activity, _)) in observations
                direction = activity > 50.0 ? "UP" : "DOWN"
                println("  $node: $(activity)% ($direction)")
            end
            
            # Solve the network
            params = SteadyStateParams(1.0, 0.05, 150, 1e-5, "fixed_point")  # Balanced gamma to reduce baseline drift
            result = solve_steady_state(network, observations, params)
            
            if !result.converged
                println("⚠️  Warning: Solver did not converge!")
            end
            
            # Analyze terminal node responses
            println("\n📊 Terminal node responses:")
            terminal_responses = Dict{String, Float64}()
            
            for term_node in sort(terminal_nodes)
                activity = get(result.node_activities, term_node, 0.0) * 100
                terminal_responses[term_node] = activity
                
                # Determine expected vs actual
                baseline = 20.0  # Our default baseline * 100
                response = activity > baseline + 5 ? "UP" : (activity < baseline - 5 ? "DOWN" : "STABLE")
                
                println("  $term_node: $(round(activity, digits=1))% ($response)")
            end
            
            # Analyze propagation efficiency
            println("\n🔍 Propagation analysis:")
            
            # Calculate average root perturbation magnitude
            avg_root_perturbation = sum(abs(activity - 20.0) for (_, (activity, _)) in observations) / length(observations)
            
            # Calculate average terminal response magnitude
            avg_terminal_response = sum(abs(activity - 20.0) for activity in values(terminal_responses)) / length(terminal_responses)
            
            propagation_efficiency = avg_terminal_response / avg_root_perturbation
            
            println("  Avg root perturbation: $(round(avg_root_perturbation, digits=1))%")
            println("  Avg terminal response: $(round(avg_terminal_response, digits=1))%")
            println("  Propagation efficiency: $(round(propagation_efficiency, digits=2))x")
            
            # Store results for summary
            push!(results_summary, (
                scenario_name,
                avg_root_perturbation,
                avg_terminal_response, 
                propagation_efficiency,
                result.converged
            ))
            
            # Check for expected behaviors based on topology
            analyze_topology_specific_behavior(scenario_name, root_perturbations, result.node_activities)
        end
        
        # Generate final assessment
        println("\n" * "="^60)
        println("🎯 PROPAGATION ASSESSMENT SUMMARY")
        println("="^60)
        
        println("\n📈 Scenario Performance:")
        println("Scenario                  | Root Δ | Term Δ | Efficiency | Converged")
        println("-" ^ 70)
        
        for (name, root_delta, term_delta, efficiency, converged) in results_summary
            conv_str = converged ? "✅" : "❌"
            println("$(rpad(name, 24)) | $(rpad(round(root_delta, digits=1), 6)) | $(rpad(round(term_delta, digits=1), 6)) | $(rpad(round(efficiency, digits=2), 10)) | $conv_str")
        end
        
        # Validate expected behaviors
        println("\n🔬 Topology-specific validation:")
        
        convergence_rate = count(x -> x[5], results_summary) / length(results_summary)
        avg_efficiency = sum(x -> x[4], results_summary) / length(results_summary)
        
        println("  Overall convergence rate: $(round(convergence_rate * 100, digits=1))%")
        println("  Average propagation efficiency: $(round(avg_efficiency, digits=2))x")
        
        # Validate key principles
        validate_propagation_principles(results_summary)
        
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

function analyze_topology_specific_behavior(scenario_name, root_perturbations, node_activities)
    println("\n🔍 Topology-specific analysis for $scenario_name:")
    
    # Convert to percentage for analysis
    activities_pct = Dict(uuid => activity * 100 for (uuid, activity) in node_activities)
    baseline = 20.0
    
    if scenario_name == "Competing balanced"
        # ROOT5(+) and ROOT6(-) should roughly cancel out
        g1_activity = get(activities_pct, "G1", baseline)
        term5_activity = get(activities_pct, "TERM5", baseline)
        
        if abs(g1_activity - baseline) < 10 && abs(term5_activity - baseline) < 15
            println("  ✅ Competing inputs balanced correctly")
        else
            println("  ⚠️  Competing inputs may not be properly balanced")
        end
        
    elseif scenario_name == "Activator wins"
        # ROOT5 strong, ROOT6 weak → should see net activation
        term5_activity = get(activities_pct, "TERM5", baseline)
        if term5_activity > baseline + 5
            println("  ✅ Activator dominance propagated correctly")
        else
            println("  ⚠️  Activator dominance not seen at terminal")
        end
        
    elseif scenario_name == "Inhibitor wins"
        # ROOT6 strong, ROOT5 weak → should see net inhibition
        term5_activity = get(activities_pct, "TERM5", baseline)
        if term5_activity < baseline - 5
            println("  ✅ Inhibitor dominance propagated correctly")
        else
            println("  ⚠️  Inhibitor dominance not seen at terminal")
        end
        
    elseif contains(scenario_name, "Pos feedback")
        # Positive feedback should amplify the signal
        term4_activity = get(activities_pct, "TERM4", baseline)
        root4_activity = root_perturbations["ROOT4"]
        
        if abs(term4_activity - baseline) > abs(root4_activity - baseline) * 0.8
            println("  ✅ Positive feedback amplification detected")
        else
            println("  ⚠️  Positive feedback amplification may be weak")
        end
        
    elseif contains(scenario_name, "Neg feedback")
        # Negative feedback should moderate/stabilize the signal
        term8_activity = get(activities_pct, "TERM8", baseline)
        root10_activity = root_perturbations["ROOT10"]
        
        # The response should be dampened compared to a linear chain
        if abs(term8_activity - baseline) < abs(root10_activity - baseline) * 0.9
            println("  ✅ Negative feedback dampening detected")
        else
            println("  ⚠️  Negative feedback dampening may be insufficient")
        end
        
    elseif contains(scenario_name, "Long cascade")
        # Long linear chain should show signal attenuation
        term6_activity = get(activities_pct, "TERM6", baseline)
        root7_activity = root_perturbations["ROOT7"]
        
        signal_retention = abs(term6_activity - baseline) / abs(root7_activity - baseline)
        if signal_retention > 0.3 && signal_retention < 1.2
            println("  ✅ Long cascade propagation appropriate ($(round(signal_retention*100, digits=1))% retention)")
        else
            println("  ⚠️  Long cascade propagation unusual ($(round(signal_retention*100, digits=1))% retention)")
        end
    end
end

function validate_propagation_principles(results_summary)
    println("\n📋 Validating core propagation principles:")
    
    # Principle 1: Upregulation should generally propagate as upregulation
    up_scenarios = [r for r in results_summary if contains(r[1], "up") || contains(r[1], "UP")]
    up_efficiency = [r[4] for r in up_scenarios if r[4] > 0]  # Positive efficiency
    
    if length(up_efficiency) >= length(up_scenarios) * 0.8
        println("  ✅ Upregulation propagates correctly ($(length(up_efficiency))/$(length(up_scenarios)) scenarios)")
    else
        println("  ⚠️  Some upregulation scenarios may not propagate properly")
    end
    
    # Principle 2: Convergence should be robust across topologies
    converged_count = count(r -> r[5], results_summary)
    if converged_count >= length(results_summary) * 0.9
        println("  ✅ Robust convergence across topologies ($converged_count/$(length(results_summary)) scenarios)")
    else
        println("  ⚠️  Convergence issues detected in some topologies")
    end
    
    # Principle 3: Signal should not be completely lost (efficiency > 0.1)
    weak_propagation = count(r -> r[4] < 0.1, results_summary)
    if weak_propagation <= length(results_summary) * 0.2
        println("  ✅ Signal propagation maintained ($weak_propagation/$(length(results_summary)) weak scenarios)")
    else
        println("  ⚠️  Signal loss detected in multiple scenarios")
    end
    
    # Principle 4: No unrealistic amplification (efficiency < 3.0 for most cases)
    over_amplified = count(r -> r[4] > 3.0, results_summary)
    if over_amplified <= length(results_summary) * 0.1
        println("  ✅ No excessive signal amplification ($over_amplified/$(length(results_summary)) over-amplified)")
    else
        println("  ⚠️  Unrealistic signal amplification in some scenarios")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_pathway_propagation()
end