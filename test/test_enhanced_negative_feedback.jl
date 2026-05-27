#!/usr/bin/env julia

# Test enhanced negative feedback loop handling
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV
using DataFrames
using JSON3
using LinearAlgebra

# Add the src directory to the load path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
include("../src/io/tsv_parser.jl")
include("../src/core/sensitivity.jl")
include("../src/core/aggregators.jl") 
include("../src/core/hill_functions.jl")
include("../src/core/reaction_model.jl")
include("../src/core/feedback_enhancements.jl")
include("../src/solvers/steady_state.jl")
include("../src/solvers/enhanced_steady_state.jl")

function test_enhanced_negative_feedback()
    println("🔄 Testing Enhanced Negative Feedback Handling")
    println("=" ^ 60)
    
    try
        # Load the propagation test network
        logic_path = joinpath(@__DIR__, "..", "examples", "propagation_test_network.tsv")
        uuid_path = joinpath(@__DIR__, "..", "examples", "propagation_test_uuid_mapping.tsv")
        
        network = parse_complete_network(logic_path, uuid_path, nothing)
        
        println("📁 Network loaded with $(length(network.nodes)) nodes, $(length(network.edges)) edges")
        
        # Test scenarios focusing on negative feedback
        test_scenarios = [
            ("Standard Solver - Neg FB Up", "standard", Dict("ROOT10" => 80.0)),
            ("Enhanced Solver - Neg FB Up", "enhanced", Dict("ROOT10" => 80.0)),
            ("Standard Solver - Neg FB Down", "standard", Dict("ROOT10" => 20.0)),
            ("Enhanced Solver - Neg FB Down", "enhanced", Dict("ROOT10" => 20.0)),
            ("Standard Solver - Strong Inhibitor", "standard", Dict("ROOT5" => 30.0, "ROOT6" => 90.0)),
            ("Enhanced Solver - Strong Inhibitor", "enhanced", Dict("ROOT5" => 30.0, "ROOT6" => 90.0))
        ]
        
        results_comparison = []
        
        for (scenario_name, solver_type, perturbations) in test_scenarios
            println("\n" * "="^50)
            println("SCENARIO: $scenario_name")
            println("="^50)
            
            # Set up observations
            observations = Dict{String, Tuple{Float64, Float64}}()
            for (root, activity) in perturbations
                observations[root] = (activity, 0.9)
            end
            
            println("🎯 Perturbations:")
            for (node, (activity, _)) in observations
                println("  $node: $(activity)%")
            end
            
            # Choose solver
            result = if solver_type == "enhanced"
                params = SteadyStateParams(1.0, 0.1, 200, 1e-5, "feedback_aware")
                solve_steady_state_enhanced(network, observations, params)
            else
                params = SteadyStateParams(1.0, 0.05, 150, 1e-5, "fixed_point")
                solve_steady_state(network, observations, params)
            end
            
            println("\n🧮 Solver performance:")
            println("  Converged: $(result.converged)")
            println("  Iterations: $(result.iterations)")
            println("  Solve time: $(round(result.solve_time, digits=3))s")
            if haskey(result.diagnostics, "final_damping")
                println("  Final damping: $(round(result.diagnostics["final_damping"], digits=3))")
            end
            
            # Analyze negative feedback effectiveness
            println("\n🔍 Negative feedback analysis:")
            
            if "ROOT10" in keys(perturbations)
                # ROOT10 → Q1 → R1 → S1 → [TERM8, Q1(-)]
                root10_input = perturbations["ROOT10"]
                q1_activity = get(result.node_activities, "Q1", 0.2) * 100
                r1_activity = get(result.node_activities, "R1", 0.2) * 100
                s1_activity = get(result.node_activities, "S1", 0.2) * 100
                term8_activity = get(result.node_activities, "TERM8", 0.2) * 100
                
                println("  ROOT10 input: $(root10_input)%")
                println("  Q1 (target of feedback): $(round(q1_activity, digits=1))%")
                println("  R1 (intermediate): $(round(r1_activity, digits=1))%")
                println("  S1 (feedback source): $(round(s1_activity, digits=1))%")
                println("  TERM8 (terminal output): $(round(term8_activity, digits=1))%")
                
                # Calculate dampening effectiveness
                input_deviation = abs(root10_input - 20.0)  # 20% is baseline
                output_deviation = abs(term8_activity - 20.0)
                dampening_factor = output_deviation / input_deviation
                
                println("  Input deviation from baseline: $(round(input_deviation, digits=1))%")
                println("  Output deviation from baseline: $(round(output_deviation, digits=1))%")
                println("  Dampening factor: $(round(dampening_factor, digits=2))x")
                
                if dampening_factor < 0.7
                    println("  ✅ Good negative feedback dampening")
                elseif dampening_factor < 0.9
                    println("  ⚠️  Moderate negative feedback dampening")
                else
                    println("  ❌ Poor negative feedback dampening")
                end
                
                push!(results_comparison, (scenario_name, solver_type, input_deviation, output_deviation, dampening_factor, result.converged))
                
            elseif "ROOT6" in keys(perturbations)
                # ROOT5(+), ROOT6(-) → G1 → H1 → TERM5
                root5_input = get(perturbations, "ROOT5", 20.0)
                root6_input = get(perturbations, "ROOT6", 20.0)
                g1_activity = get(result.node_activities, "G1", 0.2) * 100
                h1_activity = get(result.node_activities, "H1", 0.2) * 100
                term5_activity = get(result.node_activities, "TERM5", 0.2) * 100
                
                println("  ROOT5 (activator): $(root5_input)%")
                println("  ROOT6 (inhibitor): $(root6_input)%")
                println("  G1 (competition site): $(round(g1_activity, digits=1))%")
                println("  H1 (downstream): $(round(h1_activity, digits=1))%")
                println("  TERM5 (terminal): $(round(term5_activity, digits=1))%")
                
                # For strong inhibitor case, terminal should be low
                if root6_input > root5_input + 20  # Strong inhibitor scenario
                    if term5_activity < 35.0
                        println("  ✅ Strong inhibition effective")
                    else
                        println("  ⚠️  Strong inhibition insufficient")
                    end
                end
                
                net_input = root5_input - root6_input  # Net effect
                output_deviation = abs(term5_activity - 20.0)
                push!(results_comparison, (scenario_name, solver_type, abs(net_input), output_deviation, output_deviation / max(abs(net_input), 1), result.converged))
            end
        end
        
        # Compare results
        println("\n" * "="^60)
        println("🎯 ENHANCED vs STANDARD SOLVER COMPARISON")
        println("="^60)
        
        println("\nScenario                              | Solver   | Input Δ | Output Δ | Damping | Conv")
        println(repeat("-", 85))
        
        for (name, solver, input_dev, output_dev, damping, converged) in results_comparison
            conv_str = converged ? "✅" : "❌"
            name_truncated = length(name) > 35 ? name[1:35] : name
            println("$(rpad(name_truncated, 36)) | $(rpad(solver, 8)) | $(rpad(round(input_dev, digits=1), 7)) | $(rpad(round(output_dev, digits=1), 8)) | $(rpad(round(damping, digits=2), 7)) | $conv_str")
        end
        
        # Analyze improvements
        println("\n🔬 Performance Analysis:")
        
        standard_results = [r for r in results_comparison if r[2] == "standard"]
        enhanced_results = [r for r in results_comparison if r[2] == "enhanced"]
        
        if length(standard_results) > 0 && length(enhanced_results) > 0
            avg_standard_damping = sum(r[5] for r in standard_results) / length(standard_results)
            avg_enhanced_damping = sum(r[5] for r in enhanced_results) / length(enhanced_results)
            
            improvement = ((avg_standard_damping - avg_enhanced_damping) / avg_standard_damping) * 100
            
            println("  Average dampening factor:")
            println("    Standard solver: $(round(avg_standard_damping, digits=2))x")
            println("    Enhanced solver: $(round(avg_enhanced_damping, digits=2))x")
            println("    Improvement: $(round(improvement, digits=1))%")
            
            if improvement > 10
                println("  ✅ Significant improvement in negative feedback effectiveness")
            elseif improvement > 0
                println("  ⚠️  Modest improvement in negative feedback effectiveness")
            else
                println("  ❌ No improvement in negative feedback effectiveness")
            end
        end
        
        # Test individual enhancements
        println("\n🧪 Testing individual enhancement components:")
        
        # Test feedback loop detection
        feedback_loops = detect_feedback_loops(network)
        println("  Feedback loop detection:")
        for (loop_type, loops) in feedback_loops
            println("    $loop_type feedback: $(length(loops)) loops detected")
        end
        
        # Test parameter specialization
        test_reactions = convert_to_reaction_network_with_feedback_enhancement(network)
        standard_reactions = convert_to_reaction_network(network)
        
        println("  Parameter specialization:")
        println("    Standard reactions: $(length(standard_reactions))")
        println("    Enhanced reactions: $(length(test_reactions))")
        
        # Compare inhibition strengths
        standard_inhibition_strengths = Float64[]
        enhanced_inhibition_strengths = Float64[]
        
        for reaction in standard_reactions
            if length(reaction.params.inhibitor_betas) > 0
                append!(standard_inhibition_strengths, reaction.params.inhibitor_betas)
            end
        end
        
        for reaction in test_reactions
            if length(reaction.params.inhibitor_betas) > 0
                append!(enhanced_inhibition_strengths, reaction.params.inhibitor_betas)
            end
        end
        
        if !isempty(standard_inhibition_strengths) && !isempty(enhanced_inhibition_strengths)
            println("    Average inhibition strength:")
            println("      Standard: $(round(mean(standard_inhibition_strengths), digits=1))")
            println("      Enhanced: $(round(mean(enhanced_inhibition_strengths), digits=1))")
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
    test_enhanced_negative_feedback()
end