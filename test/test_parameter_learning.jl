#!/usr/bin/env julia

# Test parameter learning pipeline
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Random
using Statistics
using CSV
using DataFrames

# Add the src directory to the load path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
include("../src/io/tsv_parser.jl")
include("../src/core/sensitivity.jl")
include("../src/core/aggregators.jl") 
include("../src/core/hill_functions.jl")
include("../src/core/reaction_model.jl")
include("../src/solvers/steady_state.jl")
include("../src/learning/parameter_learning.jl")

function test_parameter_learning_pipeline()
    println("🎓 Testing Parameter Learning Pipeline")
    println("=" ^ 50)
    
    try
        # Load a network for testing
        logic_path = joinpath(@__DIR__, "..", "examples", "propagation_test_network.tsv")
        uuid_path = joinpath(@__DIR__, "..", "examples", "propagation_test_uuid_mapping.tsv")
        
        original_network = parse_complete_network(logic_path, uuid_path, nothing)
        
        println("📁 Network loaded with $(length(original_network.nodes)) nodes, $(length(original_network.edges)) edges")
        
        # Test 1: Synthetic data generation
        println("\\n" * "="^50)
        println("TEST 1: Synthetic Data Generation")
        println("="^50)
        
        Random.seed!(42)  # For reproducible results
        
        training_data = generate_synthetic_training_data(original_network, 15)
        
        println("📊 Generated training data statistics:")
        println("  Data points: $(length(training_data))")
        
        # Analyze data diversity
        all_inputs = Set{String}()
        all_outputs = Set{String}()
        
        for data_point in training_data
            for (input_node, _) in data_point.input_observations
                push!(all_inputs, input_node)
            end
            for (output_node, _) in data_point.target_outputs
                push!(all_outputs, output_node)
            end
        end
        
        println("  Input nodes covered: $(length(all_inputs))")
        println("  Output nodes covered: $(length(all_outputs))")
        
        # Show sample data point
        if !isempty(training_data)
            sample_point = training_data[1]
            println("\\n📝 Sample training data point:")
            println("  Inputs:")
            for (node, (activity, confidence)) in sample_point.input_observations
                println("    $node: $(round(activity, digits=1))% (confidence: $confidence)")
            end
            println("  Expected outputs:")
            for (node, activity) in sample_point.target_outputs
                println("    $node: $(round(activity, digits=1))%")
            end
        end
        
        # Test 2: Parameter extraction and conversion
        println("\\n" * "="^50) 
        println("TEST 2: Parameter Representation")
        println("="^50)
        
        learnable_params = LearnableParameters(original_network)
        param_vector = parameters_to_vector(learnable_params)
        
        println("🔧 Learnable parameters:")
        println("  Hill h parameters: $(length(learnable_params.hill_h))")
        println("  Hill K parameters: $(length(learnable_params.hill_K))")
        println("  Inhibitor betas: $(length(learnable_params.inhibitor_betas))")
        println("  Inhibitor ms: $(length(learnable_params.inhibitor_ms))")
        println("  Sensitivity s: $(length(learnable_params.sensitivity_s))")
        println("  Sensitivity n: $(length(learnable_params.sensitivity_n))")
        println("  Total parameters: $(length(param_vector))")
        
        # Test round-trip conversion
        learnable_params_copy = LearnableParameters(original_network)
        vector_to_parameters!(learnable_params_copy, param_vector)
        param_vector_2 = parameters_to_vector(learnable_params_copy)
        
        conversion_error = norm(param_vector - param_vector_2)
        println("  Round-trip conversion error: $(round(conversion_error, digits=8))")
        
        if conversion_error < 1e-10
            println("  ✅ Parameter conversion working correctly")
        else
            println("  ⚠️  Parameter conversion has errors")
        end
        
        # Test 3: Baseline performance (before learning)
        println("\\n" * "="^50)
        println("TEST 3: Baseline Performance")
        println("="^50)
        
        # Compute baseline loss with original parameters
        baseline_loss = 0.0
        n_evaluated = 0
        
        for data_point in training_data[1:min(5, length(training_data))]  # Test on first 5 points
            ss_params = SteadyStateParams(1.0, 0.05, 150, 1e-5, "fixed_point")
            result = solve_steady_state(original_network, data_point.input_observations, ss_params)
            
            if result.converged
                point_loss = 0.0
                n_targets = 0
                
                for (target_node, expected_activity) in data_point.target_outputs
                    if haskey(result.node_activities, target_node)
                        predicted_activity = result.node_activities[target_node] * 100.0
                        error = (predicted_activity - expected_activity)^2
                        point_loss += error
                        n_targets += 1
                    end
                end
                
                if n_targets > 0
                    baseline_loss += point_loss / n_targets
                    n_evaluated += 1
                end
            end
        end
        
        if n_evaluated > 0
            baseline_loss = baseline_loss / n_evaluated
            println("📊 Baseline performance (MSE): $(round(baseline_loss, digits=4))")
        end
        
        # Test 4: Small-scale learning test
        println("\\n" * "="^50)
        println("TEST 4: Small-Scale Learning")
        println("="^50)
        
        # Use a smaller subset for faster testing
        small_training_data = training_data[1:min(8, length(training_data))]
        
        # Create learning parameters for quick test
        learning_params = LearningParams(
            50,     # max_iters - fewer iterations for testing
            1e-4,   # tolerance - looser tolerance
            0.01,   # learning_rate
            0.1,    # regularization - higher regularization for stability
            5,      # cv_folds
            0.2     # validation_split
        )
        
        println("🚀 Running parameter learning ($(length(small_training_data)) data points, $(learning_params.max_iters) max iters)...")
        
        learning_result = learn_parameters(original_network, small_training_data, learning_params)
        
        # Analyze results
        println("\\n🎯 Learning Results:")
        println("  Converged: $(learning_result.converged)")
        println("  Iterations: $(learning_result.iterations)")
        println("  Training loss: $(round(learning_result.training_loss, digits=4))")
        println("  Validation loss: $(round(learning_result.validation_loss, digits=4))")
        println("  Solve time: $(round(learning_result.solve_time, digits=2))s")
        
        # Compare with baseline
        if n_evaluated > 0 && baseline_loss > 0
            improvement = ((baseline_loss - learning_result.training_loss) / baseline_loss) * 100
            println("  Improvement over baseline: $(round(improvement, digits=1))%")
            
            if improvement > 5.0
                println("  ✅ Significant improvement detected!")
            elseif improvement > 0.0
                println("  ✅ Some improvement detected")
            else
                println("  ⚠️  No improvement (may need more iterations/data)")
            end
        end
        
        # Test 5: Parameter validation
        println("\\n" * "="^50)
        println("TEST 5: Learned Parameter Validation")
        println("="^50)
        
        # Check if learned parameters are reasonable
        learned_network = learning_result.learned_network
        learned_reactions = convert_to_reaction_network(learned_network)
        
        println("🔍 Learned parameter analysis:")
        
        # Hill coefficients
        hill_hs = [reaction.params.h for reaction in learned_reactions]
        println("  Hill h: mean=$(round(mean(hill_hs), digits=2)), range=$(round(minimum(hill_hs), digits=2))-$(round(maximum(hill_hs), digits=2))")
        
        # Hill thresholds
        hill_Ks = [reaction.params.K for reaction in learned_reactions]  
        println("  Hill K: mean=$(round(mean(hill_Ks), digits=2)), range=$(round(minimum(hill_Ks), digits=2))-$(round(maximum(hill_Ks), digits=2))")
        
        # Inhibition strengths
        all_betas = Float64[]
        for reaction in learned_reactions
            if !isempty(reaction.params.inhibitor_betas)
                append!(all_betas, reaction.params.inhibitor_betas)
            end
        end
        
        if !isempty(all_betas)
            println("  Inhibitor β: mean=$(round(mean(all_betas), digits=2)), range=$(round(minimum(all_betas), digits=2))-$(round(maximum(all_betas), digits=2))")
        end
        
        # Check for reasonable bounds
        params_in_bounds = true
        if any(h -> h < 0.1 || h > 10.0, hill_hs)
            println("  ⚠️  Some Hill h parameters outside reasonable bounds")
            params_in_bounds = false
        end
        if any(K -> K < 0.01 || K > 1.0, hill_Ks)
            println("  ⚠️  Some Hill K parameters outside reasonable bounds")
            params_in_bounds = false
        end
        
        if params_in_bounds
            println("  ✅ All parameters within reasonable bounds")
        end
        
        # Test 6: Prediction accuracy
        println("\\n" * "="^50)
        println("TEST 6: Prediction Accuracy Test")
        println("="^50)
        
        # Test learned network on a fresh data point
        if length(training_data) > 8
            test_point = training_data[end]  # Use last point as test
            
            ss_params = SteadyStateParams(1.0, 0.05, 150, 1e-5, "fixed_point")
            
            # Original network prediction
            original_result = solve_steady_state(original_network, test_point.input_observations, ss_params)
            
            # Learned network prediction  
            learned_result = solve_steady_state(learned_network, test_point.input_observations, ss_params)
            
            if original_result.converged && learned_result.converged
                println("🔬 Prediction comparison:")
                println("Node\\tTarget\\tOriginal\\tLearned\\tError")
                
                for (target_node, expected_activity) in test_point.target_outputs
                    if haskey(original_result.node_activities, target_node) && haskey(learned_result.node_activities, target_node)
                        orig_pred = original_result.node_activities[target_node] * 100.0
                        learned_pred = learned_result.node_activities[target_node] * 100.0
                        error = abs(learned_pred - expected_activity)
                        
                        println("$(rpad(target_node, 8))\\t$(round(expected_activity, digits=1))%\\t$(round(orig_pred, digits=1))%\\t\\t$(round(learned_pred, digits=1))%\\t$(round(error, digits=1))%")
                    end
                end
            end
        end
        
        println("\\n" * "="^60)
        println("🎉 PARAMETER LEARNING PIPELINE TEST COMPLETE")
        println("="^60)
        
        # Summary
        success_indicators = 0
        total_indicators = 5
        
        if length(training_data) >= 8
            success_indicators += 1
            println("✅ Synthetic data generation: PASS")
        else
            println("⚠️  Synthetic data generation: LIMITED")
        end
        
        if conversion_error < 1e-10
            success_indicators += 1
            println("✅ Parameter representation: PASS")
        else
            println("⚠️  Parameter representation: ISSUES")
        end
        
        if learning_result.converged || learning_result.iterations >= 25
            success_indicators += 1
            println("✅ Optimization execution: PASS")
        else
            println("⚠️  Optimization execution: EARLY_TERMINATION")
        end
        
        if learning_result.training_loss < 1000.0  # Reasonable loss value
            success_indicators += 1
            println("✅ Training loss: REASONABLE")
        else
            println("⚠️  Training loss: HIGH")
        end
        
        if params_in_bounds
            success_indicators += 1
            println("✅ Parameter bounds: PASS")
        else
            println("⚠️  Parameter bounds: SOME_OUT_OF_BOUNDS")
        end
        
        success_rate = success_indicators / total_indicators * 100
        println("\\n📊 Overall success rate: $(round(success_rate, digits=1))% ($(success_indicators)/$(total_indicators))")
        
        if success_rate >= 80.0
            println("🎉 Parameter learning pipeline is WORKING WELL!")
            return true
        elseif success_rate >= 60.0
            println("⚠️  Parameter learning pipeline is PARTIALLY WORKING")
            return true
        else
            println("❌ Parameter learning pipeline needs SIGNIFICANT IMPROVEMENT")
            return false
        end
        
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
    test_parameter_learning_pipeline()
end