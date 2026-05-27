# Parameter Learning Pipeline for Pathway Dynamics

using Optim
using Statistics
using LinearAlgebra
using Random

"""
Parameters for the learning algorithm.
"""
struct LearningParams
    max_iters::Int              # Maximum optimization iterations
    tolerance::Float64          # Convergence tolerance
    learning_rate::Float64      # Learning rate for gradient-based methods
    regularization::Float64     # L2 regularization weight
    
    # Cross-validation
    cv_folds::Int              # Number of cross-validation folds
    validation_split::Float64   # Fraction of data for validation
    
    # Parameter bounds
    min_hill_h::Float64        # Minimum Hill coefficient
    max_hill_h::Float64        # Maximum Hill coefficient
    min_hill_K::Float64        # Minimum Hill threshold
    max_hill_K::Float64        # Maximum Hill threshold
    min_beta::Float64          # Minimum inhibition strength
    max_beta::Float64          # Maximum inhibition strength
    
    function LearningParams(
        max_iters::Int = 500,
        tolerance::Float64 = 1e-6,
        learning_rate::Float64 = 0.01,
        regularization::Float64 = 0.01,
        cv_folds::Int = 5,
        validation_split::Float64 = 0.2,
        min_hill_h::Float64 = 0.5,
        max_hill_h::Float64 = 5.0,
        min_hill_K::Float64 = 0.05,
        max_hill_K::Float64 = 0.8,
        min_beta::Float64 = 0.1,
        max_beta::Float64 = 10.0
    )
        new(max_iters, tolerance, learning_rate, regularization, cv_folds, validation_split,
            min_hill_h, max_hill_h, min_hill_K, max_hill_K, min_beta, max_beta)
    end
end

"""
Training data point with input observations and expected outputs.
"""
struct TrainingDataPoint
    input_observations::Dict{String, Tuple{Float64, Float64}}  # node_id -> (activity, confidence)
    target_outputs::Dict{String, Float64}                      # node_id -> expected_activity
    weight::Float64                                            # Importance weight for this data point
    
    function TrainingDataPoint(
        inputs::Dict{String, Tuple{Float64, Float64}},
        outputs::Dict{String, Float64},
        weight::Float64 = 1.0
    )
        new(inputs, outputs, weight)
    end
end

"""
Result of parameter learning including learned parameters and performance metrics.
"""
struct LearningResult
    learned_network::ReactionNetwork          # Network with learned parameters
    training_loss::Float64                   # Final training loss
    validation_loss::Float64                 # Validation loss
    parameter_history::Vector{Dict{String, Float64}}  # Parameter evolution
    converged::Bool                          # Whether optimization converged
    iterations::Int                          # Number of iterations used
    solve_time::Float64                      # Time taken for learning
    diagnostics::Dict{String, Any}           # Additional diagnostics
end

"""
Learnable parameter vector representation.
"""
mutable struct LearnableParameters
    # Hill function parameters (per reaction)
    hill_h::Dict{String, Float64}      # target_id -> hill coefficient
    hill_K::Dict{String, Float64}      # target_id -> hill threshold
    
    # Inhibition parameters (per reaction-inhibitor pair)
    inhibitor_betas::Dict{String, Float64}  # "reaction_target:inhibitor_id" -> beta
    inhibitor_ms::Dict{String, Float64}     # "reaction_target:inhibitor_id" -> m
    
    # Sensitivity parameters (per reaction-activator pair)  
    sensitivity_s::Dict{String, Float64}    # "reaction_target:activator_id" -> s
    sensitivity_n::Dict{String, Float64}    # "reaction_target:activator_id" -> n
    
    function LearnableParameters(network::ReactionNetwork)
        reactions = convert_to_reaction_network(network)
        
        hill_h = Dict{String, Float64}()
        hill_K = Dict{String, Float64}()
        inhibitor_betas = Dict{String, Float64}()
        inhibitor_ms = Dict{String, Float64}()
        sensitivity_s = Dict{String, Float64}()
        sensitivity_n = Dict{String, Float64}()
        
        # Initialize from current network parameters
        for reaction in reactions
            target = reaction.target_uuid
            
            # Hill parameters
            hill_h[target] = reaction.params.h
            hill_K[target] = reaction.params.K
            
            # Inhibitor parameters
            for (i, inhibitor) in enumerate(reaction.inhibitor_uuids)
                key = "$target:$inhibitor"
                if i <= length(reaction.params.inhibitor_betas)
                    inhibitor_betas[key] = reaction.params.inhibitor_betas[i]
                    inhibitor_ms[key] = reaction.params.inhibitor_ms[i]
                else
                    inhibitor_betas[key] = 2.0  # Default
                    inhibitor_ms[key] = 2.0     # Default
                end
            end
            
            # Sensitivity parameters
            for (i, activator) in enumerate(reaction.activator_uuids)
                key = "$target:$activator"
                if i <= length(reaction.params.activator_sensitivity_s)
                    sensitivity_s[key] = reaction.params.activator_sensitivity_s[i]
                    sensitivity_n[key] = reaction.params.activator_sensitivity_n[i]
                else
                    sensitivity_s[key] = 0.0  # Default (neutral)
                    sensitivity_n[key] = 2.0  # Default
                end
            end
        end
        
        new(hill_h, hill_K, inhibitor_betas, inhibitor_ms, sensitivity_s, sensitivity_n)
    end
end

"""
Convert learnable parameters back to a reaction network.
"""
function apply_learned_parameters(
    base_network::ReactionNetwork,
    learned_params::LearnableParameters
)::ReactionNetwork
    
    # Create new network with updated parameters
    new_network = deepcopy(base_network)
    reactions = convert_to_reaction_network(new_network)
    
    # Create new reactions with updated parameters
    updated_reactions = Reaction[]
    
    for reaction in reactions
        target = reaction.target_uuid
        
        # Get updated Hill parameters
        new_h = haskey(learned_params.hill_h, target) ? learned_params.hill_h[target] : reaction.params.h
        new_K = haskey(learned_params.hill_K, target) ? learned_params.hill_K[target] : reaction.params.K
        
        # Get updated inhibitor parameters
        new_inhibitor_betas = copy(reaction.params.inhibitor_betas)
        new_inhibitor_ms = copy(reaction.params.inhibitor_ms)
        
        for (i, inhibitor) in enumerate(reaction.inhibitor_uuids)
            key = "$target:$inhibitor"
            if haskey(learned_params.inhibitor_betas, key) && i <= length(new_inhibitor_betas)
                new_inhibitor_betas[i] = learned_params.inhibitor_betas[key]
            end
            if haskey(learned_params.inhibitor_ms, key) && i <= length(new_inhibitor_ms)
                new_inhibitor_ms[i] = learned_params.inhibitor_ms[key]
            end
        end
        
        # Get updated sensitivity parameters  
        new_sensitivity_s = copy(reaction.params.activator_sensitivity_s)
        new_sensitivity_n = copy(reaction.params.activator_sensitivity_n)
        
        for (i, activator) in enumerate(reaction.activator_uuids)
            key = "$target:$activator"
            if haskey(learned_params.sensitivity_s, key) && i <= length(new_sensitivity_s)
                new_sensitivity_s[i] = learned_params.sensitivity_s[key]
            end
            if haskey(learned_params.sensitivity_n, key) && i <= length(new_sensitivity_n)
                new_sensitivity_n[i] = learned_params.sensitivity_n[key]
            end
        end
        
        # Create new ReactionParams with updated values
        new_params = ReactionParams(
            new_h, new_K,
            reaction.params.activator_weights, new_sensitivity_s, new_sensitivity_n, reaction.params.activator_sensitivity_K,
            new_inhibitor_betas, new_inhibitor_ms,
            reaction.params.substrate_weights,
            reaction.params.consumption_lambdas, reaction.params.production_etas,
            reaction.params.replenishment_rho, reaction.params.decay_delta
        )
        
        # Create new reaction with updated parameters
        new_reaction = Reaction(
            reaction.target_uuid,
            reaction.activator_uuids,
            reaction.inhibitor_uuids,
            reaction.substrate_uuids,
            reaction.product_uuids,
            new_params,
            reaction.is_and_gate
        )
        
        push!(updated_reactions, new_reaction)
    end
    
    return new_network
end

"""
Convert learnable parameters to a flat vector for optimization.
"""
function parameters_to_vector(params::LearnableParameters)::Vector{Float64}
    vec = Float64[]
    
    # Hill parameters
    for (target, h) in params.hill_h
        push!(vec, h)
    end
    for (target, K) in params.hill_K
        push!(vec, K)
    end
    
    # Inhibitor parameters
    for (key, beta) in params.inhibitor_betas
        push!(vec, beta)
    end
    for (key, m) in params.inhibitor_ms
        push!(vec, m)
    end
    
    # Sensitivity parameters
    for (key, s) in params.sensitivity_s
        push!(vec, s)
    end
    for (key, n) in params.sensitivity_n
        push!(vec, n)
    end
    
    return vec
end

"""
Update learnable parameters from a flat vector.
"""
function vector_to_parameters!(params::LearnableParameters, vec::Vector{Float64})
    idx = 1
    
    # Hill parameters
    for target in keys(params.hill_h)
        params.hill_h[target] = vec[idx]
        idx += 1
    end
    for target in keys(params.hill_K)
        params.hill_K[target] = vec[idx]
        idx += 1
    end
    
    # Inhibitor parameters
    for key in keys(params.inhibitor_betas)
        params.inhibitor_betas[key] = vec[idx]
        idx += 1
    end
    for key in keys(params.inhibitor_ms)
        params.inhibitor_ms[key] = vec[idx]
        idx += 1
    end
    
    # Sensitivity parameters
    for key in keys(params.sensitivity_s)
        params.sensitivity_s[key] = vec[idx]
        idx += 1
    end
    for key in keys(params.sensitivity_n)
        params.sensitivity_n[key] = vec[idx]
        idx += 1
    end
end

"""
Main parameter learning function using gradient-based optimization.
"""
function learn_parameters(
    network::ReactionNetwork,
    training_data::Vector{TrainingDataPoint},
    learning_params::LearningParams = LearningParams()
)::LearningResult
    
    start_time = time()
    
    println("🎓 Starting parameter learning")
    println("📊 Training data: $(length(training_data)) data points")
    println("⚙️  Parameters: $(learning_params.max_iters) max iters, $(learning_params.tolerance) tolerance")
    
    # Initialize learnable parameters
    learnable_params = LearnableParameters(network)
    n_params = length(parameters_to_vector(learnable_params))
    println("🔧 Learning $(n_params) parameters")
    
    # Split training/validation data
    n_train = length(training_data)
    n_val = Int(round(n_train * learning_params.validation_split))
    n_train_actual = n_train - n_val
    
    # Shuffle and split
    shuffled_data = shuffle(training_data)
    train_data = shuffled_data[1:n_train_actual]
    val_data = shuffled_data[(n_train_actual+1):end]
    
    println("📈 Split: $(n_train_actual) training, $(n_val) validation")
    
    # Parameter bounds
    bounds_lower = Float64[]
    bounds_upper = Float64[]
    
    # Hill h bounds
    for _ in keys(learnable_params.hill_h)
        push!(bounds_lower, learning_params.min_hill_h)
        push!(bounds_upper, learning_params.max_hill_h)
    end
    # Hill K bounds  
    for _ in keys(learnable_params.hill_K)
        push!(bounds_lower, learning_params.min_hill_K)
        push!(bounds_upper, learning_params.max_hill_K)
    end
    # Inhibitor beta bounds
    for _ in keys(learnable_params.inhibitor_betas)
        push!(bounds_lower, learning_params.min_beta)
        push!(bounds_upper, learning_params.max_beta)
    end
    # Inhibitor m bounds
    for _ in keys(learnable_params.inhibitor_ms)
        push!(bounds_lower, 0.5)
        push!(bounds_upper, 5.0)
    end
    # Sensitivity s bounds
    for _ in keys(learnable_params.sensitivity_s)
        push!(bounds_lower, -2.0)
        push!(bounds_upper, 2.0)
    end
    # Sensitivity n bounds
    for _ in keys(learnable_params.sensitivity_n)
        push!(bounds_lower, 0.5)
        push!(bounds_upper, 5.0)
    end
    
    # Track parameter history
    parameter_history = Dict{String, Float64}[]
    
    # Objective function
    function objective(param_vec::Vector{Float64})
        # Update parameters
        vector_to_parameters!(learnable_params, param_vec)
        
        # Use a simplified forward model with current parameters
        reactions = convert_to_reaction_network(network)
        
        # Compute training loss
        total_loss = 0.0
        
        for data_point in train_data
            # Use simplified forward model with learned parameters
            predictions = predict_with_learned_params(
                reactions, data_point.input_observations, learnable_params
            )
            
            if isnothing(predictions)
                return 1e6  # Large penalty for failed prediction
            end
            
            # Compute prediction error
            point_loss = 0.0
            n_targets = 0
            
            for (target_node, expected_activity) in data_point.target_outputs
                if haskey(predictions, target_node)
                    predicted_activity = predictions[target_node] * 100.0  # Convert to percentage
                    error = (predicted_activity - expected_activity)^2
                    point_loss += error
                    n_targets += 1
                end
            end
            
            if n_targets > 0
                point_loss = point_loss / n_targets  # Average error per target
                total_loss += data_point.weight * point_loss
            end
        end
        
        # Normalize by number of data points
        if length(train_data) > 0
            total_loss = total_loss / length(train_data)
        end
        
        # Add L2 regularization
        regularization_loss = 0.0
        for param in param_vec
            regularization_loss += param^2
        end
        total_loss += learning_params.regularization * regularization_loss
        
        # Save parameter snapshot periodically
        if length(parameter_history) < 100 && length(parameter_history) % 10 == 0
            param_dict = Dict{String, Float64}()
            param_dict["iteration"] = Float64(length(parameter_history))
            param_dict["loss"] = total_loss
            push!(parameter_history, param_dict)
        end
        
        return total_loss
    end
    
    # Run optimization
    println("🚀 Starting optimization...")
    initial_params = parameters_to_vector(learnable_params)
    
    # Apply bounds by clamping
    initial_params = max.(min.(initial_params, bounds_upper), bounds_lower)
    
    # Use Nelder-Mead for robustness (derivative-free)
    optimization_result = optimize(
        objective,
        initial_params,
        NelderMead(),
        Optim.Options(
            iterations = learning_params.max_iters,
            g_tol = learning_params.tolerance,
            show_trace = true,
            show_every = 50
        )
    )
    
    # Extract final parameters
    final_params = Optim.minimizer(optimization_result)
    vector_to_parameters!(learnable_params, final_params)
    final_network = apply_learned_parameters(network, learnable_params)
    
    # Compute final losses
    training_loss = objective(final_params)
    
    # Validation loss
    validation_loss = 0.0
    if !isempty(val_data)
        for data_point in val_data
            ss_params = SteadyStateParams(1.0, 0.05, 150, 1e-5, "fixed_point")
            result = solve_steady_state(final_network, data_point.input_observations, ss_params)
            
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
                    validation_loss += data_point.weight * (point_loss / n_targets)
                end
            end
        end
        
        if length(val_data) > 0
            validation_loss = validation_loss / length(val_data)
        end
    end
    
    solve_time = time() - start_time
    converged = Optim.converged(optimization_result)
    iterations = Optim.iterations(optimization_result)
    
    println("✅ Parameter learning completed!")
    println("🎯 Training loss: $(round(training_loss, digits=4))")
    println("📊 Validation loss: $(round(validation_loss, digits=4))")
    println("🔄 Iterations: $iterations (converged: $converged)")
    println("⏱️  Time: $(round(solve_time, digits=2))s")
    
    return LearningResult(
        final_network,
        training_loss,
        validation_loss,
        parameter_history,
        converged,
        iterations,
        solve_time,
        Dict(
            "method" => "nelder_mead",
            "n_parameters" => n_params,
            "n_training_points" => n_train_actual,
            "n_validation_points" => n_val,
            "optimization_result" => optimization_result
        )
    )
end

"""
Generate synthetic training data for testing the learning pipeline.
"""
function generate_synthetic_training_data(
    network::ReactionNetwork,
    n_data_points::Int = 20
)::Vector{TrainingDataPoint}
    
    training_data = TrainingDataPoint[]
    
    # Get root nodes (sources of perturbations)
    all_children = Set(edge.child_uuid for edge in network.edges)
    all_parents = Set(edge.parent_uuid for edge in network.edges)
    root_nodes = [uuid for uuid in all_parents if !(uuid in all_children)]
    
    # Get terminal nodes (outputs we want to predict)
    terminal_nodes = [uuid for uuid in all_children if !(uuid in all_parents)]
    
    println("🧪 Generating synthetic data from $(length(root_nodes)) roots → $(length(terminal_nodes)) terminals")
    
    for i in 1:n_data_points
        # Random perturbation of root nodes
        observations = Dict{String, Tuple{Float64, Float64}}()
        
        # Randomly perturb 1-3 root nodes
        n_perturb = rand(1:min(3, length(root_nodes)))
        perturbed_roots = shuffle(root_nodes)[1:n_perturb]
        
        for root in perturbed_roots
            activity = rand(20.0:80.0)  # Random activity level
            confidence = 0.9
            observations[root] = (activity, confidence)
        end
        
        # Simulate with current network to get "ground truth"
        ss_params = SteadyStateParams(1.0, 0.05, 150, 1e-5, "fixed_point")
        result = solve_steady_state(network, observations, ss_params)
        
        if result.converged
            # Create target outputs (terminal node activities)
            target_outputs = Dict{String, Float64}()
            for terminal in terminal_nodes
                if haskey(result.node_activities, terminal)
                    activity = result.node_activities[terminal] * 100.0
                    target_outputs[terminal] = activity
                end
            end
            
            # Add some intermediate nodes as targets too
            intermediate_nodes = setdiff(collect(keys(network.nodes)), vcat(root_nodes, terminal_nodes))
            if !isempty(intermediate_nodes)
                # Sample a few intermediate nodes
                n_intermediate = min(2, length(intermediate_nodes))
                sampled_intermediate = shuffle(intermediate_nodes)[1:n_intermediate]
                
                for node in sampled_intermediate
                    if haskey(result.node_activities, node)
                        activity = result.node_activities[node] * 100.0
                        target_outputs[node] = activity
                    end
                end
            end
            
            if !isempty(target_outputs)
                push!(training_data, TrainingDataPoint(observations, target_outputs, 1.0))
            end
        end
    end
    
    println("✅ Generated $(length(training_data)) training data points")
    return training_data
end

"""
Simplified prediction function using learned parameters.
"""
function predict_with_learned_params(
    reactions::Vector{Reaction},
    observations::Dict{String, Tuple{Float64, Float64}},
    learned_params::LearnableParameters
)::Union{Dict{String, Float64}, Nothing}
    
    try
        # Initialize activities
        activities = Dict{String, Float64}()
        all_nodes = Set{String}()
        
        for reaction in reactions
            push!(all_nodes, reaction.target_uuid)
            for uuid in reaction.activator_uuids
                push!(all_nodes, uuid)
            end
            for uuid in reaction.inhibitor_uuids
                push!(all_nodes, uuid)
            end
        end
        
        # Initialize to baseline
        for node in all_nodes
            activities[node] = 0.2
        end
        
        # Apply observations
        for (node_id, (activity, confidence)) in observations
            if haskey(activities, node_id)
                activities[node_id] = activity / 100.0
            end
        end
        
        # Simple fixed-point iteration with learned parameters
        for iter in 1:50  # Limited iterations for efficiency
            max_change = 0.0
            new_activities = copy(activities)
            
            for reaction in reactions
                target = reaction.target_uuid
                
                # Get activities
                activator_activities = [get(activities, uuid, 0.0) for uuid in reaction.activator_uuids]
                inhibitor_activities = [get(activities, uuid, 0.0) for uuid in reaction.inhibitor_uuids]
                
                if isempty(activator_activities)
                    continue
                end
                
                # Use learned Hill parameters
                h = get(learned_params.hill_h, target, 1.5)
                K = get(learned_params.hill_K, target, 0.3)
                
                # Simplified Hill activation (geometric mean)
                if !isempty(activator_activities)
                    A_geom = exp(sum(log.(max.(activator_activities, 1e-12))) / length(activator_activities))
                else
                    A_geom = 0.0
                end
                
                # Hill activation
                activation = (A_geom^h) / (A_geom^h + K^h)
                
                # Simple inhibition (multiplicative)
                inhibition = 1.0
                if !isempty(inhibitor_activities)
                    # Use learned inhibition strength
                    beta = 2.0  # Default strength
                    if !isempty(reaction.inhibitor_uuids)
                        key = "$(target):$(reaction.inhibitor_uuids[1])"
                        beta = get(learned_params.inhibitor_betas, key, 2.0)
                    end
                    
                    max_inhibitor = maximum(inhibitor_activities)
                    inhibition = 1.0 / (1.0 + beta * max_inhibitor)
                end
                
                # Final activity
                new_activity = activation * inhibition
                new_activities[target] = clamp(new_activity, 0.0, 1.0)
                
                # Track convergence
                change = abs(new_activities[target] - activities[target])
                max_change = max(max_change, change)
            end
            
            activities = new_activities
            
            if max_change < 1e-4
                break
            end
        end
        
        return activities
        
    catch e
        return nothing
    end
end

"""
Default learning parameters.
"""
function default_learning_params()::LearningParams
    return LearningParams()
end