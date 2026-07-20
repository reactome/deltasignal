# Enhanced Steady-State Solver with Feedback-Aware Dynamics

using Statistics
using LinearAlgebra

"""
Enhanced steady-state solver that applies different strategies for different network motifs.
"""
function solve_steady_state_enhanced(
    network::ReactionNetwork,
    observations::Dict{String, Tuple{Float64, Float64}},
    params::SteadyStateParams = default_steady_state_params()
)::SolverResult
    
    start_time = time()
    
    # Convert network to reactions with feedback enhancement
    reactions = convert_to_reaction_network_with_feedback_enhancement(network)
    
    # Initialize node activities with better baselines
    all_nodes = collect(keys(network.nodes))
    n_nodes = length(all_nodes)
    
    # Smarter baseline initialization
    x0 = Dict{String, Float64}()
    baseline_activities = Dict{String, Float64}()
    
    # Find root and terminal nodes for better baseline assignment
    all_parents = Set(edge.parent_uuid for edge in network.edges)
    all_children = Set(edge.child_uuid for edge in network.edges)
    root_nodes = [uuid for uuid in all_parents if !(uuid in all_children)]
    terminal_nodes = [uuid for uuid in all_children if !(uuid in all_parents)]
    
    for (uuid, node) in network.nodes
        # Assign baselines based on network position
        baseline_val = if uuid in root_nodes
            0.3  # Higher baseline for roots (sources of activity)
        elseif uuid in terminal_nodes
            0.15  # Lower baseline for terminals (sinks)
        else
            0.2   # Medium baseline for intermediates
        end
        
        if haskey(observations, uuid)
            x0[uuid] = observations[uuid][1] / 100.0  # Convert from UI scale
        else
            x0[uuid] = baseline_val
        end
        baseline_activities[uuid] = baseline_val
    end
    
    # Use enhanced fixed-point solver
    return solve_steady_state_feedback_aware(reactions, observations, x0, baseline_activities, params, start_time)
end

"""
Feedback-aware fixed-point solver with adaptive damping.
"""
function solve_steady_state_feedback_aware(
    reactions::Vector{Reaction},
    observations::Dict{String, Tuple{Float64, Float64}},
    x0::Dict{String, Float64},
    baseline_activities::Dict{String, Float64},
    params::SteadyStateParams,
    start_time::Float64
)::SolverResult
    
    x_current = copy(x0)
    damping = 0.4  # Start with moderate damping
    
    # Track convergence history for adaptive damping
    change_history = Float64[]
    
    for iter in 1:params.max_iters
        # Apply forward model
        x_forward = forward_model(x_current, reactions)
        
        # Adaptive damping based on convergence behavior
        if length(change_history) >= 5
            recent_changes = change_history[end-4:end]
            if std(recent_changes) > mean(recent_changes) * 0.5
                # Oscillatory behavior - increase damping
                damping = min(0.7, damping * 1.1)
            elseif all(recent_changes[i] > recent_changes[i+1] for i in 1:4)
                # Monotonic convergence - decrease damping for speed
                damping = max(0.2, damping * 0.95)
            end
        end
        
        # Update with adaptive damping
        x_new = Dict{String, Float64}()
        max_change = 0.0
        
        for uuid in keys(x_current)
            if haskey(x_forward, uuid)
                # Standard fixed-point update
                x_new[uuid] = (1.0 - damping) * x_current[uuid] + damping * x_forward[uuid]
            else
                x_new[uuid] = x_current[uuid]
            end
            
            # Apply observations with feedback-aware blending
            if haskey(observations, uuid)
                obs_value = observations[uuid][1] / 100.0
                confidence = observations[uuid][2]
                
                # For negative feedback systems, use gentler observation constraint
                # to allow the feedback to work
                is_in_feedback = any(r -> uuid in [r.activator_uuids; r.inhibitor_uuids; r.target_uuid], reactions)
                
                blend_factor = if is_in_feedback
                    min(confidence * 0.2, 0.25)  # Very gentle for feedback nodes
                else
                    min(confidence * 0.4, 0.5)   # Normal for other nodes
                end
                
                x_new[uuid] = (1.0 - blend_factor) * x_new[uuid] + blend_factor * obs_value
            end
            
            # Apply baseline prior with feedback awareness
            baseline_pull = params.gamma * 0.5  # Weaker baseline pull
            x_new[uuid] = (1.0 - baseline_pull) * x_new[uuid] + baseline_pull * baseline_activities[uuid]
            
            # Clamp to valid range
            x_new[uuid] = clamp(x_new[uuid], 0.0, 1.0)
            
            # Track convergence
            change = abs(x_new[uuid] - x_current[uuid])
            max_change = max(max_change, change)
        end
        
        push!(change_history, max_change)
        x_current = x_new
        
        # Check convergence
        if max_change < params.tolerance
            solve_time = time() - start_time
            
            return SolverResult(
                x_current,
                true,
                iter,
                max_change,
                solve_time,
                Dict("method" => "feedback_aware", "final_damping" => damping, "avg_damping_changes" => length(change_history) > 1 ? std(change_history) : 0.0)
            )
        end
        
        if iter % 25 == 0
            println("Enhanced iteration $iter, max_change = $max_change, damping = $(round(damping, digits=2))")
        end
    end
    
    # Max iterations reached
    solve_time = time() - start_time
    
    return SolverResult(
        x_current,
        false,
        params.max_iters,
        NaN,
        solve_time,
        Dict("method" => "feedback_aware", "warning" => "max_iterations_reached", "final_damping" => damping)
    )
end

"""
Alternative approach: Quasi-Newton method for feedback systems.
Uses approximated Jacobian to better handle feedback dynamics.
"""
function solve_steady_state_quasi_newton(
    reactions::Vector{Reaction},
    observations::Dict{String, Tuple{Float64, Float64}},
    x0::Dict{String, Float64},
    baseline_activities::Dict{String, Float64},
    params::SteadyStateParams,
    start_time::Float64
)::SolverResult
    
    all_nodes = collect(keys(x0))
    n_nodes = length(all_nodes)
    
    # Convert to vector form
    x_vec = [x0[uuid] for uuid in all_nodes]
    
    # Define the system F(x) = x - G(x) where G is the forward model
    # We want to solve F(x) = 0
    function system_function(x_vec)
        x_dict = Dict(all_nodes[i] => clamp(x_vec[i], 0.0, 1.0) for i in 1:n_nodes)
        x_forward = forward_model(x_dict, reactions)
        
        residual = zeros(n_nodes)
        for i in 1:n_nodes
            uuid = all_nodes[i]
            forward_val = get(x_forward, uuid, x_dict[uuid])
            
            # System equation: x - G(x) = 0, modified with observations and baseline
            residual[i] = x_dict[uuid] - forward_val
            
            # Add observation constraints
            if haskey(observations, uuid)
                obs_value = observations[uuid][1] / 100.0
                confidence = observations[uuid][2]
                obs_weight = confidence * 0.3  # Moderate observation weight
                residual[i] += obs_weight * (x_dict[uuid] - obs_value)
            end
            
            # Add baseline prior
            baseline_weight = params.gamma * 0.1
            residual[i] += baseline_weight * (x_dict[uuid] - baseline_activities[uuid])
        end
        
        return residual
    end
    
    # Simple Newton iteration with finite difference Jacobian
    x_current = copy(x_vec)
    
    for iter in 1:params.max_iters
        # Compute function value
        f_current = system_function(x_current)
        residual_norm = norm(f_current)
        
        if residual_norm < params.tolerance
            # Converged
            solve_time = time() - start_time
            result_dict = Dict(all_nodes[i] => clamp(x_current[i], 0.0, 1.0) for i in 1:n_nodes)
            
            return SolverResult(
                result_dict,
                true,
                iter,
                residual_norm,
                solve_time,
                Dict("method" => "quasi_newton")
            )
        end
        
        # Compute finite difference Jacobian
        h = 1e-6
        J = zeros(n_nodes, n_nodes)
        
        for j in 1:n_nodes
            x_plus = copy(x_current)
            x_plus[j] += h
            f_plus = system_function(x_plus)
            
            J[:, j] = (f_plus - f_current) / h
        end
        
        # Solve linear system J * delta_x = -f_current
        try
            delta_x = J \ (-f_current)
            
            # Line search with damping
            alpha = 1.0
            for _ in 1:5
                x_new = x_current + alpha * delta_x
                x_new = clamp.(x_new, 0.0, 1.0)  # Project to feasible region
                
                f_new = system_function(x_new)
                if norm(f_new) < residual_norm  # Armijo condition (simplified)
                    x_current = x_new
                    break
                end
                alpha *= 0.5
            end
        catch e
            # Jacobian singular - fall back to steepest descent
            x_current = x_current - 0.1 * f_current
            x_current = clamp.(x_current, 0.0, 1.0)
        end
        
        if iter % 25 == 0
            println("Quasi-Newton iteration $iter, residual = $(round(residual_norm, digits=6))")
        end
    end
    
    # Max iterations reached
    solve_time = time() - start_time
    result_dict = Dict(all_nodes[i] => clamp(x_current[i], 0.0, 1.0) for i in 1:n_nodes)
    
    return SolverResult(
        result_dict,
        false,
        params.max_iters,
        NaN,
        solve_time,
        Dict("method" => "quasi_newton", "warning" => "max_iterations_reached")
    )
end