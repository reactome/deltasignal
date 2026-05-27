# Time-Dynamic Solver with Substrate Consumption and Product Accumulation

using Statistics

"""
Parameters for time-dynamic simulation.
"""
struct TimeDynamicParams
    t_start::Float64          # Start time
    t_end::Float64            # End time  
    dt::Float64               # Time step
    tolerance::Float64        # ODE solver tolerance
    max_iters::Int           # Maximum ODE iterations
    
    # Substrate dynamics
    substrate_consumption_rate::Float64  # How fast substrates are consumed
    substrate_replenishment_rate::Float64  # How fast substrates are replenished
    
    # Product dynamics  
    product_accumulation_rate::Float64   # How fast products accumulate
    product_decay_rate::Float64         # How fast products decay
    
    # Observation handling
    observation_weight::Float64         # How strongly to enforce observations over time
    
    function TimeDynamicParams(
        t_start::Float64 = 0.0,
        t_end::Float64 = 10.0,
        dt::Float64 = 0.1,
        tolerance::Float64 = 1e-6,
        max_iters::Int = 1000,
        substrate_consumption_rate::Float64 = 0.1,
        substrate_replenishment_rate::Float64 = 0.05,
        product_accumulation_rate::Float64 = 0.1,
        product_decay_rate::Float64 = 0.02,
        observation_weight::Float64 = 0.3
    )
        new(t_start, t_end, dt, tolerance, max_iters,
            substrate_consumption_rate, substrate_replenishment_rate,
            product_accumulation_rate, product_decay_rate, observation_weight)
    end
end

"""
Result of time-dynamic simulation including full time series.
"""
struct TimeDynamicResult
    time_points::Vector{Float64}
    node_activities::Dict{String, Vector{Float64}}  # node_id -> time series
    substrate_levels::Dict{String, Vector{Float64}}  # substrate_id -> time series  
    product_levels::Dict{String, Vector{Float64}}   # product_id -> time series
    converged::Bool
    final_time::Float64
    solve_time::Float64
    diagnostics::Dict{String, Any}
end

"""
Internal state for time-dynamic simulation.
"""
mutable struct TimeDynamicState
    # Node activities (main state variables)
    x::Dict{String, Float64}
    
    # Substrate levels
    substrates::Dict{String, Float64}
    
    # Product levels  
    products::Dict{String, Float64}
    
    # Track consumption/production rates
    consumption_rates::Dict{String, Float64}
    production_rates::Dict{String, Float64}
    
    function TimeDynamicState(all_nodes::Vector{String}, reactions::Vector{Reaction})
        x = Dict{String, Float64}()
        substrates = Dict{String, Float64}()
        products = Dict{String, Float64}()
        consumption_rates = Dict{String, Float64}()
        production_rates = Dict{String, Float64}()
        
        # Initialize node activities to baseline
        for node in all_nodes
            x[node] = 0.2  # Default baseline
        end
        
        # Initialize substrates and products from reactions
        for reaction in reactions
            for substrate in reaction.substrate_uuids
                substrates[substrate] = get(substrates, substrate, 1.0)  # Full substrate initially
                consumption_rates[substrate] = 0.0
            end
            
            for product in reaction.product_uuids
                products[product] = get(products, product, 0.0)  # No products initially
                production_rates[product] = 0.0
            end
        end
        
        new(x, substrates, products, consumption_rates, production_rates)
    end
end

"""
Main time-dynamic solver with substrate consumption and product accumulation.
"""
function rollout_time_dynamic(
    network::ReactionNetwork,
    observations::Dict{String, Tuple{Float64, Float64}},
    params::TimeDynamicParams = TimeDynamicParams()
)::TimeDynamicResult
    
    start_time = time()
    
    println("🕐 Starting time-dynamic simulation (t=$(params.t_start) → $(params.t_end), dt=$(params.dt))")
    
    # Convert network to reactions
    reactions = convert_to_reaction_network(network)
    all_nodes = collect(keys(network.nodes))
    
    println("📊 Simulating $(length(all_nodes)) nodes, $(length(reactions)) reactions")
    
    # Initialize state
    state = TimeDynamicState(all_nodes, reactions)
    
    # Apply initial observations  
    for (node_id, (activity, confidence)) in observations
        if haskey(state.x, node_id)
            state.x[node_id] = activity / 100.0  # Convert from percentage
        end
    end
    
    # Time integration setup
    time_points = collect(params.t_start:params.dt:params.t_end)
    n_steps = length(time_points)
    
    # Storage for results
    node_timeseries = Dict{String, Vector{Float64}}()
    substrate_timeseries = Dict{String, Vector{Float64}}()
    product_timeseries = Dict{String, Vector{Float64}}()
    
    for node in all_nodes
        node_timeseries[node] = zeros(n_steps)
        node_timeseries[node][1] = state.x[node]
    end
    
    for substrate in keys(state.substrates)
        substrate_timeseries[substrate] = zeros(n_steps)  
        substrate_timeseries[substrate][1] = state.substrates[substrate]
    end
    
    for product in keys(state.products)
        product_timeseries[product] = zeros(n_steps)
        product_timeseries[product][1] = state.products[product]
    end
    
    # Time stepping loop
    converged = true
    
    for step in 2:n_steps
        current_time = time_points[step]
        dt = params.dt
        
        # Compute derivatives
        dx_dt = compute_activity_derivatives(state, reactions, observations, params, current_time)
        ds_dt = compute_substrate_derivatives(state, reactions, params)
        dp_dt = compute_product_derivatives(state, reactions, params)
        
        # Update state (simple forward Euler)
        for node in all_nodes
            state.x[node] = clamp(state.x[node] + dt * dx_dt[node], 0.0, 1.0)
            node_timeseries[node][step] = state.x[node]
        end
        
        for substrate in keys(state.substrates)
            state.substrates[substrate] = clamp(state.substrates[substrate] + dt * ds_dt[substrate], 0.0, 2.0)
            substrate_timeseries[substrate][step] = state.substrates[substrate]
        end
        
        for product in keys(state.products)
            state.products[product] = clamp(state.products[product] + dt * dp_dt[product], 0.0, 5.0)
            product_timeseries[product][step] = state.products[product]
        end
        
        # Check for convergence (steady state reached)
        if step > 10 && step % 10 == 0
            recent_changes = Float64[]
            for node in all_nodes
                if length(node_timeseries[node]) >= 5
                    recent_vals = node_timeseries[node][end-4:end]
                    push!(recent_changes, maximum(recent_vals) - minimum(recent_vals))
                end
            end
            
            if !isempty(recent_changes) && maximum(recent_changes) < params.tolerance * 10
                println("✅ Converged to steady state at t=$(round(current_time, digits=2))")
                # Fill remaining time points with steady state
                for remaining_step in (step+1):n_steps
                    for node in all_nodes
                        node_timeseries[node][remaining_step] = state.x[node]
                    end
                    for substrate in keys(state.substrates)
                        substrate_timeseries[substrate][remaining_step] = state.substrates[substrate] 
                    end
                    for product in keys(state.products)
                        product_timeseries[product][remaining_step] = state.products[product]
                    end
                end
                break
            end
        end
        
        if step % max(1, div(n_steps, 20)) == 0
            println("⏳ Progress: t=$(round(current_time, digits=2)) ($(round(step/n_steps*100, digits=1))%)")
        end
    end
    
    solve_time = time() - start_time
    
    println("🏁 Time-dynamic simulation completed in $(round(solve_time, digits=3))s")
    
    return TimeDynamicResult(
        time_points,
        node_timeseries,
        substrate_timeseries,
        product_timeseries,
        converged,
        time_points[end],
        solve_time,
        Dict("method" => "time_dynamic_euler", "n_steps" => n_steps)
    )
end

"""
Compute time derivatives of node activities.
"""
function compute_activity_derivatives(
    state::TimeDynamicState,
    reactions::Vector{Reaction},
    observations::Dict{String, Tuple{Float64, Float64}},
    params::TimeDynamicParams,
    current_time::Float64
)::Dict{String, Float64}
    
    dx_dt = Dict{String, Float64}()
    
    # Initialize all derivatives to baseline decay
    for node in keys(state.x)
        dx_dt[node] = -0.01 * (state.x[node] - 0.2)  # Decay toward baseline
    end
    
    # Apply reaction dynamics
    for reaction in reactions
        target = reaction.target_uuid
        
        # Compute forward reaction rate (same as steady-state forward model)
        forward_activity = compute_reaction_activity(
            state.x, reaction, state.substrates, state.products
        )
        
        # Add reaction contribution (activity change toward forward value)  
        current_activity = get(state.x, target, 0.2)
        reaction_drive = 2.0 * (forward_activity - current_activity)  # Reaction rate coefficient
        
        dx_dt[target] = get(dx_dt, target, 0.0) + reaction_drive
    end
    
    # Apply observations as time-dependent forcing
    obs_weight = params.observation_weight
    for (node_id, (target_activity, confidence)) in observations
        if haskey(dx_dt, node_id)
            target_val = target_activity / 100.0
            current_val = state.x[node_id]
            observation_force = obs_weight * confidence * (target_val - current_val)
            dx_dt[node_id] += observation_force
        end
    end
    
    return dx_dt
end

"""
Compute time derivatives of substrate levels.
"""
function compute_substrate_derivatives(
    state::TimeDynamicState,
    reactions::Vector{Reaction},
    params::TimeDynamicParams
)::Dict{String, Float64}
    
    ds_dt = Dict{String, Float64}()
    
    # Initialize with replenishment
    for substrate in keys(state.substrates)
        # Replenishment toward full level (1.0)
        replenishment = params.substrate_replenishment_rate * (1.0 - state.substrates[substrate])
        ds_dt[substrate] = replenishment
    end
    
    # Apply consumption from active reactions
    for reaction in reactions
        if !isempty(reaction.substrate_uuids)
            # Compute reaction activity level
            target_activity = get(state.x, reaction.target_uuid, 0.2)
            
            # Consumption proportional to reaction activity and substrate availability
            for substrate in reaction.substrate_uuids
                consumption_rate = params.substrate_consumption_rate * target_activity * 
                                state.substrates[substrate] * reaction.params.consumption_lambdas[1]  # Use first lambda
                
                ds_dt[substrate] = get(ds_dt, substrate, 0.0) - consumption_rate
                state.consumption_rates[substrate] = consumption_rate
            end
        end
    end
    
    return ds_dt
end

"""
Compute time derivatives of product levels.
"""
function compute_product_derivatives(
    state::TimeDynamicState,
    reactions::Vector{Reaction},
    params::TimeDynamicParams
)::Dict{String, Float64}
    
    dp_dt = Dict{String, Float64}()
    
    # Initialize with decay
    for product in keys(state.products)
        decay = params.product_decay_rate * state.products[product]
        dp_dt[product] = -decay
    end
    
    # Apply production from active reactions  
    for reaction in reactions
        if !isempty(reaction.product_uuids)
            # Compute reaction activity level
            target_activity = get(state.x, reaction.target_uuid, 0.2)
            
            # Production proportional to reaction activity
            for product in reaction.product_uuids
                production_rate = params.product_accumulation_rate * target_activity
                
                dp_dt[product] = get(dp_dt, product, 0.0) + production_rate
                state.production_rates[product] = production_rate
            end
        end
    end
    
    return dp_dt
end

"""
Compute single reaction activity with substrate/product effects.
"""
function compute_reaction_activity(
    x::Dict{String, Float64},
    reaction::Reaction,
    substrates::Dict{String, Float64},
    products::Dict{String, Float64}
)::Float64
    
    # Base activity from activators/inhibitors (same as steady-state)
    base_activity = compute_reaction_output(reaction, x)
    
    # Substrate availability multiplier
    substrate_multiplier = 1.0
    if !isempty(reaction.substrate_uuids)
        substrate_levels = [get(substrates, sub, 1.0) for sub in reaction.substrate_uuids]
        substrate_multiplier = minimum(substrate_levels)  # Limiting substrate
    end
    
    # Product inhibition (negative feedback from products)
    product_inhibition = 1.0
    if !isempty(reaction.product_uuids)
        product_levels = [get(products, prod, 0.0) for prod in reaction.product_uuids]
        max_product = maximum(vcat(product_levels, [0.0]))
        product_inhibition = 1.0 / (1.0 + max_product)  # Product inhibition
    end
    
    return base_activity * substrate_multiplier * product_inhibition
end

"""
Default parameters for time-dynamic simulation.
"""
function default_time_dynamic_params()::TimeDynamicParams
    return TimeDynamicParams()
end