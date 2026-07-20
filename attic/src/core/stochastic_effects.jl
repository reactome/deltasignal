# Stochastic Effects Module
# Models biological noise, cell-to-cell variability, and transcriptional bursting

using Random
using Statistics

"""
Types of biological noise sources.
"""
@enum NoiseType begin
    TRANSCRIPTIONAL_BURSTING    # Gene expression bursts
    POISSON_NOISE              # Low copy number fluctuations
    THERMAL_NOISE              # Molecular motion effects
    MEASUREMENT_NOISE          # Experimental measurement error
    EXTRINSIC_NOISE           # Environmental fluctuations
    INTRINSIC_NOISE           # Inherent biochemical randomness
end

"""
Stochastic process parameters.
"""
struct StochasticParameters
    noise_type::NoiseType
    amplitude::Float64          # Noise strength
    correlation_time::Float64   # Correlation time scale (minutes)
    burst_frequency::Float64    # For transcriptional bursting (events/hour)
    burst_size_mean::Float64    # Mean burst size
    burst_size_variance::Float64 # Variance in burst size
end

"""
Cell population heterogeneity model.
"""
struct PopulationHeterogeneity
    cell_cycle_variability::Float64    # Cell cycle effects
    protein_expression_cv::Float64     # Coefficient of variation in expression
    metabolic_state_variance::Float64  # Metabolic heterogeneity
    age_dependent_factor::Float64      # Cell age effects
    spatial_gradient::Float64          # Position effects in tissue
end

"""
Correlated noise generator for pathway components.
"""
mutable struct CorrelatedNoiseGenerator
    correlation_matrix::Matrix{Float64}
    noise_history::Matrix{Float64}
    time_step::Float64
    correlation_length::Float64
end

"""
Get realistic stochastic parameters for different biological processes.
"""
function get_stochastic_parameters(process_type::ProcessType)::StochasticParameters
    
    if process_type == TRANSCRIPTION
        # Transcriptional bursting is well-documented
        return StochasticParameters(
            TRANSCRIPTIONAL_BURSTING,
            0.8,    # High amplitude - gene expression is noisy
            15.0,   # ~15 minute correlation time
            4.0,    # ~4 bursts per hour
            50.0,   # Mean burst size (mRNA molecules)
            25.0    # High variance in burst size
        )
        
    elseif process_type == PROTEIN_SYNTHESIS
        # Translation noise (lower than transcription)
        return StochasticParameters(
            POISSON_NOISE,
            0.4,    # Moderate amplitude
            5.0,    # ~5 minute correlation
            0.0,    # No bursting at translation level
            0.0, 0.0  # Not applicable
        )
        
    elseif process_type == ENZYME_CATALYSIS
        # Thermal fluctuations dominate
        return StochasticParameters(
            THERMAL_NOISE,
            0.15,   # Low amplitude - enzymes are stable
            0.1,    # Very short correlation (~6 seconds)
            0.0, 0.0, 0.0  # Not applicable
        )
        
    elseif process_type == SIGNALING_CASCADE
        # Mix of intrinsic and extrinsic noise
        return StochasticParameters(
            INTRINSIC_NOISE,
            0.25,   # Moderate amplitude
            2.0,    # ~2 minute correlation
            0.0, 0.0, 0.0
        )
        
    elseif process_type == METABOLIC_FLUX
        # Low noise due to homeostatic mechanisms
        return StochasticParameters(
            THERMAL_NOISE,
            0.1,    # Very low amplitude
            1.0,    # ~1 minute correlation
            0.0, 0.0, 0.0
        )
        
    else
        # Default moderate noise
        return StochasticParameters(
            INTRINSIC_NOISE,
            0.2,    # Default amplitude
            1.0,    # Default correlation
            0.0, 0.0, 0.0
        )
    end
end

"""
Generate transcriptional bursting noise.
"""
function generate_transcriptional_burst(
    params::StochasticParameters,
    current_time::Float64,
    dt::Float64
)::Float64
    
    # Probability of burst in this time step
    burst_prob = params.burst_frequency * (dt / 60.0)  # Convert to per-minute
    
    if rand() < burst_prob
        # Generate burst
        burst_size = max(0.0, params.burst_size_mean + 
                        params.burst_size_variance * randn())
        
        # Convert to normalized activity boost
        normalized_burst = burst_size / 100.0  # Assume max ~100 molecules -> activity 1.0
        return clamp(normalized_burst, 0.0, 1.0)
    else
        return 0.0
    end
end

"""
Generate Ornstein-Uhlenbeck colored noise.
"""
function generate_ou_noise(
    previous_noise::Float64,
    params::StochasticParameters,
    dt::Float64
)::Float64
    
    # Ornstein-Uhlenbeck process: dx = -x/τ dt + σ√(2/τ) dW
    τ = params.correlation_time
    σ = params.amplitude
    
    # Exact solution for discrete time step
    exp_factor = exp(-dt / τ)
    
    noise = previous_noise * exp_factor + 
            σ * sqrt(1.0 - exp_factor^2) * randn()
    
    return noise
end

"""
Generate Poisson noise for low copy numbers.
"""
function generate_poisson_noise(
    mean_activity::Float64,
    params::StochasticParameters
)::Float64
    
    # Convert activity to molecule count (assume max activity = 1000 molecules)
    molecule_count = mean_activity * 1000.0
    
    if molecule_count < 50.0  # Low copy number regime
        # Poisson fluctuations
        actual_count = max(0.0, molecule_count + sqrt(molecule_count) * randn())
        return actual_count / 1000.0
    else
        # High copy number - use Gaussian approximation
        return mean_activity + params.amplitude * randn() * mean_activity
    end
end

"""
Create correlated noise generator for pathway components.
"""
function create_correlated_noise_generator(
    node_ids::Vector{String},
    network::ReactionNetwork;
    base_correlation::Float64 = 0.3
)::CorrelatedNoiseGenerator
    
    n_nodes = length(node_ids)
    correlation_matrix = Matrix{Float64}(I, n_nodes, n_nodes)
    
    # Set correlations based on network topology
    for i in 1:n_nodes
        for j in (i+1):n_nodes
            node_i = node_ids[i]
            node_j = node_ids[j]
            
            # Check if nodes are connected
            connected = false
            for edge in network.edges
                if (edge.parent_uuid == node_i && edge.child_uuid == node_j) ||
                   (edge.parent_uuid == node_j && edge.child_uuid == node_i)
                    connected = true
                    break
                end
            end
            
            # Set correlation
            if connected
                correlation = base_correlation + 0.2  # Higher correlation for connected nodes
            else
                correlation = base_correlation * 0.3  # Lower for unconnected
            end
            
            correlation_matrix[i, j] = correlation
            correlation_matrix[j, i] = correlation
        end
    end
    
    # Initialize noise history
    history_length = 10
    noise_history = randn(n_nodes, history_length) * 0.1
    
    return CorrelatedNoiseGenerator(
        correlation_matrix,
        noise_history,
        1.0,  # Default time step
        5.0   # Default correlation length
    )
end

"""
Generate correlated noise for multiple nodes.
"""
function generate_correlated_noise!(
    noise_generator::CorrelatedNoiseGenerator,
    current_time::Float64
)::Vector{Float64}
    
    n_nodes = size(noise_generator.noise_history, 1)
    
    # Generate independent Gaussian noise
    independent_noise = randn(n_nodes)
    
    # Apply correlation matrix (Cholesky decomposition for efficiency)
    try
        L = cholesky(noise_generator.correlation_matrix).L
        correlated_noise = L * independent_noise
    catch
        # If correlation matrix is not positive definite, use independent noise
        correlated_noise = independent_noise
    end
    
    # Apply temporal correlation (simple AR(1) model)
    α = exp(-noise_generator.time_step / noise_generator.correlation_length)
    
    if size(noise_generator.noise_history, 2) > 1
        previous_noise = noise_generator.noise_history[:, end]
        temporal_correlated_noise = α * previous_noise + sqrt(1 - α^2) * correlated_noise
    else
        temporal_correlated_noise = correlated_noise
    end
    
    # Update history
    noise_generator.noise_history = hcat(noise_generator.noise_history[:, 2:end], temporal_correlated_noise)
    
    return temporal_correlated_noise
end

"""
Apply population heterogeneity to base activities.
"""
function apply_population_heterogeneity(
    base_activities::Dict{String, Float64},
    heterogeneity::PopulationHeterogeneity,
    cell_id::Int = 1
)::Dict{String, Float64}
    
    # Set random seed based on cell ID for reproducible heterogeneity
    Random.seed!(cell_id * 12345)
    
    heterogeneous_activities = Dict{String, Float64}()
    
    for (node_id, activity) in base_activities
        # Protein expression variability (log-normal distribution)
        log_mean = log(activity + 0.01)  # Avoid log(0)
        log_std = heterogeneity.protein_expression_cv
        
        varied_activity = exp(log_mean + log_std * randn())
        
        # Cell cycle effects (some proteins vary with cell cycle)
        cell_cycle_phase = (cell_id % 100) / 100.0  # Pseudo cell cycle position
        cycle_factor = 1.0 + heterogeneity.cell_cycle_variability * sin(2π * cell_cycle_phase)
        
        # Metabolic state effects
        metabolic_factor = 1.0 + heterogeneity.metabolic_state_variance * randn()
        
        # Age-dependent effects (older cells may have different expression)
        cell_age_factor = 1.0 - heterogeneity.age_dependent_factor * (cell_id % 50) / 50.0
        
        # Combine all factors
        final_activity = varied_activity * cycle_factor * metabolic_factor * cell_age_factor
        
        # Clamp to valid range
        heterogeneous_activities[node_id] = clamp(final_activity, 0.0, 1.0)
    end
    
    # Reset random seed
    Random.seed!()
    
    return heterogeneous_activities
end

"""
Stochastic simulation with realistic biological noise.
"""
function simulate_with_stochastic_effects(
    network::ReactionNetwork,
    initial_conditions::Dict{String, Float64},
    time_span::Tuple{Float64, Float64};
    n_cells::Int = 1,
    enable_bursting::Bool = true,
    enable_correlation::Bool = true,
    population_heterogeneity::PopulationHeterogeneity = PopulationHeterogeneity(0.2, 0.3, 0.1, 0.05, 0.02)
)::Dict{String, Any}
    
    println("🎲 Running stochastic simulation with biological noise...")
    println("Number of cells: $n_cells")
    println("Time span: $(time_span[1]) to $(time_span[2]) minutes")
    
    # Set up reactions and biological context
    reactions = create_biologically_enhanced_reaction_network(network)
    biological_context = infer_biological_processes(network)
    
    # Node setup
    node_ids = collect(keys(network.nodes))
    n_nodes = length(node_ids)
    
    # Time setup
    dt = 0.1  # 6 seconds time step
    t_eval = range(time_span[1], time_span[2], step=dt)
    n_steps = length(t_eval)
    
    # Create noise generators
    noise_generators = Dict{String, StochasticParameters}()
    for (node_id, process_type) in biological_context
        noise_generators[node_id] = get_stochastic_parameters(process_type)
    end
    
    # Correlated noise generator
    corr_noise_gen = enable_correlation ? create_correlated_noise_generator(node_ids, network) : nothing
    
    # Storage for results
    cell_trajectories = Dict{Int, Dict{String, Vector{Float64}}}()
    
    # Simulate each cell
    for cell_id in 1:n_cells
        if cell_id % max(1, div(n_cells, 10)) == 0
            println("  Simulating cell $cell_id / $n_cells")
        end
        
        # Apply population heterogeneity to initial conditions
        cell_initial = apply_population_heterogeneity(initial_conditions, population_heterogeneity, cell_id)
        
        # Initialize trajectory storage
        cell_trajectories[cell_id] = Dict{String, Vector{Float64}}()
        for node_id in node_ids
            cell_trajectories[cell_id][node_id] = Float64[]
        end
        
        # Current state
        current_activities = copy(cell_initial)
        noise_states = Dict{String, Float64}()  # For colored noise
        
        # Time evolution
        for (step, t) in enumerate(t_eval)
            # Generate correlated noise if enabled
            correlated_noise = enable_correlation && corr_noise_gen !== nothing ? 
                              generate_correlated_noise!(corr_noise_gen, t) : zeros(n_nodes)
            
            # Update each node
            updated_activities = Dict{String, Float64}()
            
            for (node_idx, node_id) in enumerate(node_ids)
                # Get deterministic activity from reactions
                deterministic_activity = 0.0
                reaction_count = 0
                
                for reaction in reactions
                    if reaction.target_uuid == node_id
                        deterministic_activity += compute_enhanced_reaction_output(reaction, current_activities)
                        reaction_count += 1
                    end
                end
                
                if reaction_count > 0
                    deterministic_activity /= reaction_count
                else
                    deterministic_activity = get(current_activities, node_id, 0.1)
                end
                
                # Add stochastic effects
                total_noise = 0.0
                
                if haskey(noise_generators, node_id)
                    noise_params = noise_generators[node_id]
                    
                    # Generate noise based on type
                    if noise_params.noise_type == TRANSCRIPTIONAL_BURSTING && enable_bursting
                        burst_noise = generate_transcriptional_burst(noise_params, t, dt)
                        total_noise += burst_noise
                        
                    elseif noise_params.noise_type == POISSON_NOISE
                        poisson_component = generate_poisson_noise(deterministic_activity, noise_params)
                        total_noise += poisson_component - deterministic_activity
                        
                    else
                        # Colored noise (Ornstein-Uhlenbeck)
                        previous_noise = get(noise_states, node_id, 0.0)
                        ou_noise = generate_ou_noise(previous_noise, noise_params, dt)
                        noise_states[node_id] = ou_noise
                        total_noise += ou_noise
                    end
                end
                
                # Add correlated noise component
                if enable_correlation
                    total_noise += 0.1 * correlated_noise[node_idx]  # Scale correlated component
                end
                
                # Combine deterministic + stochastic
                noisy_activity = deterministic_activity + total_noise
                
                # Simple dynamics (first-order)
                current_val = get(current_activities, node_id, 0.1)
                τ = 2.0  # ~2 minute time constant
                new_activity = current_val + dt * (noisy_activity - current_val) / τ
                
                # Clamp to valid range
                updated_activities[node_id] = clamp(new_activity, 0.0, 1.0)
            end
            
            # Update current state
            current_activities = updated_activities
            
            # Store trajectory
            for node_id in node_ids
                push!(cell_trajectories[cell_id][node_id], current_activities[node_id])
            end
        end
    end
    
    return Dict(
        "time_points" => collect(t_eval),
        "cell_trajectories" => cell_trajectories,
        "node_ids" => node_ids,
        "simulation_parameters" => Dict(
            "n_cells" => n_cells,
            "dt" => dt,
            "enable_bursting" => enable_bursting,
            "enable_correlation" => enable_correlation,
            "population_heterogeneity" => population_heterogeneity
        )
    )
end

"""
Analyze cell-to-cell variability from stochastic simulation.
"""
function analyze_cell_variability(
    stochastic_result::Dict{String, Any}
)::Dict{String, Any}
    
    time_points = stochastic_result["time_points"]
    cell_trajectories = stochastic_result["cell_trajectories"]
    node_ids = stochastic_result["node_ids"]
    n_cells = length(cell_trajectories)
    
    analysis = Dict{String, Any}()
    
    # For each time point, compute population statistics
    variability_over_time = Dict{String, Vector{Float64}}()
    
    for node_id in node_ids
        cv_over_time = Float64[]
        
        for t_idx in 1:length(time_points)
            # Collect activities across all cells at this time point
            activities_at_t = Float64[]
            
            for cell_id in 1:n_cells
                if haskey(cell_trajectories[cell_id], node_id) && 
                   t_idx <= length(cell_trajectories[cell_id][node_id])
                    push!(activities_at_t, cell_trajectories[cell_id][node_id][t_idx])
                end
            end
            
            # Compute coefficient of variation
            if length(activities_at_t) > 1 && mean(activities_at_t) > 0.01
                cv = std(activities_at_t) / mean(activities_at_t)
                push!(cv_over_time, cv)
            else
                push!(cv_over_time, 0.0)
            end
        end
        
        variability_over_time[node_id] = cv_over_time
    end
    
    analysis["variability_over_time"] = variability_over_time
    
    # Summary statistics
    mean_variability = Dict{String, Float64}()
    max_variability = Dict{String, Float64}()
    
    for (node_id, cv_trace) in variability_over_time
        mean_variability[node_id] = mean(cv_trace)
        max_variability[node_id] = maximum(cv_trace)
    end
    
    analysis["mean_variability"] = mean_variability
    analysis["max_variability"] = max_variability
    
    # Population correlations
    if n_cells > 2
        final_time_activities = Dict{String, Vector{Float64}}()
        
        for node_id in node_ids
            activities = Float64[]
            for cell_id in 1:n_cells
                if haskey(cell_trajectories[cell_id], node_id) && 
                   !isempty(cell_trajectories[cell_id][node_id])
                    push!(activities, cell_trajectories[cell_id][node_id][end])
                end
            end
            final_time_activities[node_id] = activities
        end
        
        # Compute correlation matrix
        node_list = collect(keys(final_time_activities))
        correlation_matrix = zeros(length(node_list), length(node_list))
        
        for (i, node_i) in enumerate(node_list)
            for (j, node_j) in enumerate(node_list)
                if i != j && length(final_time_activities[node_i]) == length(final_time_activities[node_j])
                    correlation_matrix[i, j] = cor(final_time_activities[node_i], final_time_activities[node_j])
                elseif i == j
                    correlation_matrix[i, j] = 1.0
                end
            end
        end
        
        analysis["population_correlations"] = Dict(
            "nodes" => node_list,
            "correlation_matrix" => correlation_matrix
        )
    end
    
    return analysis
end