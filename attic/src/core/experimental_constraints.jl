# Experimental Data Integration Module
# Constrains parameters with real biological measurements

using Statistics

"""
Experimental data types from different measurement techniques.
"""
@enum ExperimentType begin
    PROTEIN_ABUNDANCE     # Mass spectrometry, Western blot
    BINDING_AFFINITY      # Surface plasmon resonance, ITC
    ENZYME_KINETICS       # In vitro enzyme assays
    KNOCKOUT_EFFECT       # CRISPR screens, gene deletion
    DRUG_RESPONSE         # Dose-response curves
    TIME_COURSE          # Time-series measurements
    PHOSPHORYLATION      # Phosphoproteomics
    LOCALIZATION         # Microscopy, fractionation
end

"""
Experimental measurement with confidence intervals.
"""
struct ExperimentalMeasurement
    value::Float64                    # Measured value
    std_error::Float64               # Standard error
    confidence_interval::Tuple{Float64, Float64}  # 95% CI
    sample_size::Int                 # Number of replicates
    measurement_technique::String    # Method used
    experimental_conditions::String # Cell type, treatment, etc.
    publication_pmid::String         # PubMed ID for reference
end

"""
Comprehensive experimental dataset.
"""
struct ExperimentalDataset
    protein_abundances::Dict{String, ExperimentalMeasurement}
    binding_affinities::Dict{String, ExperimentalMeasurement}  # Kd values
    enzyme_kinetics::Dict{String, NamedTuple{(:km, :vmax), Tuple{ExperimentalMeasurement, ExperimentalMeasurement}}}
    knockout_effects::Dict{String, ExperimentalMeasurement}    # Fold-change in target
    drug_responses::Dict{String, Vector{Tuple{Float64, ExperimentalMeasurement}}}  # Dose-response pairs
    time_courses::Dict{String, Vector{Tuple{Float64, ExperimentalMeasurement}}}    # Time-activity pairs
    phosphorylation_sites::Dict{String, Vector{ExperimentalMeasurement}}
    subcellular_localization::Dict{String, Dict{Compartment, Float64}}
end

"""
Convert experimental binding affinity (Kd) to Hill function threshold (K).
"""
function kd_to_hill_threshold(kd_measurement::ExperimentalMeasurement, max_concentration::Float64 = 1.0)::Float64
    # Kd is typically in nM or μM, convert to normalized [0,1] scale
    # Assume max_concentration represents the scale (e.g., 1 μM = 1.0)
    
    kd_value = kd_measurement.value
    
    # Hill threshold K should be around Kd for half-maximal response
    # But need to normalize to [0,1] scale
    normalized_k = kd_value / max_concentration
    
    # Clamp to reasonable range
    return clamp(normalized_k, 0.01, 0.9)
end

"""
Convert enzyme kinetics (Km, Vmax) to reaction parameters.
"""
function enzyme_kinetics_to_params(
    km_measurement::ExperimentalMeasurement,
    vmax_measurement::ExperimentalMeasurement
)::NamedTuple{(:K, :h), Tuple{Float64, Float64}}
    
    # Km maps to Hill threshold K (half-maximal concentration)
    K = clamp(km_measurement.value / 1000.0, 0.05, 0.8)  # Assume μM scale
    
    # Vmax influences Hill coefficient (higher Vmax -> more cooperative)
    # This is a heuristic - could be refined with more data
    h_base = 2.0
    vmax_factor = log10(max(vmax_measurement.value, 0.1))
    h = clamp(h_base + vmax_factor * 0.5, 1.2, 8.0)
    
    return (K=K, h=h)
end

"""
Convert knockout effect to influence score constraint.
"""
function knockout_to_influence_constraint(knockout_measurement::ExperimentalMeasurement)::Float64
    # Knockout effect as fold-change maps to influence score
    fold_change = abs(knockout_measurement.value)
    
    # Log-scale conversion with confidence weighting
    confidence_weight = 1.0 / (1.0 + knockout_measurement.std_error)
    base_influence = log2(max(fold_change, 1.1)) * 10.0  # Scale to typical influence range
    
    return base_influence * confidence_weight
end

"""
Parse drug dose-response curve to inhibition parameters.
"""
function dose_response_to_inhibition_params(
    dose_response::Vector{Tuple{Float64, ExperimentalMeasurement}}
)::NamedTuple{(:ic50, :hill_slope, :max_inhibition), Tuple{Float64, Float64, Float64}}
    
    if length(dose_response) < 3
        # Insufficient data, return defaults
        return (ic50=0.1, hill_slope=2.0, max_inhibition=0.9)
    end
    
    doses = [pair[1] for pair in dose_response]
    responses = [pair[2].value for pair in dose_response]
    
    # Simple Hill curve fitting (in practice, would use more sophisticated fitting)
    # Find IC50 (dose giving 50% inhibition)
    max_response = maximum(responses)
    min_response = minimum(responses)
    half_response = min_response + (max_response - min_response) / 2
    
    # Find dose closest to half response
    ic50_idx = argmin(abs.(responses .- half_response))
    ic50 = doses[ic50_idx]
    
    # Estimate Hill slope from steepest part of curve
    if length(responses) > 2
        slopes = Float64[]
        for i in 2:length(responses)-1
            slope = (responses[i+1] - responses[i-1]) / (doses[i+1] - doses[i-1])
            push!(slopes, abs(slope))
        end
        hill_slope = mean(slopes) * 2.0  # Scaling factor
    else
        hill_slope = 2.0
    end
    
    max_inhibition = (max_response - min_response) / max_response
    
    return (
        ic50=clamp(ic50/1000.0, 0.01, 1.0),  # Normalize to [0,1] scale
        hill_slope=clamp(hill_slope, 1.0, 6.0),
        max_inhibition=clamp(max_inhibition, 0.1, 1.0)
    )
end

"""
Create parameter constraints from experimental dataset.
"""
function create_experimental_constraints(
    network::ReactionNetwork,
    exp_data::ExperimentalDataset
)::Dict{String, NamedTuple}
    
    println("🔬 Creating experimental parameter constraints...")
    
    constraints = Dict{String, NamedTuple}()
    constraint_count = 0
    
    # Process each node in the network
    for (node_id, node) in network.nodes
        node_constraints = []
        
        # Protein abundance constrains baseline activity
        if haskey(exp_data.protein_abundances, node_id)
            abundance = exp_data.protein_abundances[node_id]
            # Higher abundance -> higher baseline activity (log scale)
            baseline_activity = clamp(log10(max(abundance.value, 0.01)) / 3.0 + 0.5, 0.0, 1.0)
            push!(node_constraints, (:baseline_activity, baseline_activity, abundance.std_error))
        end
        
        # Binding affinity constrains Hill threshold
        if haskey(exp_data.binding_affinities, node_id)
            kd = exp_data.binding_affinities[node_id]
            K = kd_to_hill_threshold(kd)
            push!(node_constraints, (:hill_threshold, K, kd.std_error / kd.value))
        end
        
        # Enzyme kinetics constrains multiple parameters
        if haskey(exp_data.enzyme_kinetics, node_id)
            kinetics = exp_data.enzyme_kinetics[node_id]
            params = enzyme_kinetics_to_params(kinetics.km, kinetics.vmax)
            
            push!(node_constraints, (:hill_threshold, params.K, kinetics.km.std_error / kinetics.km.value))
            push!(node_constraints, (:hill_coefficient, params.h, kinetics.vmax.std_error / kinetics.vmax.value))
        end
        
        # Knockout effects constrain influence scores
        if haskey(exp_data.knockout_effects, node_id)
            knockout = exp_data.knockout_effects[node_id]
            influence = knockout_to_influence_constraint(knockout)
            push!(node_constraints, (:target_influence, influence, knockout.std_error))
        end
        
        # Drug response constrains inhibition parameters
        if haskey(exp_data.drug_responses, node_id)
            dose_response = exp_data.drug_responses[node_id]
            inhibition_params = dose_response_to_inhibition_params(dose_response)
            
            push!(node_constraints, (:inhibition_ic50, inhibition_params.ic50, 0.1))
            push!(node_constraints, (:inhibition_hill_slope, inhibition_params.hill_slope, 0.2))
            push!(node_constraints, (:max_inhibition, inhibition_params.max_inhibition, 0.1))
        end
        
        # Store constraints for this node if any exist
        if !isempty(node_constraints)
            constraints[node_id] = (constraints=node_constraints,)
            constraint_count += length(node_constraints)
        end
    end
    
    println("📊 Applied $(constraint_count) experimental constraints to $(length(constraints)) nodes")
    return constraints
end

"""
Apply experimental constraints to reaction parameters with confidence weighting.
"""
function apply_experimental_constraints(
    base_params::ReactionParams,
    node_id::String,
    constraints::Dict{String, NamedTuple},
    confidence_weight::Float64 = 0.7
)::ReactionParams
    
    if !haskey(constraints, node_id)
        return base_params
    end
    
    # Start with base parameters
    h = base_params.h
    K = base_params.K
    inhibitor_betas = copy(base_params.inhibitor_betas)
    
    # Apply each constraint with confidence weighting
    for constraint in constraints[node_id].constraints
        constraint_type = constraint[1]
        constraint_value = constraint[2]
        constraint_uncertainty = constraint[3]
        
        # Weight based on experimental confidence
        weight = confidence_weight / (1.0 + constraint_uncertainty)
        
        if constraint_type == :hill_threshold
            # Weighted average between base and experimental value
            K = K * (1 - weight) + constraint_value * weight
            
        elseif constraint_type == :hill_coefficient
            h = h * (1 - weight) + constraint_value * weight
            
        elseif constraint_type == :inhibition_ic50
            # Adjust inhibition strength based on IC50
            if !isempty(inhibitor_betas)
                # Lower IC50 -> stronger inhibition
                ic50_factor = 0.1 / max(constraint_value, 0.01)
                inhibitor_betas = inhibitor_betas .* (1 - weight) .+ (inhibitor_betas .* ic50_factor) .* weight
            end
        end
    end
    
    # Ensure parameters stay in valid ranges
    h = clamp(h, 1.0, 10.0)
    K = clamp(K, 0.01, 0.95)
    inhibitor_betas = clamp.(inhibitor_betas, 1.0, 50.0)
    
    return ReactionParams(
        h, K,
        base_params.activator_weights,
        base_params.activator_sensitivity_s,
        base_params.activator_sensitivity_n,
        base_params.activator_sensitivity_K,
        inhibitor_betas,
        base_params.inhibitor_ms,
        base_params.substrate_weights,
        base_params.consumption_lambdas,
        base_params.production_etas,
        base_params.replenishment_rho,
        base_params.decay_delta
    )
end

"""
Create experimentally-constrained reaction network.
"""
function create_experimentally_constrained_network(
    network::ReactionNetwork,
    exp_data::ExperimentalDataset;
    confidence_weight::Float64 = 0.7
)::Vector{Reaction}
    
    println("🧪 Creating experimentally-constrained reaction network...")
    
    # Get experimental constraints
    constraints = create_experimental_constraints(network, exp_data)
    
    # Start with biologically enhanced reactions
    base_reactions = create_biologically_enhanced_reaction_network(network)
    
    # Apply experimental constraints to each reaction
    constrained_reactions = Reaction[]
    
    for reaction in base_reactions
        # Apply constraints to target node
        constrained_params = apply_experimental_constraints(
            reaction.params,
            reaction.target_uuid,
            constraints,
            confidence_weight
        )
        
        constrained_reaction = Reaction(
            reaction.target_uuid,
            reaction.activator_uuids,
            reaction.inhibitor_uuids,
            reaction.depletion_uuids,
            reaction.substrate_uuids,
            reaction.product_uuids,
            constrained_params,
            reaction.is_and_gate,
            reaction.activator_is_and,
            reaction.inhibitor_is_and,
        )
        
        push!(constrained_reactions, constrained_reaction)
    end
    
    return constrained_reactions
end

"""
Load experimental data from various file formats.
"""
function load_experimental_dataset(
    abundance_file::Union{String, Nothing} = nothing,
    kinetics_file::Union{String, Nothing} = nothing,
    knockout_file::Union{String, Nothing} = nothing,
    drug_response_file::Union{String, Nothing} = nothing
)::ExperimentalDataset
    
    # Initialize empty dataset
    dataset = ExperimentalDataset(
        Dict{String, ExperimentalMeasurement}(),  # protein_abundances
        Dict{String, ExperimentalMeasurement}(),  # binding_affinities
        Dict{String, NamedTuple{(:km, :vmax), Tuple{ExperimentalMeasurement, ExperimentalMeasurement}}}(),  # enzyme_kinetics
        Dict{String, ExperimentalMeasurement}(),  # knockout_effects
        Dict{String, Vector{Tuple{Float64, ExperimentalMeasurement}}}(),  # drug_responses
        Dict{String, Vector{Tuple{Float64, ExperimentalMeasurement}}}(),  # time_courses
        Dict{String, Vector{ExperimentalMeasurement}}(),  # phosphorylation_sites
        Dict{String, Dict{Compartment, Float64}}()  # subcellular_localization
    )
    
    # Load protein abundances (simple CSV format)
    if abundance_file !== nothing && isfile(abundance_file)
        println("📊 Loading protein abundance data from $abundance_file")
        # Would implement CSV parsing here
        # Format: protein_id, abundance_value, std_error
    end
    
    # Load enzyme kinetics (CSV format)
    if kinetics_file !== nothing && isfile(kinetics_file)
        println("⚗️ Loading enzyme kinetics data from $kinetics_file")
        # Would implement CSV parsing here
        # Format: enzyme_id, km_value, km_error, vmax_value, vmax_error
    end
    
    # Load knockout effects
    if knockout_file !== nothing && isfile(knockout_file)
        println("🔬 Loading knockout effect data from $knockout_file")
        # Would implement parsing here
        # Format: gene_id, fold_change, p_value, cell_type
    end
    
    # Load drug response curves
    if drug_response_file !== nothing && isfile(drug_response_file)
        println("💊 Loading drug response data from $drug_response_file")
        # Would implement parsing here
        # Format: compound_id, target_id, dose, response, error
    end
    
    return dataset
end

"""
Validate experimental constraints against network structure.
"""
function validate_experimental_constraints(
    network::ReactionNetwork,
    constraints::Dict{String, NamedTuple}
)::Dict{String, Any}
    
    validation_report = Dict{String, Any}()
    
    # Check coverage
    network_nodes = keys(network.nodes)
    constrained_nodes = keys(constraints)
    
    coverage = length(constrained_nodes) / length(network_nodes)
    validation_report["coverage"] = coverage
    validation_report["constrained_nodes"] = length(constrained_nodes)
    validation_report["total_nodes"] = length(network_nodes)
    
    # Check constraint types distribution
    constraint_types = Dict{Symbol, Int}()
    for (node_id, node_constraints) in constraints
        for constraint in node_constraints.constraints
            constraint_type = constraint[1]
            constraint_types[constraint_type] = get(constraint_types, constraint_type, 0) + 1
        end
    end
    
    validation_report["constraint_types"] = constraint_types
    
    # Check parameter ranges
    parameter_ranges = Dict{String, Tuple{Float64, Float64}}()
    
    for (node_id, node_constraints) in constraints
        for constraint in node_constraints.constraints
            constraint_type = string(constraint[1])
            constraint_value = constraint[2]
            
            if haskey(parameter_ranges, constraint_type)
                current_min, current_max = parameter_ranges[constraint_type]
                parameter_ranges[constraint_type] = (min(current_min, constraint_value), max(current_max, constraint_value))
            else
                parameter_ranges[constraint_type] = (constraint_value, constraint_value)
            end
        end
    end
    
    validation_report["parameter_ranges"] = parameter_ranges
    
    return validation_report
end