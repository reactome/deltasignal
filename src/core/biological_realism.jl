# Enhanced Biological Realism Module
# Pathway-specific parameters and competitive inhibition improvements

using Statistics

"""
Biological pathway types with different kinetic characteristics.
"""
@enum PathwayType begin
    SIGNALING          # Fast, switch-like, amplification
    METABOLIC          # Steady flux, enzyme kinetics
    TRANSCRIPTIONAL    # Slow, highly cooperative
    STRESS_RESPONSE    # Threshold-based, all-or-nothing
    CELL_CYCLE        # Ordered, checkpoint-controlled
    APOPTOSIS         # Irreversible, high cooperativity
end

"""
Enhanced competitive inhibition parameters based on biological context.
"""
struct CompetitiveInhibitionParams
    beta::Float64      # Inhibition strength
    m::Float64         # Hill coefficient for inhibition
    ki::Float64        # Inhibition constant
    competitive_factor::Float64  # Competitive vs non-competitive ratio
end

"""
Pathway-specific parameter sets based on biological literature.
"""
struct PathwayTypeParams
    # Basic Hill parameters
    h_range::Tuple{Float64, Float64}
    K_range::Tuple{Float64, Float64}
    
    # Inhibition parameters
    competitive_inhibition::CompetitiveInhibitionParams
    
    # Time constants
    response_time::Float64
    decay_time::Float64
    
    # Cooperativity characteristics
    typical_cooperativity::Float64
    noise_level::Float64
end

"""
Create pathway-specific parameters based on biological knowledge.
"""
function get_pathway_params(pathway_type::PathwayType)::PathwayTypeParams
    
    if pathway_type == SIGNALING
        return PathwayTypeParams(
            (1.8, 3.5),  # h_range: moderate to high cooperativity
            (0.15, 0.4), # K_range: sensitive thresholds
            CompetitiveInhibitionParams(12.0, 4.0, 0.05, 0.8),  # Strong competitive inhibition
            0.5,   # Fast response (minutes)
            2.0,   # Moderate decay
            2.5,   # Typical cooperativity
            0.15   # Moderate noise
        )
        
    elseif pathway_type == METABOLIC
        return PathwayTypeParams(
            (1.2, 2.2),  # h_range: enzyme-like kinetics
            (0.2, 0.6),  # K_range: Km-like values
            CompetitiveInhibitionParams(15.0, 3.0, 0.1, 0.9),  # Very strong, highly competitive
            5.0,   # Slower response (minutes to hours)
            10.0,  # Slow decay - steady state
            1.8,   # Lower cooperativity
            0.1    # Low noise - metabolic homeostasis
        )
        
    elseif pathway_type == TRANSCRIPTIONAL
        return PathwayTypeParams(
            (3.0, 6.0),  # h_range: very high cooperativity
            (0.1, 0.3),  # K_range: sharp thresholds
            CompetitiveInhibitionParams(20.0, 5.0, 0.02, 0.6), # Very strong, less competitive
            15.0,  # Slow response (hours)
            30.0,  # Very slow decay
            4.0,   # High cooperativity
            0.25   # Higher noise - transcriptional bursting
        )
        
    elseif pathway_type == STRESS_RESPONSE
        return PathwayTypeParams(
            (4.0, 8.0),  # h_range: ultra-high cooperativity
            (0.05, 0.15), # K_range: very sharp thresholds
            CompetitiveInhibitionParams(25.0, 6.0, 0.01, 0.7), # Extremely strong
            1.0,   # Very fast response (seconds to minutes)
            5.0,   # Fast decay when stress removed
            6.0,   # Very high cooperativity
            0.3    # High noise - stress variability
        )
        
    elseif pathway_type == CELL_CYCLE
        return PathwayTypeParams(
            (2.5, 4.5),  # h_range: checkpoint-like
            (0.3, 0.5),  # K_range: robust thresholds
            CompetitiveInhibitionParams(18.0, 4.5, 0.03, 0.8), # Strong, ordered inhibition
            30.0,  # Slow, ordered progression
            60.0,  # Very slow - cell cycle length
            3.5,   # High cooperativity for checkpoints
            0.2    # Moderate noise
        )
        
    elseif pathway_type == APOPTOSIS
        return PathwayTypeParams(
            (5.0, 10.0), # h_range: irreversible switch
            (0.02, 0.1), # K_range: very sharp
            CompetitiveInhibitionParams(30.0, 8.0, 0.005, 0.5), # Extremely strong
            2.0,   # Fast commitment
            1000.0, # No decay - irreversible
            8.0,   # Ultra-high cooperativity
            0.4    # High noise - life/death decisions
        )
    end
end

"""
Enhanced competitive inhibition model with biological realism.
"""
function enhanced_competitive_inhibition(
    inhibitor_concentrations::Vector{Float64},
    substrate_concentration::Float64,
    params::CompetitiveInhibitionParams
)::Float64
    
    if isempty(inhibitor_concentrations)
        return 1.0
    end
    
    # Competitive inhibition: Ki / (Ki + I) for each inhibitor
    competitive_factors = Float64[]
    
    for I in inhibitor_concentrations
        I_safe = clamp(I, 0.0, 1.0)
        
        # Competitive inhibition with substrate competition
        if substrate_concentration > 0.0
            # Classic competitive inhibition: 1 / (1 + I/Ki * (1 + S/Km))
            # Simplified: stronger when substrate is present
            substrate_factor = 1.0 + substrate_concentration / 0.2  # Assume Km = 0.2
            effective_inhibition = params.beta * I_safe^params.m * substrate_factor
        else
            effective_inhibition = params.beta * I_safe^params.m
        end
        
        inhibition_factor = 1.0 / (1.0 + effective_inhibition)
        push!(competitive_factors, inhibition_factor)
    end
    
    # Multiple inhibitors can be competitive or non-competitive
    if length(competitive_factors) == 1
        return competitive_factors[1]
    else
        # Mixed inhibition model
        competitive_product = prod(competitive_factors)
        non_competitive_product = prod([f^(1-params.competitive_factor) for f in competitive_factors])
        
        return competitive_product * params.competitive_factor + 
               non_competitive_product * (1 - params.competitive_factor)
    end
end

"""
Create biologically realistic reaction parameters for a specific pathway type.
"""
function create_biological_reaction_params(
    pathway_type::PathwayType,
    n_activators::Int,
    n_inhibitors::Int,
    n_substrates::Int;
    biological_context::String = ""
)::ReactionParams
    
    params = get_pathway_params(pathway_type)
    
    # Sample Hill parameters from biological range
    h = params.h_range[1] + rand() * (params.h_range[2] - params.h_range[1])
    K = params.K_range[1] + rand() * (params.K_range[2] - params.K_range[1])
    
    # Activator parameters with pathway-appropriate characteristics
    activator_weights = n_activators > 0 ? fill(1.0/n_activators, n_activators) : Float64[]
    
    # Pathway-specific sensitivity
    if pathway_type == SIGNALING
        # Signaling pathways often have switch-like behavior
        activator_sensitivity_s = fill(0.8 + rand() * 1.2, n_activators)  # 0.8-2.0
    elseif pathway_type == METABOLIC
        # Metabolic pathways more linear
        activator_sensitivity_s = fill(0.0 + rand() * 0.5, n_activators)  # 0.0-0.5
    else
        # Other pathways moderate sensitivity
        activator_sensitivity_s = fill(0.3 + rand() * 0.7, n_activators)  # 0.3-1.0
    end
    
    activator_sensitivity_n = fill(params.typical_cooperativity + randn() * 0.3, n_activators)
    activator_sensitivity_K = fill(0.1 + rand() * 0.2, n_activators)
    
    # Enhanced competitive inhibition parameters
    inhibitor_betas = fill(params.competitive_inhibition.beta + randn() * 2.0, n_inhibitors)
    inhibitor_ms = fill(params.competitive_inhibition.m + randn() * 0.5, n_inhibitors)
    
    # Ensure minimum strength for competitive inhibition
    inhibitor_betas = max.(inhibitor_betas, 10.0)
    inhibitor_ms = max.(inhibitor_ms, 2.5)
    
    # Substrate parameters
    substrate_weights = n_substrates > 0 ? fill(1.0/n_substrates, n_substrates) : Float64[]
    
    # Time-dynamic parameters based on pathway type
    consumption_lambdas = fill(0.1 / params.response_time, n_substrates)
    production_etas = Float64[]  # Will be filled based on products
    replenishment_rho = 0.01 / params.response_time
    decay_delta = 1.0 / params.decay_time
    
    return ReactionParams(
        h, K,
        activator_weights, activator_sensitivity_s, activator_sensitivity_n, activator_sensitivity_K,
        inhibitor_betas, inhibitor_ms,
        substrate_weights,
        consumption_lambdas, production_etas,
        replenishment_rho, decay_delta
    )
end

"""
Infer pathway type from network topology and node names.
"""
function infer_pathway_type(target_node::NetworkNode, edges::Vector{LogicNetworkEdge})::PathwayType
    
    # Use node names and topology to infer pathway type
    node_name = lowercase(get(target_node.name, target_node.uuid, ""))
    
    # Transcriptional markers
    if contains(node_name, "transcription") || contains(node_name, "tf") || 
       contains(node_name, "promoter") || contains(node_name, "gene")
        return TRANSCRIPTIONAL
    end
    
    # Metabolic markers
    if contains(node_name, "enzyme") || contains(node_name, "metabolite") ||
       contains(node_name, "kinase") && contains(node_name, "metabolic")
        return METABOLIC
    end
    
    # Stress response markers
    if contains(node_name, "stress") || contains(node_name, "heat") ||
       contains(node_name, "dna_damage") || contains(node_name, "p53")
        return STRESS_RESPONSE
    end
    
    # Cell cycle markers
    if contains(node_name, "cyclin") || contains(node_name, "cdk") ||
       contains(node_name, "checkpoint") || contains(node_name, "cell_cycle")
        return CELL_CYCLE
    end
    
    # Apoptosis markers
    if contains(node_name, "apoptosis") || contains(node_name, "caspase") ||
       contains(node_name, "death") || contains(node_name, "bax")
        return APOPTOSIS
    end
    
    # Topology-based inference
    n_inhibitors = count(e -> !e.is_positive, edges)
    n_activators = count(e -> e.is_positive, edges)
    
    # High inhibition ratio suggests regulatory pathway
    if n_inhibitors > n_activators && n_inhibitors > 2
        return TRANSCRIPTIONAL
    end
    
    # Many inputs suggests signaling hub
    if length(edges) > 5
        return SIGNALING
    end
    
    # Default to signaling
    return SIGNALING
end

"""
Create enhanced reaction network with biological pathway-specific parameters.
"""
function create_biologically_enhanced_reaction_network(network::ReactionNetwork)::Vector{Reaction}
    
    println("🧬 Creating biologically enhanced reaction network...")
    
    # Group edges by target
    target_groups = Dict{String, Vector{LogicNetworkEdge}}()
    for edge in network.edges
        target = edge.child_uuid
        if !haskey(target_groups, target)
            target_groups[target] = LogicNetworkEdge[]
        end
        push!(target_groups[target], edge)
    end
    
    reactions = Reaction[]
    pathway_type_counts = Dict{PathwayType, Int}()
    
    for (target_uuid, edges) in target_groups
        # Get target node
        target_node = get(network.nodes, target_uuid, 
                         NetworkNode(target_uuid, nothing, nothing, nothing, nothing))
        
        # Infer pathway type
        pathway_type = infer_pathway_type(target_node, edges)
        pathway_type_counts[pathway_type] = get(pathway_type_counts, pathway_type, 0) + 1
        
        # Separate edge types
        activators = [e.parent_uuid for e in edges if e.is_positive]
        inhibitors = [e.parent_uuid for e in edges if !e.is_positive]
        substrates = String[]  # Could be inferred from stoichiometry > 1
        products = String[]
        
        # Create biologically realistic parameters
        params = create_biological_reaction_params(
            pathway_type,
            length(activators),
            length(inhibitors),
            length(substrates)
        )
        
        # Determine gate type
        is_and_gate = any(edge.is_and for edge in edges)
        
        reaction = Reaction(
            target_uuid,
            activators,
            inhibitors,
            substrates,
            products,
            params,
            is_and_gate
        )
        
        push!(reactions, reaction)
    end
    
    # Print pathway type distribution
    println("📊 Pathway type distribution:")
    for (ptype, count) in pathway_type_counts
        println("  $(ptype): $(count) reactions")
    end
    
    return reactions
end

"""
Enhanced reaction computation with competitive inhibition.
"""
function compute_enhanced_reaction_output(
    reaction::Reaction,
    node_activities::Dict{String, Float64}
)::Float64
    
    # Get activator inputs with sensitivity transforms
    activator_inputs = Float64[]
    for (i, uuid) in enumerate(reaction.activator_uuids)
        x = get(node_activities, uuid, 0.0)
        if i <= length(reaction.params.activator_sensitivity_s)
            x_transformed = apply_sensitivity_transform(
                x;
                s=reaction.params.activator_sensitivity_s[i],
                n=reaction.params.activator_sensitivity_n[i],
                K_α=reaction.params.activator_sensitivity_K[i]
            )
            push!(activator_inputs, x_transformed)
        else
            push!(activator_inputs, x)
        end
    end
    
    # Get inhibitor inputs
    inhibitor_inputs = Float64[]
    for uuid in reaction.inhibitor_uuids
        x = get(node_activities, uuid, 0.0)
        push!(inhibitor_inputs, x)
    end
    
    # Get substrate inputs
    substrate_inputs = Float64[]
    for uuid in reaction.substrate_uuids
        x = get(node_activities, uuid, 0.0)
        push!(substrate_inputs, x)
    end
    
    # Compute activator aggregation
    A = if isempty(activator_inputs)
        reaction.is_and_gate ? 1.0 : 0.0
    else
        apply_logic_aggregation(activator_inputs, reaction.is_and_gate, reaction.params.activator_weights)
    end
    
    # Enhanced competitive inhibition
    H = if isempty(inhibitor_inputs)
        1.0
    else
        # Use first substrate as representative substrate concentration
        substrate_conc = isempty(substrate_inputs) ? 0.0 : substrate_inputs[1]
        
        # Create competitive inhibition params from reaction params
        if !isempty(reaction.params.inhibitor_betas) && !isempty(reaction.params.inhibitor_ms)
            comp_params = CompetitiveInhibitionParams(
                reaction.params.inhibitor_betas[1],  # Use first as representative
                reaction.params.inhibitor_ms[1],
                0.05,  # Default Ki
                0.8    # Default competitive factor
            )
            enhanced_competitive_inhibition(inhibitor_inputs, substrate_conc, comp_params)
        else
            # Fallback to standard inhibition
            inhibition_aggregator(inhibitor_inputs, reaction.params.inhibitor_betas, reaction.params.inhibitor_ms)
        end
    end
    
    # Substrate availability
    L = if isempty(substrate_inputs)
        1.0
    else
        substrate_availability_aggregator(substrate_inputs, reaction.params.substrate_weights)
    end
    
    # Pre-activation signal
    s = A * H * L
    
    # Final Hill activation
    h = reaction.params.h
    K = reaction.params.K
    
    output = s^h / (s^h + K^h)
    
    return clamp(output, 0.0, 1.0)
end

"""
Export pathway-specific parameter statistics.
"""
function analyze_pathway_parameters(reactions::Vector{Reaction})::Dict{String, Any}
    
    analysis = Dict{String, Any}()
    
    # Collect parameter distributions
    hill_coeffs = [r.params.h for r in reactions]
    thresholds = [r.params.K for r in reactions]
    
    # Inhibition strength analysis
    all_betas = Float64[]
    all_ms = Float64[]
    
    for reaction in reactions
        append!(all_betas, reaction.params.inhibitor_betas)
        append!(all_ms, reaction.params.inhibitor_ms)
    end
    
    analysis["hill_coefficient"] = Dict(
        "mean" => mean(hill_coeffs),
        "std" => std(hill_coeffs),
        "range" => (minimum(hill_coeffs), maximum(hill_coeffs))
    )
    
    analysis["threshold"] = Dict(
        "mean" => mean(thresholds),
        "std" => std(thresholds),
        "range" => (minimum(thresholds), maximum(thresholds))
    )
    
    if !isempty(all_betas)
        analysis["inhibition_strength"] = Dict(
            "mean_beta" => mean(all_betas),
            "std_beta" => std(all_betas),
            "mean_m" => mean(all_ms),
            "std_m" => std(all_ms),
            "range_beta" => (minimum(all_betas), maximum(all_betas))
        )
    end
    
    analysis["reaction_count"] = length(reactions)
    analysis["total_inhibitors"] = sum(length(r.inhibitor_uuids) for r in reactions)
    analysis["total_activators"] = sum(length(r.activator_uuids) for r in reactions)
    
    return analysis
end