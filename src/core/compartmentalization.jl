# Cellular Compartmentalization Module
# Models how subcellular localization affects reaction kinetics

"""
Cellular compartments with different biochemical environments.
"""
@enum Compartment begin
    CYTOPLASM
    NUCLEUS
    MITOCHONDRIA
    ENDOPLASMIC_RETICULUM
    GOLGI
    LYSOSOME
    PEROXISOME
    MEMBRANE
    EXTRACELLULAR
end

"""
Compartment-specific environment parameters affecting reaction kinetics.
"""
struct CompartmentEnvironment
    ph::Float64                    # pH affects ionization states
    ionic_strength::Float64        # Affects protein stability
    volume_fraction::Float64       # Relative volume in cell
    protein_crowding::Float64      # Macromolecular crowding effects
    atp_concentration::Float64     # Energy availability
    ca_concentration::Float64      # Ca2+ signaling
    redox_potential::Float64       # Oxidizing/reducing environment
    temperature_factor::Float64    # Relative temperature effects
end

"""
Transport parameters between compartments.
"""
struct TransportParams
    passive_permeability::Float64     # Passive diffusion rate
    active_transport_rate::Float64    # Active transport Vmax
    transport_km::Float64              # Michaelis constant for transport
    energy_requirement::Float64       # ATP cost per molecule
    selectivity_factor::Float64       # Selectivity vs other molecules
end

"""
Compartmentalized reaction with location-specific parameters.
"""
struct CompartmentalizedReaction
    base_reaction::Reaction
    compartment::Compartment
    environment::CompartmentEnvironment
    transport_map::Dict{String, TransportParams}  # node_id -> transport params
end

"""
Get standard compartment environments based on cell biology literature.
"""
function get_compartment_environment(compartment::Compartment)::CompartmentEnvironment
    
    if compartment == CYTOPLASM
        return CompartmentEnvironment(
            7.2,    # pH - slightly basic
            0.15,   # ionic strength (M)
            0.65,   # ~65% of cell volume
            0.3,    # 30% macromolecular crowding
            5e-3,   # 5 mM ATP
            1e-7,   # 100 nM free Ca2+
            -0.25,  # reducing environment
            1.0     # baseline temperature
        )
        
    elseif compartment == NUCLEUS
        return CompartmentEnvironment(
            7.3,    # Slightly more basic than cytoplasm
            0.12,   # Lower ionic strength
            0.10,   # ~10% of cell volume
            0.4,    # High protein crowding (chromatin, etc.)
            3e-3,   # Lower ATP than cytoplasm
            1e-7,   # Similar Ca2+
            -0.2,   # Slightly more oxidizing
            1.0
        )
        
    elseif compartment == MITOCHONDRIA
        return CompartmentEnvironment(
            7.8,    # More basic (matrix)
            0.2,    # High ionic strength
            0.15,   # ~15% of cell volume
            0.5,    # Very high protein density
            8e-3,   # High ATP synthesis
            1e-6,   # Higher Ca2+ (signaling hub)
            -0.3,   # Very reducing (NADH/FADH2)
            1.1     # Slightly higher temperature from metabolism
        )
        
    elseif compartment == ENDOPLASMIC_RETICULUM
        return CompartmentEnvironment(
            7.1,    # Slightly acidic
            0.1,    # Lower ionic strength
            0.05,   # ~5% of cell volume
            0.3,    # Moderate crowding
            2e-3,   # Low ATP
            5e-4,   # High Ca2+ storage
            0.1,    # Oxidizing (disulfide bond formation)
            1.0
        )
        
    elseif compartment == GOLGI
        return CompartmentEnvironment(
            6.7,    # More acidic
            0.1,    
            0.02,   # Small volume
            0.4,    # High enzyme density
            1e-3,   # Low ATP
            1e-6,   
            0.2,    # Oxidizing
            1.0
        )
        
    elseif compartment == LYSOSOME
        return CompartmentEnvironment(
            4.5,    # Very acidic
            0.05,   
            0.01,   # Small volume
            0.6,    # Very high enzyme density
            5e-4,   # Very low ATP
            1e-5,   # Variable Ca2+
            0.3,    # Oxidizing
            1.0
        )
        
    elseif compartment == MEMBRANE
        return CompartmentEnvironment(
            7.0,    # Interface effects
            0.2,    # High local ionic strength
            0.05,   # Membrane area
            0.7,    # Very high protein density
            3e-3,   
            1e-6,   
            0.0,    # Mixed redox
            1.0
        )
        
    elseif compartment == EXTRACELLULAR
        return CompartmentEnvironment(
            7.4,    # Physiological pH
            0.15,   # Physiological ionic strength
            1000.0, # "Infinite" dilution
            0.05,   # Low crowding
            1e-6,   # Very low ATP
            2e-3,   # mM Ca2+ in blood
            0.1,    # Oxidizing
            1.0
        )
        
    else
        # Default to cytoplasm
        return get_compartment_environment(CYTOPLASM)
    end
end

"""
Adjust reaction parameters based on compartment environment.
"""
function adjust_for_compartment(
    params::ReactionParams, 
    environment::CompartmentEnvironment
)::ReactionParams
    
    # pH effects on Hill coefficient (ionization changes cooperativity)
    ph_factor = exp(-abs(environment.ph - 7.2) * 0.2)  # Optimal at pH 7.2
    adjusted_h = params.h * ph_factor
    
    # Ionic strength affects binding affinity (Debye-Huckel effects)
    ionic_factor = 1.0 + environment.ionic_strength * 2.0
    adjusted_K = params.K * ionic_factor
    
    # Protein crowding enhances effective concentrations
    crowding_factor = 1.0 + environment.protein_crowding * 3.0
    adjusted_betas = params.inhibitor_betas .* crowding_factor
    
    # ATP availability affects energy-dependent processes
    atp_factor = environment.atp_concentration / 5e-3  # Relative to cytoplasm
    adjusted_decay = params.decay_delta * atp_factor
    
    # Redox potential affects protein stability
    redox_factor = exp(environment.redox_potential * 0.5)
    adjusted_rho = params.replenishment_rho * redox_factor
    
    # Temperature effects on all rates
    temp_factor = environment.temperature_factor
    adjusted_lambdas = params.consumption_lambdas .* temp_factor
    
    return ReactionParams(
        adjusted_h,
        adjusted_K,
        params.activator_weights,
        params.activator_sensitivity_s,
        params.activator_sensitivity_n,
        params.activator_sensitivity_K,
        adjusted_betas,
        params.inhibitor_ms,
        params.substrate_weights,
        adjusted_lambdas,
        params.production_etas,
        adjusted_rho,
        adjusted_decay
    )
end

"""
Reactome compartment ID patterns for proper compartmentalization.
"""
const REACTOME_COMPARTMENT_PATTERNS = Dict{Compartment, Vector{String}}(
    CYTOPLASM => ["R-HSA-", "R-MMU-", "R-CEL-", "R-DME-", "R-SCE-", "R-SPO-", "R-ATH-"],  # Default
    NUCLEUS => ["R-HSA-.*[Nn]ucleus", "R-HSA-.*[Nn]uclear", "R-HSA-.*[Cc]hromatin"],
    MITOCHONDRIA => ["R-HSA-.*[Mm]itochondria", "R-HSA-.*[Mm]itochondrial"],
    ENDOPLASMIC_RETICULUM => ["R-HSA-.*[Ee][Rr]", "R-HSA-.*[Ee]ndoplasmic"],
    GOLGI => ["R-HSA-.*[Gg]olgi"],
    LYSOSOME => ["R-HSA-.*[Ll]ysosom", "R-HSA-.*[Vv]acuol"],
    PEROXISOME => ["R-HSA-.*[Pp]eroxisom"],
    MEMBRANE => ["R-HSA-.*[Mm]embrane", "R-HSA-.*[Pp]lasma"],
    EXTRACELLULAR => ["R-HSA-.*[Ee]xtracellular", "R-HSA-.*[Ss]ecret"]
)

"""
Infer compartment from Reactome ID and node information.
The same protein will have different Reactome IDs in different compartments.
"""
function infer_node_compartment(node::NetworkNode)::Compartment
    
    reactome_id = get(node.reactome_id, "", "")
    
    if isempty(reactome_id)
        # Fallback to name-based inference if no Reactome ID
        return infer_compartment_from_name(node)
    end
    
    # Check Reactome ID patterns for compartment-specific identifiers
    for (compartment, patterns) in REACTOME_COMPARTMENT_PATTERNS
        for pattern in patterns
            if occursin(Regex(pattern), reactome_id)
                return compartment
            end
        end
    end
    
    # For generic Reactome IDs without compartment info, use name/type inference
    return infer_compartment_from_name(node)
end

"""
Fallback compartment inference from node name and entity type.
Used when Reactome ID doesn't contain compartment information.
"""
function infer_compartment_from_name(node::NetworkNode)::Compartment
    
    name_lower = lowercase(get(node.name, node.uuid, ""))
    entity_type = get(node.entity_type, "unknown", "")
    
    # Explicit compartment markers in names
    if contains(name_lower, "nuclear") || contains(name_lower, "nucleus")
        return NUCLEUS
    elseif contains(name_lower, "mitochondrial") || contains(name_lower, "mito")
        return MITOCHONDRIA
    elseif contains(name_lower, "er_") || contains(name_lower, "endoplasmic")
        return ENDOPLASMIC_RETICULUM
    elseif contains(name_lower, "golgi")
        return GOLGI
    elseif contains(name_lower, "lysosomal") || contains(name_lower, "lysosome")
        return LYSOSOME
    elseif contains(name_lower, "membrane") || contains(name_lower, "receptor")
        return MEMBRANE
    elseif contains(name_lower, "extracellular") || contains(name_lower, "secreted")
        return EXTRACELLULAR
    end
    
    # Entity type based inference
    if entity_type == "transcription_factor"
        return NUCLEUS
    elseif entity_type == "receptor"
        return MEMBRANE
    elseif contains(entity_type, "extracellular") || contains(entity_type, "hormone")
        return EXTRACELLULAR
    end
    
    # Default to cytoplasm for generic proteins
    return CYTOPLASM
end

"""
Model transport between compartments for time-dynamic simulations.
"""
function compute_transport_flux(
    source_concentration::Float64,
    target_concentration::Float64,
    transport_params::TransportParams,
    atp_level::Float64 = 1.0
)::Float64
    
    # Passive diffusion component
    passive_flux = transport_params.passive_permeability * 
                   (source_concentration - target_concentration)
    
    # Active transport component (if energy available)
    if atp_level > 0.1  # Minimum ATP threshold
        # Michaelis-Menten kinetics for active transport
        active_flux = transport_params.active_transport_rate * 
                      source_concentration / 
                      (source_concentration + transport_params.transport_km)
        
        # Energy cost
        energy_factor = min(1.0, atp_level / transport_params.energy_requirement)
        active_flux *= energy_factor
    else
        active_flux = 0.0
    end
    
    # Total flux with selectivity
    total_flux = (passive_flux + active_flux) * transport_params.selectivity_factor
    
    return total_flux
end

"""
Create compartmentalized reaction network based on Reactome compartment structure.
Each entity with a different Reactome ID in a different compartment gets 
appropriate environmental parameters.
"""
function create_compartmentalized_reactions(network::ReactionNetwork)::Vector{CompartmentalizedReaction}
    
    println("🏠 Creating compartmentalized reaction network...")
    println("Using Reactome compartment structure:")
    println("  - Same entity → Different Reactome IDs per compartment")
    println("  - Logic network UUIDs → Already compartment-specific")
    
    # First create standard reactions
    base_reactions = create_biologically_enhanced_reaction_network(network)
    
    compartmentalized_reactions = CompartmentalizedReaction[]
    compartment_counts = Dict{Compartment, Int}()
    cross_compartment_reactions = 0
    
    for reaction in base_reactions
        # Infer compartment for target node based on Reactome ID
        target_node = get(network.nodes, reaction.target_uuid, 
                         NetworkNode(reaction.target_uuid, nothing, nothing, nothing, nothing))
        
        target_compartment = infer_node_compartment(target_node)
        compartment_counts[target_compartment] = get(compartment_counts, target_compartment, 0) + 1
        
        # Check if this reaction involves cross-compartment transport
        input_compartments = Set{Compartment}()
        transport_map = Dict{String, TransportParams}()
        
        all_input_uuids = [reaction.activator_uuids; reaction.inhibitor_uuids; reaction.substrate_uuids]
        
        for input_uuid in all_input_uuids
            input_node = get(network.nodes, input_uuid,
                           NetworkNode(input_uuid, nothing, nothing, nothing, nothing))
            input_compartment = infer_node_compartment(input_node)
            push!(input_compartments, input_compartment)
            
            # If input is from different compartment, create transport parameters
            if input_compartment != target_compartment
                cross_compartment_reactions += 1
                transport_params = create_transport_params(input_compartment, target_compartment)
                transport_map[input_uuid] = transport_params
            else
                # Same compartment - no transport needed
                transport_map[input_uuid] = TransportParams(1.0, 0.0, 1.0, 0.0, 1.0)
            end
        end
        
        # Get target compartment environment
        environment = get_compartment_environment(target_compartment)
        
        # Adjust reaction parameters for target compartment
        adjusted_params = adjust_for_compartment(reaction.params, environment)
        
        adjusted_reaction = Reaction(
            reaction.target_uuid,
            reaction.activator_uuids,
            reaction.inhibitor_uuids,
            reaction.substrate_uuids,
            reaction.product_uuids,
            adjusted_params,
            reaction.is_and_gate
        )
        
        comp_reaction = CompartmentalizedReaction(
            adjusted_reaction,
            target_compartment,
            environment,
            transport_map
        )
        
        push!(compartmentalized_reactions, comp_reaction)
    end
    
    # Print compartment distribution
    println("🏠 Compartment distribution:")
    for (comp, count) in compartment_counts
        println("  $(comp): $(count) reactions")
    end
    
    if cross_compartment_reactions > 0
        println("🚚 Cross-compartment transport: $(cross_compartment_reactions) reactions")
    end
    
    return compartmentalized_reactions
end

"""
Create transport parameters based on source and target compartments.
"""
function create_transport_params(
    source_compartment::Compartment,
    target_compartment::Compartment
)::TransportParams
    
    # Define transport characteristics between compartments
    if source_compartment == EXTRACELLULAR && target_compartment == MEMBRANE
        # Receptor binding - high selectivity, no energy cost
        return TransportParams(0.8, 0.0, 0.1, 0.0, 0.9)
        
    elseif source_compartment == MEMBRANE && target_compartment == CYTOPLASM
        # Membrane transport - may require energy
        return TransportParams(0.3, 2.0, 0.5, 0.2, 0.8)
        
    elseif source_compartment == CYTOPLASM && target_compartment == NUCLEUS
        # Nuclear import - active transport, energy-dependent
        return TransportParams(0.1, 3.0, 0.3, 0.5, 0.9)
        
    elseif source_compartment == NUCLEUS && target_compartment == CYTOPLASM
        # Nuclear export - active transport
        return TransportParams(0.1, 2.0, 0.4, 0.3, 0.8)
        
    elseif source_compartment == CYTOPLASM && target_compartment == MITOCHONDRIA
        # Mitochondrial import - highly regulated
        return TransportParams(0.05, 1.5, 0.2, 0.8, 0.95)
        
    elseif source_compartment == CYTOPLASM && target_compartment == ENDOPLASMIC_RETICULUM
        # ER targeting - signal sequence dependent
        return TransportParams(0.2, 1.0, 0.3, 0.1, 0.9)
        
    else
        # Generic transport between other compartments
        return TransportParams(0.2, 1.0, 0.5, 0.2, 0.7)
    end
end

"""
Compute reaction output with compartmentalization effects.
"""
function compute_compartmentalized_reaction_output(
    comp_reaction::CompartmentalizedReaction,
    node_activities::Dict{String, Float64},
    compartment_activities::Dict{String, Dict{Compartment, Float64}} = Dict()
)::Float64
    
    # If compartment-specific activities are available, use them
    if !isempty(compartment_activities)
        # Get activities in the reaction's compartment
        local_activities = Dict{String, Float64}()
        
        for node_id in [comp_reaction.base_reaction.activator_uuids; 
                       comp_reaction.base_reaction.inhibitor_uuids; 
                       comp_reaction.base_reaction.substrate_uuids]
            
            if haskey(compartment_activities, node_id)
                # Use compartment-specific concentration
                local_conc = get(compartment_activities[node_id], comp_reaction.compartment, 0.0)
                local_activities[node_id] = local_conc
            else
                # Fallback to global activity
                local_activities[node_id] = get(node_activities, node_id, 0.0)
            end
        end
        
        # Use local activities
        return compute_enhanced_reaction_output(comp_reaction.base_reaction, local_activities)
    else
        # Use standard computation with adjusted parameters
        return compute_enhanced_reaction_output(comp_reaction.base_reaction, node_activities)
    end
end

"""
Analyze compartment distribution and effects.
"""
function analyze_compartment_effects(
    comp_reactions::Vector{CompartmentalizedReaction}
)::Dict{String, Any}
    
    analysis = Dict{String, Any}()
    
    # Count reactions per compartment
    compartment_counts = Dict{Compartment, Int}()
    for comp_reaction in comp_reactions
        compartment = comp_reaction.compartment
        compartment_counts[compartment] = get(compartment_counts, compartment, 0) + 1
    end
    
    analysis["compartment_distribution"] = compartment_counts
    
    # Analyze parameter adjustments by compartment
    compartment_params = Dict{Compartment, Vector{Float64}}()
    
    for comp_reaction in comp_reactions
        compartment = comp_reaction.compartment
        if !haskey(compartment_params, compartment)
            compartment_params[compartment] = Float64[]
        end
        
        # Collect Hill coefficients as example
        push!(compartment_params[compartment], comp_reaction.base_reaction.params.h)
    end
    
    # Statistical analysis per compartment
    compartment_stats = Dict{Compartment, Dict{String, Float64}}()
    
    for (compartment, hill_coeffs) in compartment_params
        if !isempty(hill_coeffs)
            compartment_stats[compartment] = Dict(
                "mean_hill_coefficient" => mean(hill_coeffs),
                "std_hill_coefficient" => std(hill_coeffs),
                "reaction_count" => length(hill_coeffs)
            )
        end
    end
    
    analysis["compartment_statistics"] = compartment_stats
    
    # Environment parameter summary
    environments = Dict{Compartment, CompartmentEnvironment}()
    for compartment in keys(compartment_counts)
        environments[compartment] = get_compartment_environment(compartment)
    end
    
    analysis["compartment_environments"] = environments
    
    return analysis
end