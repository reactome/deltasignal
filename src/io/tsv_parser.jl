# TSV Parser for Logic Networks with Set Expansion Handling

struct NetworkNode
    uuid::String
    reactome_id::Union{String, Nothing}
    entity_type::String
    original_set_id::Union{String, Nothing}
    display_name::String
    baseline::Float64
end

struct LogicNetworkEdge
    parent_uuid::String
    child_uuid::String
    is_and::Bool  # true for AND (1), false for OR (0)
    is_positive::Bool  # true for positive (1), false for negative (-1)
    stoichiometry::Float64
    edge_type::String  # generator: input/output/catalyst/regulator/assembly/
                       # dissociation/depletion. "" for sample format.
end

struct SetExpansionMapping
    original_set_id::String
    original_name::String
    expanded_members::Vector{String}
    reactome_pathway_coords::Union{Dict{String, Any}, Nothing}
end

struct ReactionNetwork
    nodes::Dict{String, NetworkNode}
    edges::Vector{LogicNetworkEdge}
    set_mappings::Dict{String, SetExpansionMapping}
end

"""
Parse the main logic network file. Accepts two schemas:

1. **Sample/bundled format** (TSV): `parent | child | is_and | is_positive | stoichiometry`
   with 1/0 integer encoding for is_and and is_positive.

2. **Generator format** (CSV from logic-network-generator): `source_id, target_id,
   pos_neg, and_or, edge_type, stoichiometry` with string encoding ("pos"/"neg",
   "and"/"or"/""). The extra `edge_type` column is ignored.

Schema is detected from the header. Delimiter is autodetected by CSV.jl.
"""
function parse_logic_network(filepath::String)::Vector{LogicNetworkEdge}
    df = CSV.read(filepath, DataFrame, header=true)
    cols = Set(Symbol.(names(df)))

    # Generator schema detection: source_id + target_id + pos_neg
    if :source_id in cols && :target_id in cols && :pos_neg in cols
        return parse_logic_network_generator(df)
    end

    # Sample/bundled schema
    if ncol(df) != 5
        throw(ArgumentError("Logic network must have exactly 5 columns (parent, child, is_and, is_positive, stoichiometry) for the sample format, or generator-format columns (source_id, target_id, pos_neg, and_or, ...)."))
    end

    edges = LogicNetworkEdge[]
    for row in eachrow(df)
        is_and = row.is_and == 1
        is_positive = row.is_positive == 1
        stoich = Float64(row.stoichiometry)
        push!(edges, LogicNetworkEdge(
            String(row.parent),
            String(row.child),
            is_and,
            is_positive,
            stoich,
            "",  # sample format has no edge_type
        ))
    end
    return edges
end

"""
Generator-format edge parser. and_or is "and" / "or" / missing (treat empty as OR).
pos_neg is "pos" / "neg".
"""
function parse_logic_network_generator(df::DataFrame)::Vector{LogicNetworkEdge}
    edges = LogicNetworkEdge[]
    has_edge_type = "edge_type" in names(df)
    for row in eachrow(df)
        and_raw = ismissing(row.and_or) ? "" : lowercase(strip(String(row.and_or)))
        pos_raw = lowercase(strip(String(row.pos_neg)))

        is_and = and_raw == "and"
        is_positive = pos_raw == "pos"
        stoich = ismissing(row.stoichiometry) ? 1.0 : Float64(row.stoichiometry)
        et = if has_edge_type && !ismissing(row.edge_type)
            lowercase(strip(String(row.edge_type)))
        else
            ""
        end

        push!(edges, LogicNetworkEdge(
            String(row.source_id),
            String(row.target_id),
            is_and,
            is_positive,
            stoich,
            et,
        ))
    end
    return edges
end

"""
Parse UUID-to-Reactome mapping file. Accepts two schemas:

1. **Sample/bundled format**: `uuid, name, reactome_id, entity_type[, set_id]` (4-5 cols).
2. **Generator format** (`stid_to_uuid_mapping.csv`): `uuid, stable_id` only — name and
   entity_type are absent. Display name falls back to the stable_id; entity_type is
   marked `"unknown"`. Real names should be enriched downstream (Reactome ContentService).
"""
function parse_uuid_mapping(filepath::String)::Dict{String, NetworkNode}
    df = CSV.read(filepath, DataFrame, header=true)
    cols = Set(Symbol.(names(df)))

    # Generator schema detection: uuid + stable_id only
    if :uuid in cols && :stable_id in cols && !(:name in cols)
        return parse_uuid_mapping_generator(df)
    end

    if ncol(df) < 4
        throw(ArgumentError("UUID mapping must have at least 4 columns (uuid, name, reactome_id, entity_type) for the sample format, or generator-format columns (uuid, stable_id)."))
    end

    if ncol(df) == 4
        df.set_id = [missing for _ in 1:nrow(df)]
    end

    nodes = Dict{String, NetworkNode}()
    for row in eachrow(df)
        display_name = row.name
        set_id = haskey(row, :set_id) && !ismissing(row.set_id) ? String(row.set_id) : nothing
        reactome_id = ismissing(row.reactome_id) ? nothing : String(row.reactome_id)
        nodes[String(row.uuid)] = NetworkNode(
            String(row.uuid),
            reactome_id,
            String(row.entity_type),
            set_id,
            display_name,
            0.01,  # Spec x₀ = 0.01 (= 1 on UI scale, the "normal" baseline).
        )
    end
    return nodes
end

"""
Generator-format uuid mapping parser. Only `uuid` and `stable_id` are present;
synthesize missing fields so downstream code can work uniformly. The generator
writes the literal Python string "None" (or "nan") when a stable_id is absent —
treat those as missing.
"""
function parse_uuid_mapping_generator(df::DataFrame)::Dict{String, NetworkNode}
    nodes = Dict{String, NetworkNode}()
    for row in eachrow(df)
        uuid = String(row.uuid)
        stable_raw = ismissing(row.stable_id) ? nothing : String(row.stable_id)
        stable_id = if stable_raw === nothing || stable_raw in ("None", "nan", "")
            nothing
        else
            stable_raw
        end
        nodes[uuid] = NetworkNode(
            uuid,
            stable_id,
            "unknown",
            nothing,
            stable_id === nothing ? uuid : stable_id,
            0.01,  # Spec x₀ = 0.01 (= 1 on UI scale, the "normal" baseline).
        )
    end
    return nodes
end

"""
Parse set mappings file (optional).
Expected columns: Set ID | Original Name | Member UUIDs (comma-separated)
"""
function parse_set_mappings(filepath::String)::Dict{String, SetExpansionMapping}
    if !isfile(filepath)
        return Dict{String, SetExpansionMapping}()
    end
    
    df = CSV.read(filepath, DataFrame, header=true)
    
    if ncol(df) < 4
        throw(ArgumentError("Set mappings TSV must have at least 4 columns: set_id, uuid, set_name, set_type"))
    end
    
    # Group by set_id to collect all members
    mappings = Dict{String, SetExpansionMapping}()
    
    for row in eachrow(df)
        set_id = String(row.set_id)
        uuid = String(row.uuid)
        set_name = String(row.set_name)
        
        if haskey(mappings, set_id)
            # Add to existing mapping
            push!(mappings[set_id].expanded_members, uuid)
        else
            # Create new mapping
            mapping = SetExpansionMapping(
                set_id,
                set_name,
                [uuid],
                nothing  # Pathway coordinates to be filled later if available
            )
            mappings[set_id] = mapping
        end
    end
    
    return mappings
end

"""
Create a complete ReactionNetwork from parsed components.
"""
function create_reaction_network(
    edges::Vector{LogicNetworkEdge},
    nodes::Dict{String, NetworkNode},
    set_mappings::Dict{String, SetExpansionMapping} = Dict{String, SetExpansionMapping}()
)::ReactionNetwork
    
    # Validate that all edge UUIDs exist in nodes
    all_uuids = Set(keys(nodes))
    
    for edge in edges
        if !(edge.parent_uuid in all_uuids)
            @warn "Parent UUID $(edge.parent_uuid) not found in node mappings"
        end
        if !(edge.child_uuid in all_uuids)
            @warn "Child UUID $(edge.child_uuid) not found in node mappings"
        end
    end
    
    return ReactionNetwork(nodes, edges, set_mappings)
end

"""
Main parsing function that combines all inputs.
"""
function parse_complete_network(
    logic_network_path::String,
    uuid_mapping_path::String,
    set_mapping_path::Union{String, Nothing} = nothing
)::ReactionNetwork
    
    println("Parsing logic network...")
    edges = parse_logic_network(logic_network_path)
    println("Found $(length(edges)) edges")
    
    println("Parsing UUID mappings...")
    nodes = parse_uuid_mapping(uuid_mapping_path)  
    println("Found $(length(nodes)) nodes")
    
    set_mappings = Dict{String, SetExpansionMapping}()
    if set_mapping_path !== nothing
        println("Parsing set mappings...")
        set_mappings = parse_set_mappings(set_mapping_path)
        println("Found $(length(set_mappings)) set expansions")
    end
    
    return create_reaction_network(edges, nodes, set_mappings)
end