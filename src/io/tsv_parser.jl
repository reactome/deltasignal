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
    # Cofactor stable ids that shipped WITH this network, from the generator's
    # `cofactors.csv`. Empty when the bundle predates that file, in which case
    # the solver falls back to its own built-in list. Carrying it here is what
    # lets an artifact bundle pulled from S3 answer "which of these nodes is
    # ATP" without a second, separately-versioned copy of the answer.
    cofactor_stids::Set{String}
    # stable id -> the stable ids it CONTAINS (itself excluded), from the
    # generator's `containment.csv`. Read only by DS_SELF_INHIBITOR_WEIGHT
    # (specs/022) to find inhibitors built from their own reaction's input.
    # Empty when the bundle has no such file, or when a network arrives as JSON
    # (the rule is then inert, which the solve reports).
    containment::Dict{String, Set{String}}

    ReactionNetwork(nodes, edges, set_mappings,
                    cofactor_stids = Set{String}(),
                    containment = Dict{String, Set{String}}()) =
        new(nodes, edges, set_mappings, cofactor_stids, containment)
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
        # Same guard for the sample format: anything but 0/1 is a mistake,
        # and `== 1` would quietly turn it into OR / negative.
        if !(row.is_and in (0, 1)) || !(row.is_positive in (0, 1))
            throw(ArgumentError(
                "logic network row $(row.parent) -> $(row.child): is_and and " *
                "is_positive must be 0 or 1, got $(row.is_and) and $(row.is_positive)."))
        end
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

        # An unrecognised value must FAIL, not fall through to a default.
        # `pos_raw == "pos"` silently made every unknown sign an INHIBITOR,
        # which inverts the edge and is invisible downstream — the same class
        # as the DS_* mode typos that selected a different model in silence.
        # Only pos/neg and and/or/"" occur in the catalog today, so this
        # changes nothing about current data; it makes a future typo loud.
        if pos_raw != "pos" && pos_raw != "neg"
            throw(ArgumentError(
                "logic network row $(row.source_id) -> $(row.target_id): " *
                "pos_neg must be \"pos\" or \"neg\", got \"$(pos_raw)\"."))
        end
        if and_raw != "and" && and_raw != "or" && and_raw != ""
            throw(ArgumentError(
                "logic network row $(row.source_id) -> $(row.target_id): " *
                "and_or must be \"and\", \"or\" or empty, got \"$(and_raw)\"."))
        end

        is_and = and_raw == "and"     # empty means OR, as documented
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
    set_mapping_path::Union{String, Nothing} = nothing,
    cofactor_path::Union{String, Nothing} = nothing
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
    
    network = create_reaction_network(edges, nodes, set_mappings)

    # Prefer the list that shipped with the networks; see parse_cofactor_list.
    resolved = cofactor_path === nothing ?
        default_cofactor_path(logic_network_path) : cofactor_path
    stids = parse_cofactor_list(resolved; required = cofactor_path !== nothing)
    containment = parse_containment(logic_network_path)

    if isempty(stids)
        # A file that declares nothing in-network is NOT the same as no file,
        # and the generator writes every known cofactor precisely so the two
        # are distinguishable. Falling straight through to the built-in list
        # collapses that distinction and skips the disagreement check below,
        # so the 0-of-N case — the shape a wholly mis-written file takes — was
        # the one case that could never warn.
        if resolved !== nothing
            builtin_here = count(node -> node.reactome_id !== nothing &&
                                 node.reactome_id in COFACTOR_STIDS,
                                 values(network.nodes))
            if builtin_here > 0
                @warn("The bundled cofactor list declares nothing present in " *
                      "this network, but the built-in list matches nodes here. " *
                      "Falling back to the built-in list.",
                      file = resolved, builtin_would_match = builtin_here)
            end
        end
        return ReactionNetwork(network.nodes, network.edges, network.set_mappings,
                               Set{String}(), containment)
    end
    println("Found $(length(stids)) cofactor species declared by the bundle")

    # A truncated, stale or partially-written cofactors.csv silently narrows the
    # model: the bundle wins, so declaring one cofactor where the network holds
    # seven just quietly stops treating the other six as cofactors. Nothing
    # else would ever notice, which is the silent-substitution failure the DS_*
    # guard rails exist to prevent. Compare against what the built-in list
    # would have matched and say so once, at load.
    declared_nodes = Set(uuid for (uuid, node) in network.nodes
                         if node.reactome_id !== nothing && node.reactome_id in stids)
    builtin_nodes = Set(uuid for (uuid, node) in network.nodes
                        if node.reactome_id !== nothing &&
                           node.reactome_id in COFACTOR_STIDS)
    missed = setdiff(builtin_nodes, declared_nodes)
    if !isempty(missed)
        examples = sort([string(network.nodes[u].reactome_id) for u in missed])
        @warn("The bundled cofactor list matches fewer nodes than the built-in " *
              "list would. The bundle is authoritative, so these are NOT being " *
              "treated as cofactors. Expected if the bundle predates a list " *
              "change; a truncated or partly-written file looks the same.",
              file = resolved,
              declared = length(declared_nodes),
              builtin_would_match = length(builtin_nodes),
              not_treated_as_cofactors = first(unique(examples), 5))
    end

    return ReactionNetwork(network.nodes, network.edges, network.set_mappings, stids,
                           containment)
end

"""
Read the generator's `containment.csv` beside a logic network: stable id ->
the stable ids it contains, itself excluded. Empty when there is no file.
"""
function parse_containment(logic_network_path::String)::Dict{String, Set{String}}
    out = Dict{String, Set{String}}()
    path = joinpath(dirname(logic_network_path), "containment.csv")
    isfile(path) || return out
    df = CSV.read(path, DataFrame; types = String)
    for col in ("stable_id", "contains_stable_id")
        col in names(df) || throw(ArgumentError("$path has no `$col` column"))
    end
    for row in eachrow(df)
        (ismissing(row.stable_id) || ismissing(row.contains_stable_id)) && continue
        row.stable_id == row.contains_stable_id && continue
        push!(get!(out, String(row.stable_id), Set{String}()), String(row.contains_stable_id))
    end
    return out
end

"""
Path to the `cofactors.csv` the generator writes beside a logic network, or
`nothing` when there is none. The artifacts travel together, so the list is
found the same way the network was.
"""
function default_cofactor_path(logic_network_path::String)::Union{String, Nothing}
    candidate = joinpath(dirname(logic_network_path), "cofactors.csv")
    return isfile(candidate) ? candidate : nothing
end

"""
Read the cofactor stable ids a generated bundle declares.

Returns an empty set when there is no file — bundles generated before the
generator emitted one are the common case, and the solver falls back to its
built-in list rather than silently treating the pathway as cofactor-free.

Only rows flagged `in_network` are returned: the file lists every cofactor the
release defines so that an empty intersection is distinguishable from a missing
file, but only the ones actually present here can match a node.
"""
function parse_cofactor_list(path::Union{String, Nothing};
                            required::Bool = false)::Set{String}
    path === nothing && return Set{String}()
    if !isfile(path)
        # Absent is the norm for DISCOVERY — most bundles predate the file.
        # But a caller who named a path expressed intent, and silently handing
        # them the built-in list instead is the wrong kind of forgiving.
        required && throw(ArgumentError("cofactor list not found: $(path)"))
        return Set{String}()
    end
    out = Set{String}()
    df = CSV.read(path, DataFrame)
    cols = Set(Symbol.(names(df)))
    (:stable_id in cols) || throw(ArgumentError(
        "$(path) has no `stable_id` column; it is not a generator cofactor list"))
    has_flag = :in_network in cols
    for row in eachrow(df)
        ismissing(row.stable_id) && continue
        if has_flag && !_in_network_flag(row.in_network, path)
            continue
        end
        push!(out, String(row.stable_id))
    end
    return out
end

"""
Interpret one `in_network` cell.

CSV.jl types the column from its contents, so the same generated file can
arrive as Int, Float64, Bool or String depending on what else is in it, and a
round-trip through another tool can quote it. An unrecognised value used to
reach `Int(...)` and die with a bare `MethodError` naming neither the file nor
the column; a corrupt artifact should say which artifact and which field.

A MISSING flag counts as in-network: the column is an optimisation that lets a
consumer skip rows, and a row present in the file is a cofactor either way.
"""
function _in_network_flag(value, path::String)::Bool
    ismissing(value) && return true
    value isa Bool && return value
    value isa Real && return value != 0
    if value isa AbstractString
        text = strip(value)
        isempty(text) && return true
        parsed = tryparse(Float64, text)
        parsed === nothing || return parsed != 0
        lowered = lowercase(text)
        lowered in ("true", "yes") && return true
        lowered in ("false", "no") && return false
    end
    throw(ArgumentError(
        "$(path): could not read `in_network` value $(repr(value)); " *
        "expected 0/1, true/false, or empty"))
end