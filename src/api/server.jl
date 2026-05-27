# DeltaSignal API Server
using DeltaSignal
using HTTP
using JSON3
using ArgParse
using CSV
using DataFrames

# Parse command line arguments
function parse_commandline()
    s = ArgParseSettings()
    
    @add_arg_table s begin
        "--port", "-p"
            help = "Port to run server on"
            arg_type = Int
            default = 8080
        "--host"
            help = "Host to bind to"
            arg_type = String
            default = "127.0.0.1"
        "--test"
            help = "Test mode - exit after startup check"
            action = :store_true
    end
    
    return parse_args(s)
end

# CORS middleware
function cors_middleware(handler)
    return function(req)
        # Add CORS headers
        headers = [
            "Access-Control-Allow-Origin" => "*",
            "Access-Control-Allow-Methods" => "GET, POST, PUT, DELETE, OPTIONS",
            "Access-Control-Allow-Headers" => "Content-Type, Authorization",
            "Access-Control-Max-Age" => "86400"
        ]
        
        if req.method == "OPTIONS"
            return HTTP.Response(200, headers)
        end
        
        response = handler(req)
        
        # Add CORS headers to response
        for (key, value) in headers
            HTTP.setheader(response, key => value)
        end
        
        return response
    end
end

# Request logging middleware
function logging_middleware(handler)
    return function(req)
        start_time = time()
        println("$(time()) $(req.method) $(req.target)")
        
        response = handler(req)
        
        elapsed = round((time() - start_time) * 1000, digits=2)
        println("$(time()) $(response.status) $(elapsed)ms")
        
        return response
    end
end

# Health check endpoint
function health_handler(req)
    health_info = Dict(
        "status" => "ok",
        "timestamp" => string(time()),
        "version" => "0.1.0",
        "service" => "deltasignal-api",
        "julia_version" => string(VERSION)
    )
    
    return HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(health_info))
end

# Parse network endpoint
# --- Error handling helpers ---
#
# Raw Julia exception strings are never returned to clients (they leak
# internals and read poorly to users). Map known exception types to a
# user-safe (status, message) pair; log the original server-side.

const JSON_HEADERS = ["Content-Type" => "application/json"]

function user_facing_error(e::Exception)::Tuple{Int, String}
    if isa(e, ArgumentError)
        msg = e.msg
        if occursin("not found", msg) && occursin("column", msg)
            return 400, "Uploaded file has an unexpected column schema. Verify it matches the expected TSV format."
        elseif occursin("Unknown pathway id", msg) ||
               occursin("Invalid pathway id", msg) ||
               occursin("missing required files", msg)
            return 400, msg
        elseif occursin("not found", msg) || occursin("does not exist", msg)
            return 400, "Required input not found."
        end
        return 400, "Invalid input."
    elseif isa(e, SystemError)
        return 400, "Could not read input file."
    elseif isa(e, KeyError) || isa(e, BoundsError)
        return 400, "Required field missing from input."
    end
    return 500, "Internal server error."
end

function error_response(e::Exception; context::String="")
    @error "API request failed" exception=(e, catch_backtrace()) context=context
    status, msg = user_facing_error(e)
    body = Dict("status" => "error", "message" => msg)
    return HTTP.Response(status, JSON_HEADERS, JSON3.write(body))
end

const PARSE_FIELDS = ("logic_network", "uuid_mapping", "set_mappings", "observations")

# Reconstruct a ReactionNetwork from the JSON shape parse_handler returns,
# so /api/solve can use a network the client already parsed (no need to
# re-upload TSVs on every solve).
function reaction_network_from_json(data)::DeltaSignal.ReactionNetwork
    nodes_dict = Dict{String, DeltaSignal.NetworkNode}()
    for n in data.nodes
        reactome_id_raw = get(n, :reactome_id, nothing)
        set_id_raw = get(n, :set_id, nothing)
        nodes_dict[String(n.uuid)] = DeltaSignal.NetworkNode(
            String(n.uuid),
            isnothing(reactome_id_raw) ? nothing : String(reactome_id_raw),
            String(n.entity_type),
            isnothing(set_id_raw) ? nothing : String(set_id_raw),
            String(n.name),
            Float64(n.baseline),
        )
    end

    edges = DeltaSignal.LogicNetworkEdge[]
    for e in data.edges
        push!(edges, DeltaSignal.LogicNetworkEdge(
            String(e.parent_uuid),
            String(e.child_uuid),
            Bool(e.is_and),
            Bool(e.is_positive),
            Float64(e.stoichiometry),
        ))
    end

    set_mappings = Dict{String, DeltaSignal.SetExpansionMapping}()
    pathways_raw = get(data, :pathways, nothing)
    if pathways_raw !== nothing
        for p in pathways_raw
            members = String[String(m) for m in p.members]
            set_mappings[String(p.id)] = DeltaSignal.SetExpansionMapping(
                String(p.id),
                String(p.name),
                members,
                nothing,
            )
        end
    end

    return DeltaSignal.ReactionNetwork(nodes_dict, edges, set_mappings)
end

function extract_uploaded_files(req)
    multiparts = try
        HTTP.parse_multipart_form(req)
    catch
        nothing
    end
    multiparts === nothing && return nothing

    tmp_dir = mktempdir()
    paths = Dict{String, String}()
    for part in multiparts
        part.name in PARSE_FIELDS || continue
        target = joinpath(tmp_dir, "$(part.name)_$(something(part.filename, "uploaded"))")
        open(target, "w") do io
            write(io, read(part.data))
        end
        paths[part.name] = target
    end

    if !haskey(paths, "logic_network") || !haskey(paths, "uuid_mapping")
        rm(tmp_dir, recursive=true, force=true)
        return nothing
    end

    return (paths=paths, tmp_dir=tmp_dir)
end

# --- Pathway catalog ---
#
# Pre-generated Reactome pathway networks from the logic-network-generator
# project. Volume-mounted into the container at CATALOG_DIR by
# docker-compose.dev.yml. Each subdirectory is one pathway and contains
# logic_network.csv + stid_to_uuid_mapping.csv (generator format).

const CATALOG_DIR = "/app/pathway_catalog"
const SAMPLE_DIR = "examples"

# --- Reactome ContentService enrichment ---
#
# Generator-format networks have no display names — nodes are addressed by
# Reactome stable id only. We fetch (displayName, className) from the
# ContentService and patch the parsed nodes before sending to the client.
#
# Cache lives for the container's lifetime; failures degrade gracefully
# (the API returns whatever it could parse, with stable_ids as names).

const REACTOME_BATCH_URL = "https://reactome.org/ContentService/data/query/ids"
const REACTOME_TIMEOUT_SECONDS = 15
# ContentService silently truncates the response at 20 entries per batch
# regardless of how many IDs you send. Empirically verified May 2026.
const REACTOME_BATCH_SIZE = 20

const REACTOME_CACHE = Dict{String, NamedTuple{(:name, :entity_type), Tuple{String, String}}}()

function fetch_reactome_batch!(ids::Vector{String})
    isempty(ids) && return
    try
        response = HTTP.post(
            REACTOME_BATCH_URL,
            ["Content-Type" => "text/plain"],
            join(ids, ",");
            readtimeout = REACTOME_TIMEOUT_SECONDS,
            retry = false,
        )
        response.status == 200 || return
        for entry in JSON3.read(String(response.body))
            stid = String(entry.stId)
            name = haskey(entry, :displayName) ? String(entry.displayName) : stid
            etype = haskey(entry, :className) ? lowercase(String(entry.className)) : "unknown"
            REACTOME_CACHE[stid] = (name=name, entity_type=etype)
        end
    catch e
        @warn "Reactome ContentService batch lookup failed; nodes will keep stable_id placeholder names" exception=e batch_size=length(ids)
    end
end

# Enrich nodes whose display_name still equals their reactome_id (the
# placeholder we set in parse_uuid_mapping_generator). Nodes that already
# have a real name (e.g. from the bundled sample TSVs) are left untouched.
function enrich_with_reactome_names!(nodes::Dict{String, DeltaSignal.NetworkNode})
    needed = String[]
    for node in values(nodes)
        node.reactome_id === nothing && continue
        node.display_name == node.reactome_id || continue
        haskey(REACTOME_CACHE, node.reactome_id) && continue
        push!(needed, node.reactome_id)
    end

    if !isempty(needed)
        unique!(needed)
        for chunk_start in 1:REACTOME_BATCH_SIZE:length(needed)
            chunk_end = min(chunk_start + REACTOME_BATCH_SIZE - 1, length(needed))
            fetch_reactome_batch!(needed[chunk_start:chunk_end])
        end
    end

    for (uuid, node) in nodes
        node.reactome_id === nothing && continue
        node.display_name == node.reactome_id || continue
        cached = get(REACTOME_CACHE, node.reactome_id, nothing)
        cached === nothing && continue
        nodes[uuid] = DeltaSignal.NetworkNode(
            node.uuid,
            node.reactome_id,
            cached.entity_type,
            node.original_set_id,
            cached.name,
            node.baseline,
        )
    end
end

function sample_network_paths()
    return (
        joinpath(SAMPLE_DIR, "sample_logic_network.tsv"),
        joinpath(SAMPLE_DIR, "sample_uuid_mapping.tsv"),
        joinpath(SAMPLE_DIR, "sample_set_mappings.tsv"),
    )
end

function catalog_network_paths(pathway_id::String)
    # Guard against path traversal — pathway_id must be a single component
    if occursin("/", pathway_id) || pathway_id == "." || pathway_id == ".."
        throw(ArgumentError("Invalid pathway id."))
    end
    dir = joinpath(CATALOG_DIR, pathway_id)
    if !isdir(dir)
        throw(ArgumentError("Unknown pathway id: $pathway_id"))
    end
    logic = joinpath(dir, "logic_network.csv")
    uuid = joinpath(dir, "stid_to_uuid_mapping.csv")
    if !isfile(logic) || !isfile(uuid)
        throw(ArgumentError("Pathway $pathway_id is missing required files."))
    end
    return (logic, uuid, nothing)  # generator has no set_mappings
end

# Build a deduplicated, pretty-named catalog from CATALOG_DIR contents.
# Directory names look like "Cell_Cycle_Checkpoints_69620" or
# "Cell_Cycle_Checkpoints_R-HSA-69620". Group by trailing numeric id;
# prefer the R-HSA variant when both exist.
function pathway_catalog_entries()
    isdir(CATALOG_DIR) || return Vector{Dict{String, String}}()

    groups = Dict{String, Vector{NamedTuple{(:id, :stable_id, :name, :is_rhsa)}}}()
    for entry in readdir(CATALOG_DIR; sort=true)
        isdir(joinpath(CATALOG_DIR, entry)) || continue

        m_rhsa = match(r"^(.*)_R-HSA-(\d+)$", entry)
        if m_rhsa !== nothing
            pretty = replace(m_rhsa.captures[1], "_" => " ")
            numeric = m_rhsa.captures[2]
            stable = "R-HSA-" * numeric
            push!(get!(groups, numeric, NamedTuple{(:id, :stable_id, :name, :is_rhsa)}[]),
                  (id=entry, stable_id=stable, name=pretty, is_rhsa=true))
            continue
        end

        m_plain = match(r"^(.*)_(\d+)$", entry)
        if m_plain !== nothing
            pretty = replace(m_plain.captures[1], "_" => " ")
            numeric = m_plain.captures[2]
            push!(get!(groups, numeric, NamedTuple{(:id, :stable_id, :name, :is_rhsa)}[]),
                  (id=entry, stable_id=numeric, name=pretty, is_rhsa=false))
        end
    end

    entries = Dict{String, String}[]
    for numeric in sort(collect(keys(groups)))
        candidates = groups[numeric]
        rhsa_idx = findfirst(c -> c.is_rhsa, candidates)
        chosen = rhsa_idx === nothing ? candidates[1] : candidates[rhsa_idx]
        push!(entries, Dict(
            "id" => chosen.id,
            "stable_id" => chosen.stable_id,
            "name" => chosen.name,
        ))
    end

    sort!(entries, by = e -> e["name"])
    return entries
end

function pathways_handler(req)
    try
        return HTTP.Response(200, JSON_HEADERS, JSON3.write(pathway_catalog_entries()))
    catch e
        return error_response(e; context="pathways_handler")
    end
end

function parse_handler(req)
    upload = extract_uploaded_files(req)
    try
        local logic_network_path, uuid_mapping_path, set_mapping_path

        if upload !== nothing
            logic_network_path = upload.paths["logic_network"]
            uuid_mapping_path = upload.paths["uuid_mapping"]
            set_mapping_path = get(upload.paths, "set_mappings", nothing)
        else
            # Check for JSON body with a catalog pathway_id
            body = String(req.body)
            content_type = HTTP.header(req, "Content-Type", "")
            pathway_id = nothing
            if !isempty(body) && occursin("application/json", lowercase(content_type))
                request_data = JSON3.read(body)
                pid_raw = get(request_data, :pathway_id, nothing)
                if pid_raw !== nothing
                    pathway_id = String(pid_raw)
                end
            end

            if pathway_id !== nothing
                logic_network_path, uuid_mapping_path, set_mapping_path =
                    catalog_network_paths(pathway_id)
            else
                # Final fallback: bundled sample (back-compat)
                logic_network_path, uuid_mapping_path, set_mapping_path = sample_network_paths()
            end
        end

        # Parse the network using DeltaSignal functions
        network = DeltaSignal.parse_complete_network(
            logic_network_path,
            uuid_mapping_path,
            set_mapping_path
        )

        # Replace placeholder names with real Reactome display names + types.
        # No-op if the network already has them or if ContentService is down.
        enrich_with_reactome_names!(network.nodes)

        # Convert to API response format
        nodes_array = []
        for (uuid, node) in network.nodes
            push!(nodes_array, Dict(
                "uuid" => node.uuid,
                "name" => node.display_name,
                "reactome_id" => node.reactome_id,
                "entity_type" => node.entity_type,
                "baseline" => node.baseline,
                "set_id" => node.original_set_id
            ))
        end
        
        edges_array = []
        for edge in network.edges
            push!(edges_array, Dict(
                "parent_uuid" => edge.parent_uuid,
                "child_uuid" => edge.child_uuid,
                "is_and" => edge.is_and,
                "is_positive" => edge.is_positive,
                "stoichiometry" => edge.stoichiometry
            ))
        end
        
        pathways_array = []
        for (set_id, mapping) in network.set_mappings
            push!(pathways_array, Dict(
                "id" => mapping.original_set_id,
                "name" => mapping.original_name,
                "members" => mapping.expanded_members
            ))
        end
        
        result = Dict(
            "status" => "success",
            "message" => "Network parsed successfully",
            "nodes" => nodes_array,
            "edges" => edges_array,
            "pathways" => pathways_array
        )
        
        return HTTP.Response(200, JSON_HEADERS, JSON3.write(result))

    catch e
        return error_response(e; context="parse_handler")
    finally
        if upload !== nothing
            rm(upload.tmp_dir, recursive=true, force=true)
        end
    end
end

# Solve steady-state endpoint
function solve_handler(req)
    try
        # Parse request body
        body = String(req.body)
        if isempty(body)
            return HTTP.Response(400, JSON_HEADERS,
                JSON3.write(Dict("status" => "error", "message" => "Request body is required.")))
        end
        
        # Parse the JSON request
        request_data = JSON3.read(body)

        # If the client sent the previously-parsed network in the request,
        # use it. Otherwise fall back to the bundled sample so older clients
        # / smoke tests still work.
        network_raw = get(request_data, :network, nothing)
        if network_raw !== nothing
            network = reaction_network_from_json(network_raw)
        else
            examples_dir = "examples"
            network = DeltaSignal.parse_complete_network(
                joinpath(examples_dir, "sample_logic_network.tsv"),
                joinpath(examples_dir, "sample_uuid_mapping.tsv"),
                joinpath(examples_dir, "sample_set_mappings.tsv"),
            )
        end
        
        # Create observations from request perturbations
        observations = Dict{String, Tuple{Float64, Float64}}()
        
        if haskey(request_data, :observations) && request_data.observations !== nothing
            for (node_uuid, obs_data) in request_data.observations
                # obs_data contains [activity, confidence] 
                activity = Float64(obs_data[1])
                confidence = Float64(obs_data[2])
                observations[String(node_uuid)] = (activity, confidence)
            end
        end
        
        println("Solving with ", length(observations), " observations")
        
        # Solve the steady state
        solver_result = DeltaSignal.solve_steady_state(network, observations)
        
        # Keep node activities in internal 0-1 scale for frontend processing
        node_activities_display = Dict{String, Float64}()
        for (uuid, activity) in solver_result.node_activities
            node_activities_display[uuid] = activity  # Keep in 0-1 scale
        end
        
        # Create influence scores (simplified - could be enhanced)
        reactions = DeltaSignal.convert_to_reaction_network(network)
        influence_scores = DeltaSignal.compute_influence_scores(solver_result, reactions)
        
        result = Dict(
            "status" => "success",
            "message" => "Steady-state solved successfully",
            "node_activities" => node_activities_display,
            "influence_scores" => influence_scores,
            "converged" => solver_result.converged,
            "iterations" => solver_result.iterations,
            "solve_time" => solver_result.solve_time
        )
        
        return HTTP.Response(200, JSON_HEADERS, JSON3.write(result))

    catch e
        return error_response(e; context="solve_handler")
    end
end

# Default 404 handler
function not_found_handler(req)
    body = Dict("status" => "error", "message" => "Endpoint not found: $(req.target)")
    return HTTP.Response(404, JSON_HEADERS, JSON3.write(body))
end

# Main router
function create_router()
    router = HTTP.Router()
    
    # API endpoints
    HTTP.register!(router, "GET", "/api/health", health_handler)
    HTTP.register!(router, "GET", "/api/pathways", pathways_handler)
    HTTP.register!(router, "POST", "/api/parse", parse_handler)
    HTTP.register!(router, "POST", "/api/solve", solve_handler)
    
    # Catch-all for unknown routes
    HTTP.register!(router, "*", "*", not_found_handler)
    
    return router
end

# Start server
function start_server(host="127.0.0.1", port=8080; test_mode=false)
    println("🧬 Starting DeltaSignal API server...")
    println("Julia version: $(VERSION)")
    println("DeltaSignal module loaded: $(isdefined(Main, :DeltaSignal))")
    
    if test_mode
        println("Test mode - server would start on http://$host:$port")
        println("Available endpoints:")
        println("  GET  /api/health")
        println("  POST /api/parse")
        println("  POST /api/solve")
        return
    end
    
    try
        router = create_router()
        
        # Apply middleware
        handler = cors_middleware(logging_middleware(router))
        
        println("Server starting on http://$host:$port")
        println("Available endpoints:")
        println("  GET  /api/health")
        println("  POST /api/parse")
        println("  POST /api/solve")
        println("Press Ctrl+C to stop")
        
        HTTP.serve(handler, host, port; verbose=false)
        
    catch e
        if isa(e, InterruptException)
            println("\n🛑 Server stopped")
        else
            println("❌ Server error: $e")
            rethrow(e)
        end
    end
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    args = parse_commandline()
    start_server(args["host"], args["port"]; test_mode=args["test"])
end