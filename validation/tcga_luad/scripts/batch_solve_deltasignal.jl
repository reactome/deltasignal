#!/usr/bin/env julia

# Batch DeltaSignal steady-state solves for one parsed network.
#
# It loads one parsed network JSON once, then solves many observation CSVs
# from a manifest.

import Pkg

default_deltasignal_dir = abspath(joinpath(@__DIR__, "..", "..", ".."))
deltasignal_dir = get(ENV, "DELTASIGNAL_DIR", default_deltasignal_dir)
Pkg.activate(deltasignal_dir)

using CSV
using DataFrames
using Dates
using JSON3

push!(LOAD_PATH, joinpath(deltasignal_dir, "src"))
using DeltaSignal

const DS_CONFIG_DEFAULTS = Dict(
    "DS_AND_MODE" => "hill_log",
    "DS_ASSEMBLY_LIMITING" => "1",
    "DS_DAMPING" => "0.0",
    "DS_DEPLETION_H_MAX" => "10.0",
    "DS_DROP_PASSTHROUGH" => "0",
    "DS_GENE_STIDS_FILE" => "",
    "DS_HILL_LOG_ZMAX" => "10.0",
    "DS_HILL_SAT_EPS" => "0.001",
    "DS_HILL_SAT_H_MAX" => "10.0",
    "DS_INHIBITION_MODE" => "divide",
    "DS_INHIBITOR_BETA" => "1.0",
    "DS_INHIBITOR_EPS" => "0.001",
    "DS_INHIBITOR_FLOOR" => "0.0",
    "DS_INHIBITOR_FLOOR_SCOPE" => "loops",
    "DS_INHIBITOR_K" => "0.1",
    "DS_LOOP_DEPTH" => "3",
    "DS_OR_MODE" => "mean",
    "DS_SCC_BREAK_CATALYST" => "0",
    "DS_SCC_DAMPING" => "0.5",
    "DS_SCC_NEG_FRAC" => "0.5",
    "DS_SCC_NEG_ITERS" => "1",
    "DS_SCC_NEG_MODE" => "converge",
    "DS_SCC_SOLVE" => "1",
)


function effective_solver_config()
    return Dict(
        key => Dict(
            "value" => get(ENV, key, default),
            "source" => haskey(ENV, key) ? "environment" : "code_default",
        )
        for (key, default) in DS_CONFIG_DEFAULTS
    )
end


function git_commit(repo_dir::String)
    try
        return readchomp(`git -C $repo_dir rev-parse HEAD`)
    catch
        return nothing
    end
end


function parse_cli_args(args)
    parsed = Dict{String, String}()
    i = 1
    while i <= length(args)
        key = args[i]
        if !startswith(key, "--")
            error("Unexpected positional argument: $key")
        end
        if i == length(args)
            error("Missing value for argument: $key")
        end
        parsed[key[3:end]] = args[i + 1]
        i += 2
    end

    required = ["network", "manifest", "summary"]
    missing = [key for key in required if !haskey(parsed, key)]
    if !isempty(missing)
        error("Missing required arguments: $(join(missing, ", "))")
    end
    return parsed
end


function optional_string(value)
    if value === nothing || value === missing
        return nothing
    end
    text = String(value)
    if text in ["", "None", "nan", "NaN", "<NA>"]
        return nothing
    end
    return text
end


function load_network(network_path::String)::ReactionNetwork
    network_json = JSON3.read(read(network_path, String))

    nodes = Dict{String, NetworkNode}()
    for (uuid, node_data) in network_json["nodes"]
        nodes[String(uuid)] = NetworkNode(
            String(node_data["uuid"]),
            optional_string(node_data["reactome_id"]),
            String(node_data["entity_type"]),
            optional_string(node_data["original_set_id"]),
            String(node_data["display_name"]),
            Float64(node_data["baseline"]),
        )
    end

    edges = LogicNetworkEdge[]
    for edge_data in network_json["edges"]
        edge_type = haskey(edge_data, "edge_type") && edge_data["edge_type"] !== nothing ?
            String(edge_data["edge_type"]) : ""
        push!(edges, LogicNetworkEdge(
            String(edge_data["parent_uuid"]),
            String(edge_data["child_uuid"]),
            Bool(edge_data["is_and"]),
            Bool(edge_data["is_positive"]),
            Float64(edge_data["stoichiometry"]),
            edge_type,
        ))
    end

    set_mappings = Dict{String, SetExpansionMapping}()
    if haskey(network_json, "set_mappings")
        for (set_id, mapping_data) in network_json["set_mappings"]
            set_mappings[String(set_id)] = SetExpansionMapping(
                String(mapping_data["original_set_id"]),
                String(mapping_data["original_name"]),
                [String(member) for member in mapping_data["expanded_members"]],
                nothing,
            )
        end
    end

    return ReactionNetwork(nodes, edges, set_mappings)
end


function load_observations(path::String)::Dict{String, Tuple{Float64, Float64}}
    observations_df = CSV.read(path, DataFrame)
    observations = Dict{String, Tuple{Float64, Float64}}()
    for row in eachrow(observations_df)
        uuid = hasproperty(row, :node_uuid) ? String(row.node_uuid) : String(row.node)
        activity = hasproperty(row, :activity) ? Float64(row.activity) : Float64(row.value)
        confidence = Float64(row.confidence)
        observations[uuid] = (activity, confidence)
    end
    return observations
end


function write_result(
    output_path::String,
    result::SolverResult,
    observations::Dict{String, Tuple{Float64, Float64}},
    params::SteadyStateParams,
    solver_config::Dict,
    deltasignal_commit,
)
    activities_ui = Dict(uuid => activity * 100.0 for (uuid, activity) in result.node_activities)
    output_dict = Dict(
        "node_activities" => activities_ui,
        "influence_scores" => Dict{String, Float64}(),
        "solver_info" => Dict(
            "converged" => result.converged,
            "iterations" => result.iterations,
            "final_residual" => result.final_residual,
            "solve_time" => result.solve_time,
            "diagnostics" => result.diagnostics,
        ),
        "observations" => Dict(
            uuid => Dict("activity" => act, "confidence" => conf)
            for (uuid, (act, conf)) in observations
        ),
        "parameters" => Dict(
            "mu" => params.mu,
            "gamma" => params.gamma,
            "max_iters" => params.max_iters,
            "tolerance" => params.tolerance,
            "aggregation" => "not_computed",
            "batch_solver" => true,
            "influence_scores_computed" => false,
        ),
        "provenance" => Dict(
            "deltasignal_commit" => deltasignal_commit,
            "effective_solver_config" => solver_config,
        ),
    )
    open(output_path, "w") do handle
        JSON3.pretty(handle, output_dict)
    end
end


function write_sample_log(log_path::String, sample_id::String, result::SolverResult, elapsed::Float64)
    open(log_path, "w") do handle
        println(handle, "Batch DeltaSignal solve")
        println(handle, "sample_id: $sample_id")
        println(handle, "converged: $(result.converged)")
        println(handle, "iterations: $(result.iterations)")
        println(handle, "final_residual: $(result.final_residual)")
        println(handle, "reported_solve_time: $(result.solve_time)")
        println(handle, "batch_wall_seconds: $elapsed")
        println(handle, "timestamp: $(Dates.now())")
    end
end


function main()
    args = parse_cli_args(ARGS)
    network_path = args["network"]
    manifest_path = args["manifest"]
    summary_path = args["summary"]
    mu = parse(Float64, get(args, "mu", "1.0"))
    gamma = parse(Float64, get(args, "gamma", "0.1"))
    max_iters = parse(Int, get(args, "max-iters", "2000"))
    tolerance = parse(Float64, get(args, "tolerance", "1e-6"))
    params = SteadyStateParams(mu, gamma, max_iters, tolerance, "penalty")
    solver_config = effective_solver_config()
    deltasignal_commit = git_commit(deltasignal_dir)

    batch_start = time()
    network = load_network(network_path)
    manifest = CSV.read(manifest_path, DataFrame)
    required_cols = Set(["sample_id", "observations_csv", "result_json", "solve_log"])
    missing_cols = setdiff(required_cols, Set(String.(names(manifest))))
    if !isempty(missing_cols)
        error("Manifest missing required columns: $(join(collect(missing_cols), ", "))")
    end

    rows = Vector{Dict{String, Any}}()
    for row in eachrow(manifest)
        sample_id = String(row.sample_id)
        observations_csv = String(row.observations_csv)
        result_json = String(row.result_json)
        solve_log = String(row.solve_log)
        sample_start = time()
        status = "ok"
        error_message = ""
        observations_count = 0
        converged = false
        iterations = 0
        final_residual = NaN
        solve_time = NaN

        try
            observations = load_observations(observations_csv)
            observations_count = length(observations)
            result = solve_steady_state(network, observations, params)
            write_result(
                result_json,
                result,
                observations,
                params,
                solver_config,
                deltasignal_commit,
            )
            elapsed = time() - sample_start
            write_sample_log(solve_log, sample_id, result, elapsed)
            converged = result.converged
            iterations = result.iterations
            final_residual = result.final_residual
            solve_time = result.solve_time
        catch err
            status = "error"
            error_message = sprint(showerror, err, catch_backtrace())
            open(solve_log, "w") do handle
                println(handle, error_message)
            end
        end

        push!(rows, Dict(
            "sample_id" => sample_id,
            "status" => status,
            "error" => error_message,
            "observations_count" => observations_count,
            "converged" => converged,
            "iterations" => iterations,
            "final_residual" => final_residual,
            "solve_time" => solve_time,
            "wall_seconds" => time() - sample_start,
            "result_json" => result_json,
            "observations_csv" => observations_csv,
            "solve_log" => solve_log,
        ))
    end

    summary = Dict(
        "network" => network_path,
        "manifest" => manifest_path,
        "sample_count" => nrow(manifest),
        "ok_count" => count(row -> row["status"] == "ok", rows),
        "error_count" => count(row -> row["status"] == "error", rows),
        "converged_count" => count(row -> row["converged"] == true, rows),
        "max_iters" => max_iters,
        "tolerance" => tolerance,
        "mu" => mu,
        "gamma" => gamma,
        "deltasignal_commit" => deltasignal_commit,
        "effective_solver_config" => solver_config,
        "batch_wall_seconds" => time() - batch_start,
        "results" => rows,
    )
    open(summary_path, "w") do handle
        JSON3.pretty(handle, summary)
    end
end


main()
