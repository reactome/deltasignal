#!/usr/bin/env julia

# DeltaSignal Command Line Interface

# Activate the project environment
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using ArgParse
using JSON3
using CSV
using DataFrames

# Add the src directory to the load path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using DeltaSignal

function parse_commandline()
    # Check if we have at least one argument
    if length(ARGS) == 0
        println("Error: No command provided")
        println("Available commands: parse, solve, rollout, train, validate, export, server")
        println("\nUse 'deltasignal <command> --help' for more information")
        exit(1)
    end

    # Get the command from first argument
    command = ARGS[1]

    # Create command-specific argument tables
    if command == "parse"
        return parse_parse_args()
    elseif command == "solve"
        return parse_solve_args()
    elseif command == "rollout"
        return parse_rollout_args()
    elseif command == "train"
        return parse_train_args()
    elseif command == "validate"
        return parse_validate_args()
    elseif command == "export"
        return parse_export_args()
    elseif command == "server"
        return parse_server_args()
    elseif command == "--version" || command == "-v"
        println("DeltaSignal version 0.1.0")
        exit(0)
    elseif command == "--help" || command == "-h"
        println("DeltaSignal - A Pathway Perturbation & Dynamics Engine")
        println("\nAvailable commands:")
        println("  parse      - Parse logic network from TSV files")
        println("  solve      - Solve steady-state network")
        println("  rollout    - Run time-dynamic rollout")
        println("  train      - Train parameters from data")
        println("  validate   - Validate model against test data")
        println("  export     - Export results in various formats")
        println("  server     - Start web API server")
        println("\nUse 'deltasignal <command> --help' for more information on a specific command")
        exit(0)
    else
        println("Unknown command: $command")
        println("Available commands: parse, solve, rollout, train, validate, export, server")
        println("\nUse 'deltasignal --help' for more information")
        exit(1)
    end
end

function parse_parse_args()
    s = ArgParseSettings(autofix_names=true)
    @add_arg_table! s begin
        "--logic", "-l"
            help = "Logic network TSV file"
            arg_type = String
            required = true
        "--uuid-map", "-u"
            help = "UUID to Reactome ID mapping TSV file"
            arg_type = String
            required = true
        "--set-map", "-s"
            help = "Set expansion mapping TSV file (optional)"
            arg_type = String
        "--output", "-o"
            help = "Output JSON file"
            arg_type = String
            required = true
        "--validate"
            help = "Validate the parsed network"
            action = :store_true
    end
    return merge(Dict(:command => "parse"), parse_args(ARGS[2:end], s, as_symbols=true))
end

function parse_solve_args()
    s = ArgParseSettings(autofix_names=true)
    @add_arg_table! s begin
        "--network", "-n"
            help = "Network JSON file"
            arg_type = String
            required = true
        "--observations"
            help = "Observations CSV file"
            arg_type = String
            required = true
        "--mode", "-m"
            help = "Solver mode"
            arg_type = String
            default = "steady-state"
            range_tester = x -> x in ["steady-state", "ss"]
        "--output", "-o"
            help = "Output results JSON file"
            arg_type = String
            required = true
        "--output-pathway"
            help = "Output aggregated pathway results JSON file"
            arg_type = String
        "--mu"
            help = "Model consistency penalty weight"
            arg_type = Float64
            default = 1.0
        "--gamma"
            help = "Baseline prior weight"
            arg_type = Float64
            default = 0.1
        "--max-iters"
            help = "Maximum optimization iterations"
            arg_type = Int
            default = 500
        "--tolerance"
            help = "Convergence tolerance"
            arg_type = Float64
            default = 1e-6
        "--aggregation"
            help = "Pathway aggregation method"
            arg_type = String
            default = "stoichiometry_weighted"
            range_tester = x -> x in ["mean", "max", "min", "stoichiometry_weighted", "confidence_weighted", "geometric_mean"]
    end
    return merge(Dict(:command => "solve"), parse_args(ARGS[2:end], s, as_symbols=true))
end

function parse_rollout_args()
    s = ArgParseSettings(autofix_names=true)
    @add_arg_table! s begin
        "--network", "-n"
            help = "Network JSON file"
            arg_type = String
            required = true
        "--initial", "-i"
            help = "Initial conditions CSV file"
            arg_type = String
            required = true
        "--steps", "-T"
            help = "Number of time steps"
            arg_type = Int
            default = 100
        "--output", "-o"
            help = "Output trajectory JSON file"
            arg_type = String
            required = true
        "--damping"
            help = "Global damping factor for stability"
            arg_type = Float64
            default = 0.1
    end
    return merge(Dict(:command => "rollout"), parse_args(ARGS[2:end], s, as_symbols=true))
end

function parse_train_args()
    s = ArgParseSettings(autofix_names=true)
    @add_arg_table! s begin
        "--data-dir", "-d"
            help = "Training data directory"
            arg_type = String
            required = true
        "--network", "-n"
            help = "Network JSON file"
            arg_type = String
            required = true
        "--output", "-o"
            help = "Output trained parameters JSON file"
            arg_type = String
            required = true
        "--epochs"
            help = "Number of training epochs"
            arg_type = Int
            default = 100
        "--learning-rate"
            help = "Learning rate"
            arg_type = Float64
            default = 0.01
        "--regularization"
            help = "L2 regularization weight"
            arg_type = Float64
            default = 0.01
    end
    return merge(Dict(:command => "train"), parse_args(ARGS[2:end], s, as_symbols=true))
end

function parse_validate_args()
    s = ArgParseSettings(autofix_names=true)
    @add_arg_table! s begin
        "--network", "-n"
            help = "Network JSON file"
            arg_type = String
            required = true
        "--test-data", "-t"
            help = "Test data directory"
            arg_type = String
            required = true
        "--params", "-p"
            help = "Parameters JSON file"
            arg_type = String
        "--output", "-o"
            help = "Validation report JSON file"
            arg_type = String
            required = true
    end
    return merge(Dict(:command => "validate"), parse_args(ARGS[2:end], s, as_symbols=true))
end

function parse_export_args()
    s = ArgParseSettings(autofix_names=true)
    @add_arg_table! s begin
        "--results", "-r"
            help = "Results JSON file"
            arg_type = String
            required = true
        "--format", "-f"
            help = "Export format"
            arg_type = String
            default = "pathway-browser"
            range_tester = x -> x in ["pathway-browser", "cytoscape", "csv"]
        "--output", "-o"
            help = "Output file"
            arg_type = String
            required = true
    end
    return merge(Dict(:command => "export"), parse_args(ARGS[2:end], s, as_symbols=true))
end

function parse_server_args()
    s = ArgParseSettings(autofix_names=true)
    @add_arg_table! s begin
        "--port", "-p"
            help = "Server port"
            arg_type = Int
            default = 8000
        "--host"
            help = "Server host"
            arg_type = String
            default = "localhost"
    end
    return merge(Dict(:command => "server"), parse_args(ARGS[2:end], s, as_symbols=true))
end

function execute_parse_command(args)
    println("🔍 Parsing logic network...")
    
    try
        network = parse_complete_network(
            args[:logic],
            args[:uuid_map],
            args[:set_map]
        )
        
        if args[:validate]
            println("\n📊 Validating network mapping...")
            validation_report = validate_network_pathway_mapping(network)
            
            println("Network validation:")
            println("- Total nodes: $(validation_report["total_network_nodes"])")
            println("- Nodes with Reactome IDs: $(validation_report["nodes_with_reactome_ids"])")
            println("- Orphaned nodes: $(length(validation_report["nodes_without_reactome_ids"]))")
            println("- Set expansions: $(validation_report["set_expansion_stats"]["total_sets"])")
            
            if !isempty(validation_report["warnings"])
                println("\n⚠️  Warnings:")
                for warning in validation_report["warnings"]
                    println("  - $warning")
                end
            end
        end
        
        # Save network to JSON
        network_dict = Dict(
            "nodes" => Dict(uuid => Dict(
                "uuid" => node.uuid,
                "reactome_id" => node.reactome_id,
                "entity_type" => node.entity_type,
                "original_set_id" => node.original_set_id,
                "display_name" => node.display_name,
                "baseline" => node.baseline
            ) for (uuid, node) in network.nodes),
            "edges" => [Dict(
                "parent_uuid" => edge.parent_uuid,
                "child_uuid" => edge.child_uuid,
                "is_and" => edge.is_and,
                "is_positive" => edge.is_positive,
                "stoichiometry" => edge.stoichiometry,
                "edge_type" => edge.edge_type
            ) for edge in network.edges],
            "set_mappings" => Dict(set_id => Dict(
                "original_set_id" => mapping.original_set_id,
                "original_name" => mapping.original_name,
                "expanded_members" => mapping.expanded_members
            ) for (set_id, mapping) in network.set_mappings)
        )
        
        open(args[:output], "w") do f
            JSON3.pretty(f, network_dict)
        end
        
        println("✅ Network saved to: $(args[:output])")
        
    catch e
        println("❌ Error parsing network: $e")
        exit(1)
    end
end

function execute_solve_command(args)
    println("🧮 Solving steady-state network...")

    try
        # Load network JSON
        println("Loading network from: $(args[:network])")
        network_json = JSON3.read(read(args[:network], String))

        # Reconstruct ReactionNetwork from JSON
        nodes = Dict{String, NetworkNode}()
        for (uuid, node_data) in network_json["nodes"]
            nodes[String(uuid)] = NetworkNode(
                String(node_data["uuid"]),
                node_data["reactome_id"] !== nothing ? String(node_data["reactome_id"]) : nothing,
                String(node_data["entity_type"]),
                node_data["original_set_id"] !== nothing ? String(node_data["original_set_id"]) : nothing,
                String(node_data["display_name"]),
                Float64(node_data["baseline"])
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
                edge_type
            ))
        end

        set_mappings = Dict{String, SetExpansionMapping}()
        if haskey(network_json, "set_mappings")
            for (set_id, mapping_data) in network_json["set_mappings"]
                set_mappings[String(set_id)] = SetExpansionMapping(
                    String(mapping_data["original_set_id"]),
                    String(mapping_data["original_name"]),
                    [String(m) for m in mapping_data["expanded_members"]],
                    nothing
                )
            end
        end

        network = ReactionNetwork(nodes, edges, set_mappings)
        println("✓ Network loaded: $(length(nodes)) nodes, $(length(edges)) edges")

        # Load observations CSV
        println("Loading observations from: $(args[:observations])")
        observations_df = CSV.read(args[:observations], DataFrame)

        # Parse observations into Dict{String, Tuple{Float64, Float64}}
        observations = Dict{String, Tuple{Float64, Float64}}()
        for row in eachrow(observations_df)
            # Handle different possible column names
            uuid = if hasproperty(row, :node_uuid)
                String(row.node_uuid)
            else
                String(row.node)
            end

            # Get activity value (might be called "value" or "activity")
            # Values should already be in 0-100 scale (0=none, 1=1% of normal, 100=100% of normal)
            activity = if hasproperty(row, :activity)
                Float64(row.activity)
            else
                Float64(row.value)
            end

            confidence = Float64(row.confidence)
            observations[uuid] = (activity, confidence)
        end
        println("✓ Loaded $(length(observations)) observations")

        # Create solver parameters
        params = SteadyStateParams(
            args[:mu],
            args[:gamma],
            args[:max_iters],
            args[:tolerance],
            "penalty"  # Default to penalty method
        )

        # Solve steady-state
        println("\n🔬 Running steady-state solver...")
        println("   Solver parameters:")
        println("   - mu (model consistency): $(params.mu)")
        println("   - gamma (baseline prior): $(params.gamma)")
        println("   - max iterations: $(params.max_iters)")
        println("   - tolerance: $(params.tolerance)")

        result = solve_steady_state(network, observations, params)

        # Report results
        println("\n📊 Solver completed:")
        println("   - Converged: $(result.converged)")
        println("   - Iterations: $(result.iterations)")
        println("   - Final residual: $(round(result.final_residual, digits=6))")
        println("   - Solve time: $(round(result.solve_time, digits=3))s")

        # Compute influence scores
        reactions = convert_to_reaction_network(network)
        influence_scores = compute_influence_scores(result, reactions)

        # Convert activities back to 0-100 UI scale for output
        activities_ui = Dict(uuid => activity * 100.0 for (uuid, activity) in result.node_activities)

        # Prepare output
        output_dict = Dict(
            "node_activities" => activities_ui,
            "influence_scores" => influence_scores,
            "solver_info" => Dict(
                "converged" => result.converged,
                "iterations" => result.iterations,
                "final_residual" => result.final_residual,
                "solve_time" => result.solve_time,
                "diagnostics" => result.diagnostics
            ),
            "observations" => Dict(uuid => Dict("activity" => act, "confidence" => conf)
                                  for (uuid, (act, conf)) in observations),
            "parameters" => Dict(
                "mu" => params.mu,
                "gamma" => params.gamma,
                "max_iters" => params.max_iters,
                "tolerance" => params.tolerance,
                "aggregation" => args[:aggregation]
            )
        )

        # Save results
        open(args[:output], "w") do f
            JSON3.pretty(f, output_dict)
        end
        println("\n✅ Results saved to: $(args[:output])")

        # Show top influential nodes
        sorted_influence = sort(collect(influence_scores), by=x->x[2], rev=true)
        if !isempty(sorted_influence)
            println("\n🎯 Top 5 influential nodes:")
            for (i, (uuid, score)) in enumerate(sorted_influence[1:min(5, length(sorted_influence))])
                node_name = nodes[uuid].display_name
                println("   $i. $node_name ($uuid): $(round(score, digits=4))")
            end
        end

        # Generate pathway-level aggregation if requested
        if haskey(args, :output_pathway) && args[:output_pathway] !== nothing
            println("\n📦 Aggregating to pathway view...")
            # aggregate_to_pathway_view expects Vector{NetworkResult} + the
            # ReactionNetwork, with aggregation_method as a String keyword.
            # Build NetworkResults from the solved activities (internal 0-1),
            # carrying observation confidence and influence score where known.
            network_results = DeltaSignal.NetworkResult[
                DeltaSignal.NetworkResult(
                    uuid,
                    activity,
                    haskey(observations, uuid) ? observations[uuid][2] : 1.0,
                    get(influence_scores, uuid, 0.0),
                )
                for (uuid, activity) in result.node_activities
            ]
            pathway_results = aggregate_to_pathway_view(
                network_results,
                network;
                aggregation_method = String(args[:aggregation]),
            )

            open(args[:output_pathway], "w") do f
                JSON3.pretty(f, pathway_results)
            end
            println("✅ Pathway results saved to: $(args[:output_pathway])")
        end

    catch e
        println("❌ Error in solve: $e")
        println(stacktrace(catch_backtrace()))
        exit(1)
    end
end

function execute_rollout_command(args)
    println("⏱️  Running time-dynamic rollout...")
    
    try
        # TODO: Implement time-dynamic rollout
        println("⚠️  Time-dynamic rollout not yet implemented")
        
        println("✅ Trajectory would be saved to: $(args[:output])")
        
    catch e
        println("❌ Error in rollout: $e")
        exit(1)
    end
end

function execute_train_command(args)
    println("🎓 Training parameters...")
    
    try
        # TODO: Implement parameter learning
        println("⚠️  Parameter training not yet implemented")
        
        println("✅ Trained parameters would be saved to: $(args[:output])")
        
    catch e
        println("❌ Error in training: $e")
        exit(1)
    end
end

function execute_validate_command(args)
    println("🔬 Validating model...")
    
    try
        # TODO: Implement validation
        println("⚠️  Model validation not yet implemented")
        
        println("✅ Validation report would be saved to: $(args[:output])")
        
    catch e
        println("❌ Error in validation: $e")
        exit(1)
    end
end

function execute_export_command(args)
    println("📤 Exporting results...")

    try
        # Load results JSON
        println("Loading results from: $(args[:results])")
        results_json = JSON3.read(read(args[:results], String))

        format = args[:format]
        println("Export format: $format")

        if format == "pathway-browser"
            # Export for PathwayBrowser overlay
            export_pathway_browser_format(results_json, args[:output])

        elseif format == "cytoscape"
            # Export for Cytoscape.js visualization
            export_cytoscape_format(results_json, args[:output])

        elseif format == "csv"
            # Export as CSV table
            export_csv_format(results_json, args[:output])

        else
            println("❌ Unknown export format: $format")
            println("Available formats: pathway-browser, cytoscape, csv")
            exit(1)
        end

        println("✅ Exported to: $(args[:output])")

    catch e
        println("❌ Error in export: $e")
        println(stacktrace(catch_backtrace()))
        exit(1)
    end
end

function export_pathway_browser_format(results_json, output_path::String)
    # Create PathwayBrowser overlay format
    overlay = Dict(
        "type" => "activity_overlay",
        "metadata" => Dict(
            "solver" => results_json["solver_info"],
            "parameters" => results_json["parameters"]
        ),
        "entities" => []
    )

    # Add node activities
    for (uuid, activity) in results_json["node_activities"]
        influence = get(results_json["influence_scores"], uuid, 0.0)
        observed = haskey(results_json["observations"], uuid)

        entity = Dict(
            "id" => uuid,
            "activity" => round(activity, digits=2),
            "influence" => round(influence, digits=4),
            "observed" => observed
        )

        push!(overlay["entities"], entity)
    end

    # Write JSON
    open(output_path, "w") do f
        JSON3.pretty(f, overlay)
    end

    println("✓ PathwayBrowser overlay with $(length(overlay["entities"])) entities")
end

function export_cytoscape_format(results_json, output_path::String)
    # Create Cytoscape.js compatible JSON
    cytoscape_data = Dict(
        "elements" => Dict(
            "nodes" => [],
            "edges" => []
        ),
        "data" => Dict(
            "solver_info" => results_json["solver_info"],
            "parameters" => results_json["parameters"]
        )
    )

    # Add nodes with activities
    for (uuid, activity) in results_json["node_activities"]
        influence = get(results_json["influence_scores"], uuid, 0.0)
        observed = haskey(results_json["observations"], uuid)

        node = Dict(
            "data" => Dict(
                "id" => uuid,
                "activity" => round(activity, digits=2),
                "influence" => round(influence, digits=4),
                "observed" => observed,
                "label" => uuid
            ),
            "classes" => observed ? "observed" : "inferred"
        )

        push!(cytoscape_data["elements"]["nodes"], node)
    end

    # Write JSON
    open(output_path, "w") do f
        JSON3.pretty(f, cytoscape_data)
    end

    println("✓ Cytoscape.js format with $(length(cytoscape_data["elements"]["nodes"])) nodes")
end

function export_csv_format(results_json, output_path::String)
    # Create CSV table with all node data
    rows = []

    for (uuid, activity) in results_json["node_activities"]
        influence = get(results_json["influence_scores"], uuid, 0.0)

        # Check if observed
        obs_activity = ""
        obs_confidence = ""
        if haskey(results_json["observations"], uuid)
            obs_data = results_json["observations"][uuid]
            obs_activity = obs_data["activity"]
            obs_confidence = obs_data["confidence"]
        end

        row = Dict(
            "node_uuid" => uuid,
            "activity" => round(activity, digits=4),
            "influence_score" => round(influence, digits=6),
            "observed" => !isempty(obs_activity),
            "observed_activity" => obs_activity,
            "observed_confidence" => obs_confidence
        )

        push!(rows, row)
    end

    # Convert to DataFrame and save
    df = DataFrame(rows)

    # Sort by influence score (descending)
    sort!(df, :influence_score, rev=true)

    CSV.write(output_path, df)

    println("✓ CSV export with $(nrow(df)) nodes")
end

function execute_server_command(args)
    println("🚀 Starting DeltaSignal server...")
    println("   Host: $(args[:host])")
    println("   Port: $(args[:port])")
    
    try
        # TODO: Implement web server
        println("⚠️  Web server not yet implemented")
        
    catch e
        println("❌ Error starting server: $e")
        exit(1)
    end
end

function main()
    println("🧬 DeltaSignal - Pathway Perturbation & Dynamics Engine")
    println("=" ^ 60)
    
    args = parse_commandline()
    
    if args[:command] == "parse"
        execute_parse_command(args)
    elseif args[:command] == "solve"
        execute_solve_command(args)
    elseif args[:command] == "rollout"
        execute_rollout_command(args)
    elseif args[:command] == "train"
        execute_train_command(args)
    elseif args[:command] == "validate"
        execute_validate_command(args)
    elseif args[:command] == "export"
        execute_export_command(args)
    elseif args[:command] == "server"
        execute_server_command(args)
    end
    
    println("\n✨ Complete!")
end

# Run main function if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
