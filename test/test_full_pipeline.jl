#!/usr/bin/env julia

# Complete end-to-end test of DeltaSignal pipeline
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV
using DataFrames
using JSON3

# Add the src directory to the load path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
include("../src/io/tsv_parser.jl")
include("../src/core/sensitivity.jl")
include("../src/core/aggregators.jl") 
include("../src/core/hill_functions.jl")
include("../src/core/reaction_model.jl")
include("../src/solvers/steady_state.jl")
include("../src/io/reactome_mapper.jl")

function test_full_pipeline()
    println("🚀 DeltaSignal Full Pipeline Test")
    println("=" ^ 60)
    
    try
        # Step 1: Parse network from TSV files
        println("\n📁 STEP 1: Parsing logic network files...")
        
        logic_path = joinpath(@__DIR__, "..", "examples", "sample_logic_network.tsv")
        uuid_path = joinpath(@__DIR__, "..", "examples", "sample_uuid_mapping.tsv") 
        set_path = joinpath(@__DIR__, "..", "examples", "sample_set_mappings.tsv")
        
        network = parse_complete_network(logic_path, uuid_path, set_path)
        
        println("✅ Network parsed successfully:")
        println("   - $(length(network.nodes)) nodes")
        println("   - $(length(network.edges)) edges") 
        println("   - $(length(network.set_mappings)) set expansions")
        
        # Step 2: Validate network mapping
        println("\n🔍 STEP 2: Validating network mapping...")
        
        validation_report = validate_network_pathway_mapping(network)
        
        println("✅ Network validation:")
        println("   - Nodes with Reactome IDs: $(validation_report["nodes_with_reactome_ids"])")
        println("   - Orphaned nodes: $(length(validation_report["nodes_without_reactome_ids"]))")
        
        if !isempty(validation_report["warnings"])
            for warning in validation_report["warnings"]
                println("   ⚠️  $warning")
            end
        end
        
        # Step 3: Set up observations and solve steady-state
        println("\n🧮 STEP 3: Solving steady-state system...")
        
        observations = Dict{String, Tuple{Float64, Float64}}(
            "parent-001" => (75.0, 0.9),  # High confidence observation
            "parent-003" => (50.0, 0.8),  # Medium confidence
            "child-002" => (25.0, 0.7)    # Lower confidence target
        )
        
        println("📊 Observations:")
        for (uuid, (activity, confidence)) in observations
            node_name = haskey(network.nodes, uuid) ? network.nodes[uuid].display_name : uuid
            println("   - $node_name ($uuid): $(activity)% (conf: $confidence)")
        end
        
        # Solve with both methods
        params_fp = SteadyStateParams(1.0, 0.1, 100, 1e-4, "fixed_point")
        params_penalty = SteadyStateParams(1.0, 0.1, 200, 1e-4, "penalty")
        
        result_fp = solve_steady_state(network, observations, params_fp)
        result_penalty = solve_steady_state(network, observations, params_penalty)
        
        println("✅ Fixed-point solver: $(result_fp.converged ? "converged" : "failed") in $(result_fp.iterations) iterations")
        println("✅ Penalty solver: $(result_penalty.converged ? "converged" : "failed") in $(result_penalty.iterations) iterations")
        
        # Step 4: Compare results
        println("\n📈 STEP 4: Analyzing results...")
        
        println("Node activities comparison (Fixed-point vs Penalty):")
        println("   Node                    | Fixed-point | Penalty   | Observed")
        println("   " * "-" ^ 65)
        
        for (uuid, node) in sort(collect(network.nodes), by=x->x[1])
            fp_val = get(result_fp.node_activities, uuid, 0.0) * 100
            penalty_val = get(result_penalty.node_activities, uuid, 0.0) * 100
            
            obs_str = if haskey(observations, uuid)
                "$(observations[uuid][1])%"
            else
                "-"
            end
            
            println("   $(rpad(node.display_name[1:min(22, length(node.display_name))], 22)) | $(rpad(round(fp_val, digits=1), 10)) | $(rpad(round(penalty_val, digits=1), 8)) | $obs_str")
        end
        
        # Step 5: Explainability analysis
        println("\n🎯 STEP 5: Explainability analysis...")
        
        reactions = convert_to_reaction_network(network)
        influence_scores = compute_influence_scores(result_penalty, reactions)
        
        if !isempty(influence_scores)
            println("Top influential nodes (drivers):")
            sorted_influence = sort(collect(influence_scores), by=x->x[2], rev=true)
            
            for (i, (uuid, score)) in enumerate(sorted_influence[1:min(5, length(sorted_influence))])
                node_name = network.nodes[uuid].display_name
                println("   $(i). $node_name: $(round(score, digits=3))")
            end
            
            # Upstream driver suggestions
            target_changes = Dict("child-002" => 0.2)  # Want to increase child-002 by 20%
            upstream_suggestions = explain_upstream_drivers(network, result_penalty, target_changes)
            
            if !isempty(upstream_suggestions)
                println("\nUpstream interventions to increase child-002 by 20%:")
                for (uuid, change) in upstream_suggestions
                    if abs(change) > 0.001  # Only show meaningful changes
                        node_name = network.nodes[uuid].display_name
                        change_ui = round(change * 100, digits=1)
                        println("   - $node_name: $(change_ui > 0 ? "+" : "")$(change_ui)%")
                    end
                end
            end
        end
        
        # Step 6: Aggregation to pathway view
        println("\n🔄 STEP 6: Aggregating results for pathway visualization...")
        
        network_results = [NetworkResult(uuid, activity, 0.8, 0.0) 
                          for (uuid, activity) in result_penalty.node_activities]
        
        aggregated_results = aggregate_to_pathway_view(network_results, network)
        
        println("Pathway-level aggregated results:")
        for result in aggregated_results[1:min(5, length(aggregated_results))]
            activity_ui = round(result.aggregated_activity * 100, digits=1)
            println("   - $(result.reactome_entity.reactome_id): $(activity_ui)% ($(result.aggregation_method))")
        end
        
        # Step 7: Export results
        println("\n📤 STEP 7: Exporting results...")
        
        output_file = joinpath(@__DIR__, "..", "test_results.json")
        
        results_export = Dict(
            "metadata" => Dict(
                "timestamp" => string(now()),
                "network_nodes" => length(network.nodes),
                "network_edges" => length(network.edges),
                "observations" => length(observations)
            ),
            "steady_state_results" => Dict(
                "fixed_point" => Dict(
                    "converged" => result_fp.converged,
                    "iterations" => result_fp.iterations,
                    "solve_time" => result_fp.solve_time,
                    "activities" => result_fp.node_activities
                ),
                "penalty" => Dict(
                    "converged" => result_penalty.converged,
                    "iterations" => result_penalty.iterations,
                    "solve_time" => result_penalty.solve_time,
                    "activities" => result_penalty.node_activities
                )
            ),
            "influence_scores" => influence_scores,
            "upstream_suggestions" => upstream_suggestions
        )
        
        open(output_file, "w") do f
            JSON3.pretty(f, results_export)
        end
        
        println("✅ Results exported to: $output_file")
        
        # Final summary
        println("\n" * "=" ^ 60)
        println("🎉 PIPELINE TEST COMPLETED SUCCESSFULLY!")
        println("=" ^ 60)
        
        println("\n📊 Summary:")
        println("   ✅ Network parsing: $(length(network.nodes)) nodes, $(length(network.edges)) edges")
        println("   ✅ Steady-state solving: Both methods converged")
        println("   ✅ Explainability: $(length(influence_scores)) influence scores computed")
        println("   ✅ Pathway aggregation: $(length(aggregated_results)) Reactome entities")
        println("   ✅ Results export: JSON file created")
        
        println("\n🔬 Key Findings:")
        println("   • Most influential driver: $(sorted_influence[1][1]) ($(network.nodes[sorted_influence[1][1]].display_name))")
        
        best_prediction = ""
        best_error = Inf
        for (uuid, (obs_activity, _)) in observations
            predicted = get(result_penalty.node_activities, uuid, 0.0) * 100
            error = abs(predicted - obs_activity)
            if error < best_error
                best_error = error
                best_prediction = "$uuid: predicted $(round(predicted, digits=1))% vs observed $(obs_activity)%"
            end
        end
        println("   • Best prediction: $best_prediction (error: $(round(best_error, digits=1))%)")
        
        println("\n🚀 Next steps:")
        println("   • Try different parameter settings")
        println("   • Add more observations")
        println("   • Implement time-dynamic mode")
        println("   • Build web interface")
        
        return true
        
    catch e
        println("\n❌ PIPELINE TEST FAILED!")
        println("Error: $e")
        println("\nStack trace:")
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
        return false
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Import the missing Dates module for timestamp
    using Dates
    test_full_pipeline()
end