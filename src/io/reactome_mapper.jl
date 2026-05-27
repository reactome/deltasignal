# Reactome Mapping and Set Expansion Handling

"""
Map network results back to original Reactome entities for PathwayBrowser visualization.
Handles the many-to-many mapping between expanded network nodes and original pathway entities.
"""

struct ReactomeEntity
    reactome_id::String
    entity_type::String
    display_name::String
    pathway_coordinates::Union{Dict{String, Any}, Nothing}
end

struct NetworkResult
    node_uuid::String
    activity::Float64
    confidence::Float64
    upstream_influence::Float64
end

struct AggregatedResult
    reactome_entity::ReactomeEntity
    aggregated_activity::Float64
    aggregated_confidence::Float64
    contributing_nodes::Vector{String}
    aggregation_method::String
end

"""
Aggregate network results back to original Reactome entities.
Multiple aggregation strategies are supported based on biological context.
"""
function aggregate_to_pathway_view(
    network_results::Vector{NetworkResult},
    reaction_network::ReactionNetwork;
    aggregation_method::String = "stoichiometry_weighted"
)::Vector{AggregatedResult}
    
    # Group results by Reactome entity
    reactome_groups = Dict{String, Vector{NetworkResult}}()
    reactome_entities = Dict{String, ReactomeEntity}()
    
    for result in network_results
        node = reaction_network.nodes[result.node_uuid]
        
        if node.reactome_id !== nothing
            reactome_id = node.reactome_id
            
            if !haskey(reactome_groups, reactome_id)
                reactome_groups[reactome_id] = NetworkResult[]
                reactome_entities[reactome_id] = ReactomeEntity(
                    reactome_id,
                    node.entity_type,
                    node.display_name,
                    nothing  # Coordinates to be filled from pathway data
                )
            end
            
            push!(reactome_groups[reactome_id], result)
        end
    end
    
    # Aggregate each group
    aggregated_results = AggregatedResult[]
    
    for (reactome_id, group_results) in reactome_groups
        if length(group_results) == 1
            # Single node - no aggregation needed
            result = group_results[1]
            aggregated_result = AggregatedResult(
                reactome_entities[reactome_id],
                result.activity,
                result.confidence,
                [result.node_uuid],
                "single"
            )
        else
            # Multiple nodes - apply aggregation strategy
            aggregated_result = apply_aggregation_strategy(
                reactome_entities[reactome_id],
                group_results,
                reaction_network,
                aggregation_method
            )
        end
        
        push!(aggregated_results, aggregated_result)
    end
    
    return aggregated_results
end

"""
Apply specific aggregation strategy to combine results from multiple network nodes.
"""
function apply_aggregation_strategy(
    entity::ReactomeEntity,
    results::Vector{NetworkResult},
    network::ReactionNetwork,
    method::String
)::AggregatedResult
    
    activities = [r.activity for r in results]
    confidences = [r.confidence for r in results]
    node_uuids = [r.node_uuid for r in results]
    
    if method == "mean"
        # Simple arithmetic mean
        agg_activity = mean(activities)
        agg_confidence = mean(confidences)
        
    elseif method == "max"
        # Maximum activity (dominant member)
        max_idx = argmax(activities)
        agg_activity = activities[max_idx]
        agg_confidence = confidences[max_idx]
        
    elseif method == "min"
        # Minimum activity (limiting member for complexes)
        min_idx = argmin(activities)
        agg_activity = activities[min_idx]
        agg_confidence = confidences[min_idx]
        
    elseif method == "stoichiometry_weighted"
        # Weight by stoichiometry if available
        weights = get_stoichiometry_weights(node_uuids, network)
        agg_activity = sum(activities .* weights) / sum(weights)
        agg_confidence = sum(confidences .* weights) / sum(weights)
        
    elseif method == "confidence_weighted"
        # Weight by confidence
        weight_sum = sum(confidences)
        if weight_sum > 0
            agg_activity = sum(activities .* confidences) / weight_sum
            agg_confidence = weight_sum / length(confidences)  # Average confidence
        else
            agg_activity = mean(activities)
            agg_confidence = 0.0
        end
        
    elseif method == "geometric_mean"
        # Geometric mean (good for multiplicative effects)
        agg_activity = exp(mean(log.(max.(activities, 1e-12))))
        agg_confidence = exp(mean(log.(max.(confidences, 1e-12))))
        
    else
        error("Unknown aggregation method: $method")
    end
    
    return AggregatedResult(
        entity,
        clamp(agg_activity, 0.0, 1.0),
        clamp(agg_confidence, 0.0, 1.0),
        node_uuids,
        method
    )
end

"""
Get stoichiometry weights for nodes based on their roles in reactions.
"""
function get_stoichiometry_weights(node_uuids::Vector{String}, network::ReactionNetwork)::Vector{Float64}
    weights = ones(Float64, length(node_uuids))
    
    # Extract stoichiometry information from edges
    for (i, uuid) in enumerate(node_uuids)
        # Find edges involving this node and sum stoichiometry
        total_stoich = 0.0
        count = 0
        
        for edge in network.edges
            if edge.parent_uuid == uuid || edge.child_uuid == uuid
                total_stoich += edge.stoichiometry
                count += 1
            end
        end
        
        if count > 0
            weights[i] = total_stoich / count
        end
    end
    
    # Normalize weights
    weight_sum = sum(weights)
    if weight_sum > 0
        weights ./= weight_sum
    else
        weights .= 1.0 / length(weights)
    end
    
    return weights
end

"""
Create mapping from original sets to their expanded members.
This helps identify which network nodes came from the same original complex/set.
"""
function create_set_member_mapping(network::ReactionNetwork)::Dict{String, Vector{String}}
    set_to_members = Dict{String, Vector{String}}()
    
    for (uuid, node) in network.nodes
        if node.original_set_id !== nothing
            set_id = node.original_set_id
            
            if !haskey(set_to_members, set_id)
                set_to_members[set_id] = String[]
            end
            
            push!(set_to_members[set_id], uuid)
        end
    end
    
    return set_to_members
end

"""
Get all network nodes that correspond to a given Reactome entity.
Handles both direct mappings and set expansions.
"""
function get_network_nodes_for_reactome_id(
    reactome_id::String,
    network::ReactionNetwork
)::Vector{NetworkNode}
    
    matching_nodes = NetworkNode[]
    
    for (uuid, node) in network.nodes
        if node.reactome_id == reactome_id
            push!(matching_nodes, node)
        end
    end
    
    return matching_nodes
end

"""
Export aggregated results in format suitable for PathwayBrowser overlay.
"""
function export_pathway_browser_overlay(
    aggregated_results::Vector{AggregatedResult},
    output_path::String
)
    
    # Create overlay data structure for PathwayBrowser
    overlay_data = Dict{String, Any}(
        "type" => "activity_overlay",
        "entities" => Dict{String, Any}[]
    )
    
    for result in aggregated_results
        entity_data = Dict{String, Any}(
            "id" => result.reactome_entity.reactome_id,
            "activity" => result.aggregated_activity * 100,  # Convert to 0-100 scale
            "confidence" => result.aggregated_confidence,
            "type" => result.reactome_entity.entity_type,
            "aggregation_method" => result.aggregation_method,
            "contributing_nodes" => length(result.contributing_nodes)
        )
        
        push!(overlay_data["entities"], entity_data)
    end
    
    # Write JSON file
    open(output_path, "w") do f
        JSON3.pretty(f, overlay_data)
    end
    
    println("PathwayBrowser overlay exported to: $output_path")
end

"""
Validate the mapping between network and pathway entities.
Report potential issues with the aggregation process.
"""
function validate_network_pathway_mapping(network::ReactionNetwork)::Dict{String, Any}
    
    validation_report = Dict{String, Any}(
        "total_network_nodes" => length(network.nodes),
        "nodes_with_reactome_ids" => 0,
        "nodes_without_reactome_ids" => String[],
        "reactome_id_frequency" => Dict{String, Int}(),
        "set_expansion_stats" => Dict{String, Any}(),
        "warnings" => String[]
    )
    
    # Count mappings and identify issues
    reactome_id_counts = Dict{String, Int}()
    
    for (uuid, node) in network.nodes
        if node.reactome_id !== nothing
            validation_report["nodes_with_reactome_ids"] += 1
            
            # Count frequency of each Reactome ID
            if !haskey(reactome_id_counts, node.reactome_id)
                reactome_id_counts[node.reactome_id] = 0
            end
            reactome_id_counts[node.reactome_id] += 1
        else
            push!(validation_report["nodes_without_reactome_ids"], uuid)
        end
    end
    
    validation_report["reactome_id_frequency"] = reactome_id_counts
    
    # Analyze set expansions
    set_member_mapping = create_set_member_mapping(network)
    validation_report["set_expansion_stats"] = Dict{String, Any}(
        "total_sets" => length(set_member_mapping),
        "expansion_sizes" => [length(members) for members in values(set_member_mapping)]
    )
    
    # Generate warnings
    if validation_report["nodes_with_reactome_ids"] == 0
        push!(validation_report["warnings"], "No network nodes have Reactome IDs - pathway visualization will be empty")
    end
    
    orphaned_nodes = length(validation_report["nodes_without_reactome_ids"])
    if orphaned_nodes > 0
        push!(validation_report["warnings"], "$orphaned_nodes network nodes lack Reactome IDs and will be excluded from pathway view")
    end
    
    # Check for highly expanded entities (potential issues)
    for (reactome_id, count) in reactome_id_counts
        if count > 10
            push!(validation_report["warnings"], "Reactome ID $reactome_id maps to $count network nodes - verify aggregation method")
        end
    end
    
    return validation_report
end