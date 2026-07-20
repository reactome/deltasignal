# Enhanced Feedback Loop Handling

"""
Detect feedback loops in the network and apply specialized parameter sets.
This allows different mathematical behaviors for different network motifs.
"""
function detect_feedback_loops(network::ReactionNetwork)::Dict{String, Vector{String}}
    # Simple cycle detection using depth-first search
    feedback_loops = Dict{String, Vector{String}}()
    
    # Build adjacency list: source -> [(target, is_positive), ...]
    adjacency = Dict{String, Vector{Tuple{String, Bool}}}()
    
    for edge in network.edges
        if !haskey(adjacency, edge.parent_uuid)
            adjacency[edge.parent_uuid] = []
        end
        push!(adjacency[edge.parent_uuid], (edge.child_uuid, edge.is_positive))
    end
    
    # Find cycles using DFS
    visited = Set{String}()
    rec_stack = Set{String}()
    
    function dfs_cycles(node::String, path::Vector{String}, path_edges::Vector{Bool})
        if node in rec_stack
            # Found a cycle
            cycle_start = findfirst(x -> x == node, path)
            if cycle_start !== nothing
                cycle = path[cycle_start:end]
                cycle_edge_signs = path_edges[cycle_start:end-1]  # Exclude last duplicate
                
                # Determine if positive or negative feedback
                # A feedback loop is negative if it contains an odd number of negative edges
                negative_edges = count(x -> !x, cycle_edge_signs)
                is_negative_feedback = (negative_edges % 2 == 1)
                
                # Add the closing edge back to start
                if haskey(adjacency, cycle[end])
                    for (target, is_positive) in adjacency[cycle[end]]
                        if target == cycle[1]
                            if !is_positive
                                negative_edges += 1
                            end
                            break
                        end
                    end
                end
                
                # Recompute with closing edge
                is_negative_feedback = (negative_edges % 2 == 1)
                
                loop_type = is_negative_feedback ? "negative" : "positive"
                if !haskey(feedback_loops, loop_type)
                    feedback_loops[loop_type] = []
                end
                
                # Create cycle display with closing edge
                full_cycle = vcat(cycle, [cycle[1]])
                push!(feedback_loops[loop_type], join(full_cycle, "→"))
                
                # println("  🔍 Found $loop_type feedback: $(join(full_cycle, \"→\")) (neg_edges: $negative_edges)")
            end
            return
        end
        
        if node in visited
            return
        end
        
        visited = union(visited, [node])
        rec_stack = union(rec_stack, [node])
        push!(path, node)
        
        # Visit neighbors
        if haskey(adjacency, node)
            for (neighbor, is_positive) in adjacency[node]
                new_path_edges = vcat(path_edges, [is_positive])
                dfs_cycles(neighbor, copy(path), new_path_edges)
            end
        end
        
        rec_stack = setdiff(rec_stack, [node])
        pop!(path)
    end
    
    # Start DFS from all nodes
    for node in keys(network.nodes)
        if !(node in visited)
            dfs_cycles(node, String[], Bool[])
        end
    end
    
    return feedback_loops
end

"""
Create specialized reaction parameters for negative feedback loops.
These parameters emphasize strong inhibition and stability.
"""
function create_negative_feedback_params(
    n_activators::Int,
    n_inhibitors::Int,
    n_substrates::Int
)::ReactionParams
    
    # More conservative Hill parameters for stability
    h = 2.0      # Moderate steepness to avoid sharp transitions
    K = 0.4      # Higher threshold for more gradual activation
    
    # Standard activator parameters
    activator_weights = n_activators > 0 ? fill(1.0/n_activators, n_activators) : Float64[]
    activator_sensitivity_s = fill(0.0, n_activators)  # Neutral sensitivity
    activator_sensitivity_n = fill(2.0, n_activators)
    activator_sensitivity_K = fill(0.1, n_activators)
    
    # VERY STRONG inhibitor parameters for effective negative feedback
    inhibitor_betas = fill(n_inhibitors > 0 ? 8.0 : 0.0, n_inhibitors)  # Very strong inhibition
    inhibitor_ms = fill(3.0, n_inhibitors)  # Steep inhibition curve
    
    # Standard substrate parameters
    substrate_weights = n_substrates > 0 ? fill(1.0/n_substrates, n_substrates) : Float64[]
    
    # Minimal time-dynamic effects
    consumption_lambdas = fill(0.05, n_substrates)
    production_etas = Float64[]
    replenishment_rho = 0.005
    decay_delta = 0.005
    
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
Create specialized reaction parameters for positive feedback loops.
These parameters allow amplification while maintaining stability.
"""
function create_positive_feedback_params(
    n_activators::Int,
    n_inhibitors::Int,
    n_substrates::Int
)::ReactionParams
    
    # Slightly more sensitive Hill parameters for amplification
    h = 1.8      # Slightly lower steepness for gradual amplification
    K = 0.25     # Lower threshold for easier activation
    
    # Standard activator parameters
    activator_weights = n_activators > 0 ? fill(1.0/n_activators, n_activators) : Float64[]
    activator_sensitivity_s = fill(0.0, n_activators)  # Neutral sensitivity
    activator_sensitivity_n = fill(2.0, n_activators)
    activator_sensitivity_K = fill(0.1, n_activators)
    
    # Moderate inhibitor parameters
    inhibitor_betas = fill(n_inhibitors > 0 ? 4.0 : 0.0, n_inhibitors)
    inhibitor_ms = fill(2.0, n_inhibitors)
    
    # Standard substrate parameters
    substrate_weights = n_substrates > 0 ? fill(1.0/n_substrates, n_substrates) : Float64[]
    
    # Minimal time-dynamic effects
    consumption_lambdas = fill(0.05, n_substrates)
    production_etas = Float64[]
    replenishment_rho = 0.005
    decay_delta = 0.005
    
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
Enhanced convert_to_reaction_network that applies feedback-specific parameters.
"""
function convert_to_reaction_network_with_feedback_enhancement(network::ReactionNetwork)::Vector{Reaction}
    
    # First detect feedback loops
    feedback_loops = detect_feedback_loops(network)
    
    println("🔍 Detected feedback loops:")
    for (loop_type, loops) in feedback_loops
        println("  $loop_type feedback: $(length(loops)) loops")
        for loop in loops[1:min(3, length(loops))]  # Show first few
            println("    - $loop")
        end
    end
    
    # Identify nodes involved in negative feedback
    negative_feedback_nodes = Set{String}()
    if haskey(feedback_loops, "negative")
        for loop_str in feedback_loops["negative"]
            for node in split(loop_str, "→")
                push!(negative_feedback_nodes, strip(node))
            end
        end
    end
    
    # Identify nodes involved in positive feedback
    positive_feedback_nodes = Set{String}()
    if haskey(feedback_loops, "positive")
        for loop_str in feedback_loops["positive"]
            for node in split(loop_str, "→")
                push!(positive_feedback_nodes, strip(node))
            end
        end
    end
    
    # Group edges by target (child) node
    target_groups = Dict{String, Vector{LogicNetworkEdge}}()
    
    for edge in network.edges
        target = edge.child_uuid
        if !haskey(target_groups, target)
            target_groups[target] = LogicNetworkEdge[]
        end
        push!(target_groups[target], edge)
    end
    
    reactions = Reaction[]
    
    for (target_uuid, edges) in target_groups
        # Determine context-specific parameters
        params = if target_uuid in negative_feedback_nodes
            # Node is part of negative feedback - use strong inhibition
            n_activators = count(e -> e.is_positive, edges)
            n_inhibitors = count(e -> !e.is_positive, edges)
            create_negative_feedback_params(n_activators, n_inhibitors, 0)
        elseif target_uuid in positive_feedback_nodes
            # Node is part of positive feedback - use moderate parameters
            n_activators = count(e -> e.is_positive, edges)
            n_inhibitors = count(e -> !e.is_positive, edges)
            create_positive_feedback_params(n_activators, n_inhibitors, 0)
        else
            # Regular node - use default parameters
            n_activators = count(e -> e.is_positive, edges)
            n_inhibitors = count(e -> !e.is_positive, edges)
            create_default_reaction_params(n_activators, n_inhibitors, 0)
        end
        
        reaction = create_reaction_from_edges_with_params(target_uuid, edges, network, params)
        push!(reactions, reaction)
    end
    
    return reactions
end

"""
Create a reaction from edges with specific parameters.
"""
function create_reaction_from_edges_with_params(
    target_uuid::String,
    edges::Vector{LogicNetworkEdge},
    network::ReactionNetwork,
    params::ReactionParams
)::Reaction
    
    # Separate edges by type
    activators = String[]
    inhibitors = String[]
    substrates = String[]
    products = String[]
    
    # Determine gate type
    is_and_gate = any(edge.is_and for edge in edges)
    
    for edge in edges
        if edge.is_positive
            push!(activators, edge.parent_uuid)
        else
            push!(inhibitors, edge.parent_uuid)
        end
    end
    
    return Reaction(
        target_uuid,
        activators,
        inhibitors,
        String[],   # depletion_uuids: none for feedback-enhanced reactions
        substrates,
        products,
        params,
        is_and_gate,
        fill(is_and_gate, length(activators)),  # activator_is_and
        fill(is_and_gate, length(inhibitors)),  # inhibitor_is_and
    )
end