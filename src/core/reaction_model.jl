# Complete Reaction Model Implementation

"""
Parameters for a single reaction/node in the network.
"""
struct ReactionParams
    # Hill function parameters
    h::Float64          # Hill coefficient (steepness)
    K::Float64          # Half-maximum threshold
    
    # Activator parameters (per input)
    activator_weights::Vector{Float64}      # Weights for geometric mean
    activator_sensitivity_s::Vector{Float64}  # Sensitivity strengths
    activator_sensitivity_n::Vector{Float64}  # Sensitivity cooperativities  
    activator_sensitivity_K::Vector{Float64}  # Sensitivity half-max
    
    # Inhibitor parameters (per input)
    inhibitor_betas::Vector{Float64}        # Inhibition strengths
    inhibitor_ms::Vector{Float64}           # Inhibition Hill coefficients
    
    # Substrate parameters (per input)
    substrate_weights::Vector{Float64}      # Weights for availability
    
    # Time-dynamic parameters
    consumption_lambdas::Vector{Float64}    # Substrate consumption rates
    production_etas::Vector{Float64}        # Product formation rates
    replenishment_rho::Float64             # Baseline replenishment
    decay_delta::Float64                   # Product decay
end

"""
A complete reaction in the network with all inputs and parameters.

`depletion_uuids` is a SEPARATE set of negative inputs marked with
edge_type="depletion" — these represent catalyst→substrate "consumption"
edges (e.g., PTEN → PIP3, MDM2 → TP53). Unlike normal inhibitor edges which
typically use devspec (only suppresses above baseline), depletion edges use
divide-form inhibition so a knocked-out catalyst boosts the substrate
(de-repression). This is what unblocks substrate-depletion biology without
the global -14pp regression that comes from switching all inhibitor edges
to divide form.
"""
struct Reaction
    target_uuid::String
    activator_uuids::Vector{String}
    inhibitor_uuids::Vector{String}
    depletion_uuids::Vector{String}
    substrate_uuids::Vector{String}
    product_uuids::Vector{String}
    params::ReactionParams
    is_and_gate::Bool  # legacy reaction-level flag (true if ANY edge is AND)

    # Per-edge AND/OR flags, aligned with the activator/inhibitor vectors above.
    # Distinguishes AND-clustered inputs (geomean: all required) from OR-clustered
    # inputs (max: any sufficient). Without this, e.g. P53 OR SMAD2 → P21 would
    # be aggregated with geomean and dragged down by whichever input is at baseline.
    activator_is_and::Vector{Bool}
    inhibitor_is_and::Vector{Bool}

    # Per-activator catalyst flag (edge_type=="catalyst"), aligned with
    # activator_uuids. Catalysts are regenerated (conserved moiety), not
    # consumed — used by the SCC-aware solver to relax catalyst edges that
    # close a recycling cycle (see DS_SCC_SOLVE in compute_reaction_output_vec).
    activator_is_catalyst::Vector{Bool}
end

# Back-compat outer constructor: callers predating the per-activator catalyst
# flag (the biological-realism / feedback / compartmentalization enhancement
# modules) pass the 10-arg form. Default to no catalysts so they keep working.
function Reaction(
    target_uuid::String,
    activator_uuids::Vector{String},
    inhibitor_uuids::Vector{String},
    depletion_uuids::Vector{String},
    substrate_uuids::Vector{String},
    product_uuids::Vector{String},
    params::ReactionParams,
    is_and_gate::Bool,
    activator_is_and::Vector{Bool},
    inhibitor_is_and::Vector{Bool},
)
    return Reaction(
        target_uuid, activator_uuids, inhibitor_uuids, depletion_uuids,
        substrate_uuids, product_uuids, params, is_and_gate,
        activator_is_and, inhibitor_is_and, fill(false, length(activator_uuids)),
    )
end

"""
Convert logic network to reaction-based representation.
Groups edges by target node and creates reactions.
"""
function convert_to_reaction_network(network::ReactionNetwork)::Vector{Reaction}
    
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
        reaction = create_reaction_from_edges(target_uuid, edges, network)
        push!(reactions, reaction)
    end
    
    return reactions
end

"""
Create a single reaction from grouped edges targeting the same node.
"""
function create_reaction_from_edges(
    target_uuid::String,
    edges::Vector{LogicNetworkEdge},
    network::ReactionNetwork
)::Reaction
    
    # Separate edges by role, preserving each edge's is_and flag.
    # Depletion edges (catalyst→substrate, marked with edge_type="depletion")
    # are split into their own vector so the propagator can apply divide-form
    # inhibition to them specifically while regular inhibitors stay on
    # whatever DS_INHIBITION_MODE is set to (default devspec).
    activators = String[]
    activator_is_and = Bool[]
    activator_is_catalyst = Bool[]
    inhibitors = String[]
    inhibitor_is_and = Bool[]
    depletions = String[]
    substrates = String[]
    products = String[]

    for edge in edges
        if edge.is_positive
            push!(activators, edge.parent_uuid)
            push!(activator_is_and, edge.is_and)
            push!(activator_is_catalyst, edge.edge_type == "catalyst")
        elseif edge.edge_type == "depletion"
            push!(depletions, edge.parent_uuid)
        else
            push!(inhibitors, edge.parent_uuid)
            push!(inhibitor_is_and, edge.is_and)
        end
    end

    # Reaction-level fallback flag (preserved for any legacy caller).
    is_and_gate = any(edge.is_and for edge in edges)

    params = create_default_reaction_params(
        length(activators),
        length(inhibitors),
        length(substrates),
    )

    return Reaction(
        target_uuid,
        activators,
        inhibitors,
        depletions,
        substrates,
        products,
        params,
        is_and_gate,
        activator_is_and,
        inhibitor_is_and,
        activator_is_catalyst,
    )
end

"""
Create default parameters for a reaction with specified input counts.
"""
function create_default_reaction_params(
    n_activators::Int,
    n_inhibitors::Int,
    n_substrates::Int
)::ReactionParams
    
    # Per spec section 5/15: neutral defaults.
    # K=0.1 anchors near the spec baseline x₀=0.01 such that Hill(x₀, 2, 0.1) ≈ x₀,
    # giving a self-consistent fixed point at baseline.
    h = 2.0
    K = 0.1

    # Activators: uniform weights, neutral sensitivity (s=0 ⇒ α=1)
    activator_weights = n_activators > 0 ? fill(1.0/n_activators, n_activators) : Float64[]
    activator_sensitivity_s = fill(0.0, n_activators)
    activator_sensitivity_n = fill(2.0, n_activators)
    activator_sensitivity_K = fill(0.1, n_activators)

    # Inhibitors: β=0 by default per spec (neutral). Specific networks can
    # raise β to engage inhibition; we don't apply blanket suppression here.
    inhibitor_betas = fill(0.0, n_inhibitors)
    inhibitor_ms = fill(2.0, n_inhibitors)

    # Substrates: uniform weights
    substrate_weights = n_substrates > 0 ? fill(1.0/n_substrates, n_substrates) : Float64[]

    # Time-dynamic defaults (only used in TD mode)
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
Compute reaction output given input node activities.
Implements the complete mathematical model from the specification.
"""
function compute_reaction_output(
    reaction::Reaction,
    node_activities::Dict{String, Float64}
)::Float64
    
    # Get input activities
    activator_activities = [get(node_activities, uuid, 0.0) for uuid in reaction.activator_uuids]
    inhibitor_activities = [get(node_activities, uuid, 0.0) for uuid in reaction.inhibitor_uuids]
    substrate_activities = [get(node_activities, uuid, 0.0) for uuid in reaction.substrate_uuids]
    
    # Step 1: Apply sensitivity transforms to activators
    transformed_activators = Float64[]
    for (i, activity) in enumerate(activator_activities)
        transformed = apply_sensitivity_transform(
            activity;
            s = reaction.params.activator_sensitivity_s[i],
            n = reaction.params.activator_sensitivity_n[i],
            K_α = reaction.params.activator_sensitivity_K[i]
        )
        push!(transformed_activators, transformed)
    end
    
    # Step 2: Compute activator aggregation (A).
    # Split by per-edge AND/OR (see compute_reaction_output_vec for the full
    # rationale). AND-cluster → weighted geomean, OR-cluster → max,
    # combined → max. Empty activator set falls back to target baseline.
    target_baseline = get(node_activities, reaction.target_uuid, 0.01)
    if length(transformed_activators) > 0
        and_vals = Float64[]
        and_wts = Float64[]
        or_vals = Float64[]
        for (i, v) in enumerate(transformed_activators)
            flag = i <= length(reaction.activator_is_and) ? reaction.activator_is_and[i] : reaction.is_and_gate
            if flag
                push!(and_vals, v)
                push!(and_wts, reaction.params.activator_weights[i])
            else
                push!(or_vals, v)
            end
        end
        and_A = isempty(and_vals) ? nothing : geometric_mean_aggregator(and_vals, and_wts)
        or_A = isempty(or_vals) ? nothing : maximum(or_vals)
        A = if and_A === nothing
            or_A
        elseif or_A === nothing
            and_A
        else
            max(and_A, or_A)
        end
    else
        A = target_baseline
    end
    
    # Step 3: Compute inhibition aggregation (H)  
    if length(inhibitor_activities) > 0
        H = inhibition_aggregator(
            inhibitor_activities,
            reaction.params.inhibitor_betas,
            reaction.params.inhibitor_ms
        )
    else
        H = 1.0  # No inhibition
    end
    
    # Step 4: Compute substrate availability (L)
    if length(substrate_activities) > 0
        L = substrate_availability_aggregator(
            substrate_activities,
            reaction.params.substrate_weights
        )
    else
        L = 1.0  # No substrate limitation
    end
    
    # Step 5: Reaction output = A·H·L, clamped to [0,1].
    # We deliberately drop the spec's Section 2.3 output Hill step here.
    # With Hill's default h=2, K=0.1, baseline x₀=0.01 is an UNSTABLE fixed
    # point (F'(x₀) ≈ 1.96 > 1) — there's no K that simultaneously satisfies
    # (a) Hill(baseline)=baseline and (b) F'(baseline) < 1. Without Hill,
    # F=A·H·L is exactly boundary-stable at baseline (slope 1) and amounts
    # propagate proportionally, which is what the user wants for the
    # interactive perturbation case. Hill output can be re-enabled per
    # reaction when h_r and K_r are explicitly trained for that reaction.
    return clamp(A * H * L, zero(typeof(A)), one(typeof(A)))
end

"""
Compute the forward model F(x; θ) for all reactions.
Returns new node activities based on reaction outputs.
"""
function forward_model(
    current_activities::Dict{String, Float64},
    reactions::Vector{Reaction}
)::Dict{String, Float64}

    new_activities = copy(current_activities)

    for reaction in reactions
        output = compute_reaction_output(reaction, current_activities)
        new_activities[reaction.target_uuid] = output
    end

    return new_activities
end

# --- Vector-form forward model (autodiff-friendly) ---
#
# The dict-keyed compute_reaction_output above is convenient for callers
# but doesn't play well with ForwardDiff — dicts allocate per call and
# Dict{String, Float64} can't hold Dual numbers. The vector form below
# is what the Optim+LBFGS solver actually drives. Same math.

"""
Reaction with integer node indices instead of UUIDs. Built once per solve.
`target_baseline` is the resting-state activity of the target node — used as
the implicit source for the activator term when a reaction has no explicit
activator inputs (e.g. a node modulated by inhibitors only, biologically
"constitutively expressed and suppressed by its inhibitors").
"""
struct IndexedReaction
    target_idx::Int
    target_baseline::Float64
    activator_indices::Vector{Int}
    activator_is_and::Vector{Bool}  # per-activator AND/OR flag
    activator_break::Vector{Bool}   # per-activator: relax this edge under
                                    # DS_SCC_SOLVE. True iff it is a catalyst
                                    # edge whose source sits in the SAME SCC as
                                    # the target — i.e. a recycling-cycle
                                    # closure. The catalyst is then read at
                                    # baseline (conserved-moiety modulator) so
                                    # the artifactual SCC dissolves at solve
                                    # time without editing the network.
    inhibitor_indices::Vector{Int}
    inhibitor_in_short_loop::Vector{Bool}  # per-inhibitor: is the inhibitor a
                                            # node that the reaction's target
                                            # can reach in ≤ DS_LOOP_DEPTH
                                            # forward hops? If so, it's part
                                            # of a small negative-feedback loop
                                            # and DS_INHIBITOR_FLOOR applies.
    depletion_indices::Vector{Int}          # catalyst→substrate "consumption"
                                            # inhibitor edges from
                                            # edge_type="depletion". Always
                                            # treated with divide-form
                                            # inhibition (separate from the
                                            # DS_INHIBITION_MODE selection,
                                            # so KO-of-catalyst can boost the
                                            # substrate via de-repression).
    substrate_indices::Vector{Int}
    params::ReactionParams
end

"""
Return a copy of `p` whose edge-aligned parameter vectors keep only the entries
at the surviving original positions (`act_orig`/`inh_orig`/`sub_orig`, in
increasing order). This keeps per-edge parameters aligned with the filtered
activator/inhibitor/substrate index vectors in `index_reactions`. It is a no-op
in value when nothing was dropped (positions == 1:n). Defensive: if a parameter
vector is shorter than the requested positions (e.g. reactions built by a
non-default constructor), that vector is left unchanged rather than erroring.
"""
function compact_reaction_params(
    p::ReactionParams,
    act_orig::Vector{Int},
    inh_orig::Vector{Int},
    sub_orig::Vector{Int},
)::ReactionParams
    pick(v, idx) = (isempty(idx) || maximum(idx) <= length(v)) ? v[idx] : v
    return ReactionParams(
        p.h, p.K,
        pick(p.activator_weights, act_orig),
        pick(p.activator_sensitivity_s, act_orig),
        pick(p.activator_sensitivity_n, act_orig),
        pick(p.activator_sensitivity_K, act_orig),
        pick(p.inhibitor_betas, inh_orig),
        pick(p.inhibitor_ms, inh_orig),
        pick(p.substrate_weights, sub_orig),
        pick(p.consumption_lambdas, sub_orig),
        p.production_etas,
        p.replenishment_rho,
        p.decay_delta,
    )
end

"""
Convert Vector{Reaction} → Vector{IndexedReaction} using a uuid → index map
and a baseline lookup (uuid → baseline). Reactions targeting nodes outside
the map are dropped (defensive — shouldn't happen for well-formed networks).
Preserves the per-edge AND/OR flag for each activator.
"""
function index_reactions(
    reactions::Vector{Reaction},
    uuid_to_idx::Dict{String, Int},
    baselines::Dict{String, Float64},
)::Tuple{Vector{IndexedReaction}, Vector{Int}, Int}
    # First pass: build the bare IndexedReactions (without loop flags) and the
    # forward adjacency for short-loop detection.
    raw = NamedTuple[]
    n_nodes = length(uuid_to_idx)
    fwd_adj = [Set{Int}() for _ in 1:n_nodes]

    for r in reactions
        haskey(uuid_to_idx, r.target_uuid) || continue
        target_idx = uuid_to_idx[r.target_uuid]

        # Track the ORIGINAL edge position of each surviving input so the
        # per-edge parameter vectors (weights, sensitivity, Hill m/β) can be
        # compacted in lockstep. Without this, dropping an input whose UUID is
        # absent from the node set shifts every later input onto the wrong
        # parameter slot — silent under uniform defaults, corrupting once
        # per-edge parameters are trained.
        act_indices = Int[]
        act_is_and = Bool[]
        act_is_catalyst = Bool[]
        act_orig = Int[]
        for (k, uuid) in enumerate(r.activator_uuids)
            haskey(uuid_to_idx, uuid) || continue
            push!(act_indices, uuid_to_idx[uuid])
            push!(act_is_and, k <= length(r.activator_is_and) ? r.activator_is_and[k] : r.is_and_gate)
            push!(act_is_catalyst, k <= length(r.activator_is_catalyst) ? r.activator_is_catalyst[k] : false)
            push!(act_orig, k)
        end
        inh_orig = [k for (k, u) in enumerate(r.inhibitor_uuids) if haskey(uuid_to_idx, u)]
        sub_orig = [k for (k, u) in enumerate(r.substrate_uuids) if haskey(uuid_to_idx, u)]
        inh_indices = [uuid_to_idx[u] for u in r.inhibitor_uuids if haskey(uuid_to_idx, u)]
        dep_indices = [uuid_to_idx[u] for u in r.depletion_uuids if haskey(uuid_to_idx, u)]
        sub_indices = [uuid_to_idx[u] for u in r.substrate_uuids if haskey(uuid_to_idx, u)]

        # Compact the edge-aligned parameter vectors to match the surviving
        # inputs. No-op when nothing was dropped (act_orig == 1:n), so this is
        # behaviour-preserving for well-formed networks.
        params = compact_reaction_params(r.params, act_orig, inh_orig, sub_orig)

        # Forward adjacency: every activator AND inhibitor source points to this
        # reaction's target (both define how the target's value depends on
        # upstream nodes for short-loop reachability).
        for a in act_indices
            push!(fwd_adj[a], target_idx)
        end
        for i in inh_indices
            push!(fwd_adj[i], target_idx)
        end
        for d in dep_indices
            push!(fwd_adj[d], target_idx)
        end

        push!(raw, (
            target_idx=target_idx, baseline=get(baselines, r.target_uuid, 0.01),
            act_indices=act_indices, act_is_and=act_is_and,
            act_is_catalyst=act_is_catalyst,
            inh_indices=inh_indices, dep_indices=dep_indices,
            sub_indices=sub_indices, params=params,
        ))
    end

    # SCC component ids. Computed when either SCC feature is enabled (keeps the
    # default path byte-for-byte unchanged otherwise). comp_id[i] = SCC id of
    # node i, numbered in reverse-topological order. Used for (a) the SCC-
    # condensation solver's processing order [DS_SCC_SOLVE] and (b) marking
    # recycling catalyst back-edges to relax [DS_SCC_BREAK_CATALYST].
    scc_solve = get(ENV, "DS_SCC_SOLVE", "0") != "0"
    break_catalyst = get(ENV, "DS_SCC_BREAK_CATALYST", "0") != "0"
    comp_id, n_comp = (scc_solve || break_catalyst) ?
        tarjan_scc_components(fwd_adj) : (Int[], 0)

    # Second pass: for each (target, inhibitor) pair, BFS forward from the
    # target up to max_depth hops to see if the inhibitor is reachable. If so,
    # they sit in a small negative-feedback loop and the inhibitor edge is
    # eligible for DS_INHIBITOR_FLOOR dampening.
    max_depth = parse(Int, get(ENV, "DS_LOOP_DEPTH", "3"))
    indexed = IndexedReaction[]
    for rec in raw
        in_loop = Bool[]
        for inh in rec.inh_indices
            found = false
            if inh == rec.target_idx
                # Self-loop: target is its own inhibitor. Trivially in a loop.
                found = true
            else
                visited = Set{Int}([rec.target_idx])
                frontier = Set{Int}([rec.target_idx])
                for _ in 1:max_depth
                    next_frontier = Set{Int}()
                    for u in frontier
                        for v in fwd_adj[u]
                            if v == inh
                                found = true; break
                            end
                            if !(v in visited)
                                push!(visited, v); push!(next_frontier, v)
                            end
                        end
                        found && break
                    end
                    found && break
                    isempty(next_frontier) && break
                    frontier = next_frontier
                end
            end
            push!(in_loop, found)
        end

        # Recycling-catalyst back-edges: a catalyst activator whose source is in
        # the same SCC as this reaction's target closes a recycling cycle. Mark
        # it for relaxation under DS_SCC_SOLVE. All-false when SCC solve is off.
        act_break = Bool[]
        for k in eachindex(rec.act_indices)
            brk = false
            if break_catalyst && rec.act_is_catalyst[k]
                src = rec.act_indices[k]
                if comp_id[src] == comp_id[rec.target_idx]
                    brk = true
                end
            end
            push!(act_break, brk)
        end

        push!(indexed, IndexedReaction(
            rec.target_idx, rec.baseline,
            rec.act_indices, rec.act_is_and, act_break,
            rec.inh_indices, in_loop,
            rec.dep_indices,
            rec.sub_indices, rec.params,
        ))
    end
    return indexed, comp_id, n_comp
end

"""
Iterative Tarjan SCC. `adj[u]` is the set of nodes u points to. Returns
`(comp_id, ncomp)` where comp_id[i] = SCC id of node i (1-based) and ncomp is
the number of components. Nodes in the same cycle share an id; acyclic nodes
each get their own. Components are numbered in REVERSE topological order
(sinks get low ids, sources high), so processing in DESCENDING comp_id is a
valid topological order (upstream before downstream). O(V+E).
"""
function tarjan_scc_components(adj::Vector{Set{Int}})::Tuple{Vector{Int},Int}
    n = length(adj)
    index = fill(0, n)
    low = fill(0, n)
    onstack = falses(n)
    comp_id = fill(0, n)
    stack = Int[]
    counter = 0
    ncomp = 0
    # iterative DFS: work stack holds (node, neighbor-iteration-position)
    for s in 1:n
        index[s] != 0 && continue
        work = Tuple{Int,Int}[(s, 1)]
        # neighbors as indexable vectors (Set isn't indexable)
        neigh = Dict{Int,Vector{Int}}()
        neigh[s] = collect(adj[s])
        counter += 1; index[s] = counter; low[s] = counter
        push!(stack, s); onstack[s] = true
        while !isempty(work)
            v, pi = work[end]
            vs = get!(neigh, v) do; collect(adj[v]) end
            advanced = false
            i = pi
            while i <= length(vs)
                w = vs[i]
                if index[w] == 0
                    work[end] = (v, i + 1)
                    counter += 1; index[w] = counter; low[w] = counter
                    push!(stack, w); onstack[w] = true
                    push!(work, (w, 1))
                    advanced = true
                    break
                elseif onstack[w]
                    low[v] = min(low[v], index[w])
                end
                i += 1
            end
            advanced && continue
            pop!(work)
            if !isempty(work)
                p = work[end][1]
                low[p] = min(low[p], low[v])
            end
            if low[v] == index[v]
                ncomp += 1
                while true
                    w = pop!(stack); onstack[w] = false; comp_id[w] = ncomp
                    w == v && break
                end
            end
        end
    end
    return comp_id, ncomp
end

"""
Compute reaction output from a flat activity vector. Same math as the dict
version (sensitivity transform → geomean activators × Hill-suppression
inhibitors × geomean substrates → Hill output), but generic-typed over the
element type so ForwardDiff Duals propagate.
"""
function compute_reaction_output_vec(x::AbstractVector{T}, rxn::IndexedReaction;
                                     supply::Union{Nothing,AbstractVector}=nothing) where {T<:Real}
    p = rxn.params

    # Activators (and optionally inverted inhibitors): split by per-edge AND/OR
    # flag, then combine.
    #
    #   - AND-clustered activators → weighted geomean / min / signed (DS_AND_MODE)
    #   - OR-clustered  activators → max (any sufficient)
    #   - Combined        A        → max(AND_result, OR_result), so an OR alt
    #                                 bypasses the AND group
    #
    # Empty-activator fallback: A = target's baseline (constitutive source —
    # biology: a node with only inhibitor inputs is presumed expressed at
    # baseline and modulated downward).
    #
    # Inhibition handling, selectable via DS_INHIBITION_MODE:
    #   "spec"     — multiplicative Hill suppression H = ∏ 1/(1+β·x^m). With
    #                β=0 default (untrained), H≡1 — i.e. inhibitor edges have
    #                NO effect. Faithful to the spec but means negative-
    #                regulator perturbations don't propagate out of the box.
    #   "inversion" — invert the inhibitor's source state (piecewise so
    #                baseline is preserved) and treat it as an AND-clustered
    #                activator. Continuous analog of the naive baseline's
    #                {0↔2, 1↔1} integer inversion. Default for benchmark.
    inhibition_mode = get(ENV, "DS_INHIBITION_MODE", "spec")
    bl = T(rxn.target_baseline)

    and_vals = T[]
    and_wts = Float64[]
    or_vals = T[]

    # Activator inputs (with per-input sensitivity transform + per-edge AND/OR)
    @inbounds for k in 1:length(rxn.activator_indices)
        i = rxn.activator_indices[k]
        # Recycling-catalyst back-edge (DS_SCC_BREAK_CATALYST): relax this edge
        # so the catalyst acts as a conserved-moiety modulator rather than a
        # loop-carrying signal. Read it from `supply` (its externally-supplied
        # value, frozen before the SCC was iterated) when available, else fall
        # back to baseline. The catalyst's real perturbation still reaches this
        # reaction through its forward (non-recycling) paths.
        src_val = if k <= length(rxn.activator_break) && rxn.activator_break[k]
            supply === nothing ? bl : T(supply[i])
        else
            x[i]
        end
        transformed = apply_sensitivity_transform(
            src_val;
            s=p.activator_sensitivity_s[k],
            n=p.activator_sensitivity_n[k],
            K_α=p.activator_sensitivity_K[k],
        )
        if rxn.activator_is_and[k]
            push!(and_vals, transformed)
            push!(and_wts, p.activator_weights[k])
        else
            push!(or_vals, transformed)
        end
    end

    # Inverted-inhibitor activators (when inhibition_mode == "inversion").
    # Piecewise inversion: f(0)=1, f(baseline)=baseline, f(1)=0, monotone
    # decreasing. Keeps baseline a fixed point of the dynamics.
    if inhibition_mode == "inversion"
        @inbounds for i in rxn.inhibitor_indices
            v = x[i]
            inv_v = if v <= bl
                one(T) - (v / bl) * (one(T) - bl)
            else
                bl * (one(T) - v) / (one(T) - bl)
            end
            push!(and_vals, clamp(inv_v, zero(T), one(T)))
            push!(and_wts, 1.0)
        end
    end

    A = if isempty(and_vals) && isempty(or_vals)
        bl
    else
        and_result = if isempty(and_vals)
            nothing
        else
            mode = get(ENV, "DS_AND_MODE", "geomean")
            if mode == "min"
                minimum(and_vals)
            elseif mode == "signed"
                summed = bl
                for v in and_vals
                    summed += (v - bl)
                end
                clamp(summed, zero(T), one(T))
            elseif mode == "signed_gated"
                # Hybrid: signed-AND propagation gated by a "completeness" factor.
                # gate = min over inputs of (v / baseline) capped at 1. Equal to 1
                # when every input is at-or-above baseline (so signed propagation
                # is unmodified, preserving upregulation cascades). Below 1 when
                # any input is sub-baseline; exactly 0 when any input is
                # knocked out — which forces A=0 (correct "all required" AND
                # semantics). Fixes the wrong-direction failures where signed
                # alone lets one saturated input cancel another's knockout.
                summed = bl
                gate = one(T)
                for v in and_vals
                    summed += (v - bl)
                    ratio = v / bl
                    gate = min(gate, min(ratio, one(T)))
                end
                clamp(summed * gate, zero(T), one(T))
            elseif mode == "multiplicative"
                # Multiplicative fold-change AND:  out = baseline · Π (vᵢ / baseline)
                # Symmetric (above- and below-baseline behave equivalently),
                # saturates at the clamp, composes log-additively under chains.
                # Matches "multiply for AND positives" in MP-BioPath, and gives
                # the user-expected behaviour: 2·2 ⇒ 4, ½·½ ⇒ ¼, 100·100 ⇒ 100.
                # ε keeps the division safe when an input is exactly zero
                # (knockout); the clamp then bounds the result.
                eps_T = T(1e-6)
                prod_T = bl
                for v in and_vals
                    prod_T *= (v + eps_T) / (bl + eps_T)
                end
                clamp(prod_T, zero(T), one(T))
            elseif mode == "hill_sat"
                # Multiplicative AND with smooth boundary saturation.
                #
                # Behaviour goal (per design): essentially exact multiplication
                # of fold-changes through the common operating range, with
                # smooth (autodiff-friendly) sigmoid transitions ONLY at the
                # UI=0 and UI=100 boundaries. Composes for any number of
                # inputs and any cascade depth without per-hop signal decay.
                #
                #     raw  = bl · Π(vᵢ/bl)                     # raw product
                #     hi   = ½(raw + max − √((raw−max)² + ε²)) # smooth-min cap
                #     out  = ½(hi  +  0  + √((hi −  0)² + ε²)) # smooth-max floor
                #
                # The smooth-min/max pair (Hjelmfelt softening) is identity
                # except inside a transition zone of width ~ε at each boundary.
                # With the default ε = 0.001 (DS_HILL_SAT_EPS) the offset at
                # exact boundary is ε/2 — well below benchmark resolution.
                # For raw in [3·ε, max−3·ε] (UI ≈ 0.3 to 99.7) the output is
                # within 1e-7 of pure multiplication.
                #
                # Per-edge parameters (n, K, max) are NOT plumbed into this
                # formula yet — the planned next step. See memory:
                # project-sigmoid-design-intent.
                eps_T = T(1e-6)
                sat_eps = T(parse(Float64, get(ENV, "DS_HILL_SAT_EPS", "0.001")))
                max_internal = one(T)
                prod_fold = one(T)
                for v in and_vals
                    prod_fold *= (v + eps_T) / (bl + eps_T)
                end
                raw = bl * prod_fold
                # Smooth-min with max_internal (caps at 1.0).
                diff_hi = raw - max_internal
                hi = (raw + max_internal -
                      sqrt(diff_hi * diff_hi + sat_eps * sat_eps)) / T(2.0)
                # Smooth-max with 0 (floors at 0).
                lo = (hi + sqrt(hi * hi + sat_eps * sat_eps)) / T(2.0)
                clamp(lo, zero(T), one(T))
            elseif mode == "hill_log"
                # Sigmoid AND in LOG-FOLD space — gives genuinely continuous
                # outputs, not the bimodal saturation that multiplicative or
                # hill_sat (Hjelmfelt) produce. Math:
                #     log_fold = Σ log(vᵢ / baseline)            (sum of log-folds)
                #     log_out  = z_max · tanh(log_fold / z_max)  (sigmoid in log space)
                #     out      = baseline · exp(log_out)
                # where z_max = log(max_fold) = log(100) ≈ 4.605.
                #
                # Behaviour (for two inputs both at UI):
                #   2 × 2  → 3.6 (close to mult's 4; slight smoothing)
                #   3 × 3  → 7.2 (vs mult's 9; ~80%)
                #   10 × 10 → 51 (vs mult's 100, capped; ~half because tanh saturates)
                #   100 × 1 → 33 (vs mult's 100; tanh in log-space smooths)
                #   100 × 100 → 100 (asymptotically saturated)
                #   0.5 × 0.5 → 0.41 (vs mult's 0.25; sigmoid lifts toward baseline)
                #   0 × 100  → 0 (knockout dominates)
                #
                # The cost: 2-30% deviation from multiplication for moderate
                # fold-changes. The benefit: the output landscape is smooth,
                # so cascade outputs SPAN the [0, max] range instead of
                # clustering at {0, 1, max}. This restores meaningful
                # rank-correlation between predicted and experimental
                # fold-changes (GSEA-style ranking benchmark becomes
                # possible).
                #
                # z_max controlled by DS_HILL_LOG_ZMAX (default log(100) = 4.605).
                eps_T = T(1e-6)
                zmax_default = T(log(100.0))
                z_max = T(parse(Float64, get(ENV, "DS_HILL_LOG_ZMAX",
                                              string(Float64(zmax_default)))))
                log_fold = zero(T)
                for v in and_vals
                    log_fold += log((v + eps_T) / (bl + eps_T))
                end
                log_out = z_max * tanh(log_fold / z_max)
                out = bl * exp(log_out)
                clamp(out, zero(T), one(T))
            else
                geometric_mean_aggregator(and_vals, and_wts)
            end
        end
        # OR cluster = biological alternatives ("any one of {A,B,C}").
        #   "max"    — reaction proceeds at the best-available alternative. Clean
        #              logical-OR reading, but masks single-member knockouts:
        #              max(0, x₀, x₀) = x₀ ⇒ no change registered.
        #   "mean"   — each alternative carries weight, so knocking one out
        #              lowers the cluster. Captures KOs that propagate but
        #              over-predicts when paralogs compensate (e.g. RB1 family).
        #   "median" — middle-rank alternative wins. Preserves paralog
        #              compensation (median of (0,1,1) = 1) but registers KOs
        #              when most alternatives drop (median of (0,0,1) = 0). A
        #              "majority wins" semantics — useful when many OR
        #              alternatives are themselves cascaded from the same
        #              upstream gene.
        or_result = if isempty(or_vals)
            nothing
        else
            or_mode = get(ENV, "DS_OR_MODE", "max")
            if or_mode == "mean"
                sum(or_vals) / length(or_vals)
            elseif or_mode == "median"
                # Plain median; ties are fine because clamp keeps us in [0,1].
                sorted = sort(collect(or_vals))
                n = length(sorted)
                isodd(n) ? sorted[(n + 1) ÷ 2] : (sorted[n ÷ 2] + sorted[n ÷ 2 + 1]) / 2
            else
                maximum(or_vals)
            end
        end
        if and_result === nothing
            or_result
        elseif or_result === nothing
            and_result
        else
            max(and_result, or_result)
        end
    end

    # Inhibition formula, selectable via DS_INHIBITION_MODE:
    #
    #   "spec"    — multiplicative Hill ∏ 1/(1+β·xᵐ). Asymmetric: H≈1 at
    #               baseline but drifts BELOW 1 when x is just slightly above,
    #               which combined with signed AND aggregation cascades into
    #               catastrophic network-wide collapse at any β>0. Default β=0
    #               (inert). The "publishable" β>0 setting doesn't work here.
    #
    #   "krep"    — Hill repressor ∏ Kᵐ/(Kᵐ+xᵐ). Symmetric (H=1 at x=0, full
    #               de-repression on knockout) but still has baseline drift
    #               unless K >> baseline. K controlled by DS_INHIBITOR_K.
    #
    #   "divide"  — Smooth MP-BioPath analog: H = (b+ε)/(x+ε) per inhibitor.
    #               EXACTLY 1 at x=baseline (no spurious drift through signed
    #               AND), >1 for knockout (de-repression), <1 for upregulation.
    #               Matches "doubling the inhibitor halves the target." ε
    #               controlled by DS_INHIBITOR_EPS (default 1e-3, small vs
    #               baseline 0.01). The H factor is bounded by the output
    #               clamping later.
    #
    #   "devspec" — Deviation-only spec: 1/(1+β·max(0, x-b)ᵐ). H=1 below or at
    #               baseline (no drift, no de-repression), suppresses above.
    #               β via DS_INHIBITOR_BETA.
    #
    #   "inversion" — Folded into A as an AND-clustered activator above.
    H = if inhibition_mode == "inversion" || isempty(rxn.inhibitor_indices)
        one(T)
    elseif inhibition_mode == "krep"
        K = T(parse(Float64, get(ENV, "DS_INHIBITOR_K", "0.1")))
        result = one(T)
        @inbounds for k in 1:length(rxn.inhibitor_indices)
            x_inh = clamp(x[rxn.inhibitor_indices[k]], zero(T), one(T))
            m = T(p.inhibitor_ms[k])
            Km = K^m
            result *= Km / (Km + x_inh^m)
        end
        clamp(result, zero(T), one(T))
    elseif inhibition_mode == "divide"
        eps_T = T(parse(Float64, get(ENV, "DS_INHIBITOR_EPS", "0.001")))
        # DS_INHIBITOR_FLOOR caps how strongly an inhibitor edge can suppress
        # its target. By default it applies ONLY to inhibitor edges within
        # small negative-feedback loops (where compounding through the loop
        # over-represses). DS_INHIBITOR_FLOOR_SCOPE=all forces it onto every
        # inhibitor (the global variant, for comparison).
        h_floor = T(parse(Float64, get(ENV, "DS_INHIBITOR_FLOOR", "0.0")))
        floor_scope = get(ENV, "DS_INHIBITOR_FLOOR_SCOPE", "loops")
        result = one(T)
        @inbounds for k in 1:length(rxn.inhibitor_indices)
            x_inh = clamp(x[rxn.inhibitor_indices[k]], zero(T), one(T))
            h_k = (bl + eps_T) / (x_inh + eps_T)
            apply_floor = floor_scope == "all" ||
                          (floor_scope == "loops" && rxn.inhibitor_in_short_loop[k])
            if apply_floor
                h_k = max(h_floor, h_k)
            end
            result *= h_k
        end
        # Cap H to prevent catastrophic de-repression from multiple knockouts.
        # The output is also clamped later, but keeping H bounded keeps signed
        # AND propagation well-behaved upstream.
        clamp(result, zero(T), T(10.0))
    elseif inhibition_mode == "devspec"
        beta_env = get(ENV, "DS_INHIBITOR_BETA", "1.0")
        β = T(parse(Float64, beta_env))
        result = one(T)
        @inbounds for k in 1:length(rxn.inhibitor_indices)
            x_inh = clamp(x[rxn.inhibitor_indices[k]], zero(T), one(T))
            dev = max(zero(T), x_inh - bl)
            m = T(p.inhibitor_ms[k])
            result *= one(T) / (one(T) + β * dev^m)
        end
        clamp(result, zero(T), one(T))
    elseif inhibition_mode == "hill_sat"
        # Divide-form inhibition with smooth boundary saturation.
        #
        # Behaviour goal: exact A/B-style suppression in the common operating
        # range (no compression for inhibitors at UI ≈ 0.25 to 10), smooth
        # transitions only at the H=0 and H=H_max boundaries. This matches
        # the AND-side hill_sat philosophy.
        #
        #     H_raw = Π (b+ε)/(xᵢ+ε)                         # raw A/B product
        #     H_hi  = ½(H_raw + H_max − √((H_raw−H_max)² + ε²))  # smooth cap
        #     H_out = ½(H_hi  +  0    + √(H_hi² + ε²))           # smooth floor
        #
        # Outside the transition zones (width ~ε around each boundary),
        # H_out equals H_raw exactly. The H_max default of 10 caps the
        # de-repression boost from inhibitor knockouts at 10× baseline.
        # Use DS_HILL_SAT_H_MAX to tune, DS_HILL_SAT_EPS for transition width.
        #
        # Same per-edge learnable-parameter intent as hill_sat AND. See
        # memory: project-sigmoid-design-intent.
        eps_T = T(parse(Float64, get(ENV, "DS_INHIBITOR_EPS", "0.001")))
        h_max = T(parse(Float64, get(ENV, "DS_HILL_SAT_H_MAX", "10.0")))
        sat_eps = T(parse(Float64, get(ENV, "DS_HILL_SAT_EPS", "0.001")))
        h_raw = one(T)
        @inbounds for k in 1:length(rxn.inhibitor_indices)
            x_inh = clamp(x[rxn.inhibitor_indices[k]], zero(T), one(T))
            h_raw *= (bl + eps_T) / (x_inh + eps_T)
        end
        diff_hi = h_raw - h_max
        h_hi = (h_raw + h_max -
                sqrt(diff_hi * diff_hi + sat_eps * sat_eps)) / T(2.0)
        h_out = (h_hi + sqrt(h_hi * h_hi + sat_eps * sat_eps)) / T(2.0)
        clamp(h_out, zero(T), h_max)
    else  # "spec" — multiplicative Hill ∏ 1/(1+β·xᵐ)
        n_inh = length(rxn.inhibitor_indices)
        inh = Vector{T}(undef, n_inh)
        @inbounds for k in 1:n_inh
            inh[k] = x[rxn.inhibitor_indices[k]]
        end
        beta_env = get(ENV, "DS_INHIBITOR_BETA", "")
        betas = isempty(beta_env) ? p.inhibitor_betas :
                fill(parse(Float64, beta_env), n_inh)
        inhibition_aggregator(inh, betas, p.inhibitor_ms)
    end

    # Substrate availability (soft-AND)
    L = one(T)
    if !isempty(rxn.substrate_indices)
        n_sub = length(rxn.substrate_indices)
        sub = Vector{T}(undef, n_sub)
        @inbounds for k in 1:n_sub
            sub[k] = x[rxn.substrate_indices[k]]
        end
        L = substrate_availability_aggregator(sub, p.substrate_weights)
    end

    # Depletion factor (H_dep): catalyst → substrate "consumption" edges.
    # Always divide-form `(b+ε)/(x+ε)` regardless of DS_INHIBITION_MODE — the
    # POINT is that catalyst-knockout boosts the substrate via de-repression
    # (H > 1 when catalyst < baseline). devspec and similar can't do this.
    # Capped at DS_DEPLETION_H_MAX (default 10) to bound runaway de-repression.
    H_dep = one(T)
    if !isempty(rxn.depletion_indices)
        eps_dep = T(parse(Float64, get(ENV, "DS_INHIBITOR_EPS", "0.001")))
        h_max_dep = T(parse(Float64, get(ENV, "DS_DEPLETION_H_MAX", "10.0")))
        @inbounds for k in 1:length(rxn.depletion_indices)
            x_dep = clamp(x[rxn.depletion_indices[k]], zero(T), one(T))
            H_dep *= (bl + eps_dep) / (x_dep + eps_dep)
        end
        H_dep = clamp(H_dep, zero(T), h_max_dep)
    end

    # Same rationale as the dict-version compute_reaction_output above:
    # drop the spec's output Hill step in favor of clamped linear propagation
    # to keep baseline as a stable fixed point. Hill primitives remain
    # available for explicit per-reaction use after training.
    return clamp(A * H * H_dep * L, zero(T), one(T))
end

"""
Forward model F(x; θ) in vector form. Returns a new vector where each
target node carries the output of its driving reaction; nodes with no
incoming reaction keep their current value.
"""
function forward_model_vec(x::AbstractVector{T}, reactions::Vector{IndexedReaction};
                           supply::Union{Nothing,AbstractVector}=nothing) where {T<:Real}
    y = copy(x)  # preserves element type, including Dual
    for rxn in reactions
        y[rxn.target_idx] = compute_reaction_output_vec(x, rxn; supply=supply)
    end
    return y
end

"""
Compute Jacobian matrix ∂F/∂x for sensitivity analysis.
Returns sparse matrix representation.
"""
function compute_reaction_jacobian(
    activities::Dict{String, Float64},
    reactions::Vector{Reaction}
)::SparseMatrixCSC{Float64, Int}
    
    # Create mapping from UUIDs to indices
    all_uuids = collect(keys(activities))
    uuid_to_idx = Dict(uuid => i for (i, uuid) in enumerate(all_uuids))
    n = length(all_uuids)
    
    # Initialize sparse matrix components
    I = Int[]  # Row indices
    J = Int[]  # Column indices  
    V = Float64[]  # Values
    
    for reaction in reactions
        target_idx = uuid_to_idx[reaction.target_uuid]
        
        # Compute partial derivatives w.r.t. each input
        for uuid in [reaction.activator_uuids; reaction.inhibitor_uuids; reaction.substrate_uuids]
            if haskey(uuid_to_idx, uuid)
                input_idx = uuid_to_idx[uuid]
                
                # Numerical derivative (small perturbation)
                h = 1e-8
                activities_plus = copy(activities)
                activities_plus[uuid] = min(1.0, activities[uuid] + h)
                
                output_plus = compute_reaction_output(reaction, activities_plus)
                output_current = compute_reaction_output(reaction, activities)
                
                derivative = (output_plus - output_current) / h
                
                # Add to sparse matrix if non-zero
                if abs(derivative) > 1e-12
                    push!(I, target_idx)
                    push!(J, input_idx)
                    push!(V, derivative)
                end
            end
        end
    end
    
    # Create sparse matrix
    return sparse(I, J, V, n, n)
end

"""
Time-dynamic update step for substrate consumption and product formation.
"""
function time_dynamic_update!(
    activities::Dict{String, Float64},
    reactions::Vector{Reaction},
    baseline_activities::Dict{String, Float64},
    damping_factor::Float64 = 0.1
)
    
    # First compute all reaction outputs
    reaction_outputs = Dict{String, Float64}()
    for reaction in reactions
        reaction_outputs[reaction.target_uuid] = compute_reaction_output(reaction, activities)
    end
    
    # Apply substrate consumption
    for reaction in reactions
        y_r = reaction_outputs[reaction.target_uuid]
        
        # Substrate consumption
        for (i, substrate_uuid) in enumerate(reaction.substrate_uuids)
            if i <= length(reaction.params.consumption_lambdas)
                λ = reaction.params.consumption_lambdas[i]
                ρ = reaction.params.replenishment_rho
                x0 = get(baseline_activities, substrate_uuid, 0.01)
                
                # Update: x_k = x_k - λ * y_r * x_k + ρ * (x0 - x_k)
                x_current = activities[substrate_uuid]
                x_new = x_current - λ * y_r * x_current + ρ * (x0 - x_current)
                activities[substrate_uuid] = clamp(x_new, 0.0, 1.0)
            end
        end
        
        # Product formation
        for (i, product_uuid) in enumerate(reaction.product_uuids)
            if i <= length(reaction.params.production_etas)
                η = reaction.params.production_etas[i]
                δ = reaction.params.decay_delta
                x0 = get(baseline_activities, product_uuid, 0.01)
                
                # Substrate availability for this reaction
                if length(reaction.substrate_uuids) > 0
                    substrate_activities_vec = [activities[uuid] for uuid in reaction.substrate_uuids]
                    L = substrate_availability_aggregator(
                        substrate_activities_vec,
                        reaction.params.substrate_weights
                    )
                else
                    L = 1.0
                end
                
                # Update: x_p = x_p + η * y_r * L - δ * (x_p - x0)
                x_current = activities[product_uuid]
                x_new = x_current + η * y_r * L - δ * (x_current - x0)
                activities[product_uuid] = clamp(x_new, 0.0, 1.0)
            end
        end
    end
    
    # Apply global damping for stability
    for uuid in keys(activities)
        baseline = get(baseline_activities, uuid, 0.01)
        activities[uuid] = (1.0 - damping_factor) * activities[uuid] + damping_factor * baseline
    end
    
    # Update regular node activities from reaction outputs
    for (target_uuid, output) in reaction_outputs
        activities[target_uuid] = (1.0 - damping_factor) * activities[target_uuid] + damping_factor * output
    end
end