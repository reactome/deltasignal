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

    # Per-activator assembly flag (edge_type=="assembly"), aligned with
    # activator_uuids. Assembly edges are the synthetic member→complex edges of
    # a decomposed complex: a Complex is the AND of its subunits and cannot
    # exceed the abundance of its scarcest subunit (hard stoichiometric cap).
    # Used by DS_ASSEMBLY_LIMITING to aggregate assembly inputs with a
    # limiting-reactant (min) rule instead of the permissive DS_AND_MODE, so
    # overexpressing one subunit of a many-membered complex doesn't spuriously
    # drive the whole complex up (the dominant curator OE false-positive mode).
    activator_is_assembly::Vector{Bool}

    # Per-activator GROUP id for `composition` edges, aligned with
    # activator_uuids; 0 = ungrouped. Inputs sharing a group are node copies of
    # ONE entity (the same base stable id -- HDR's BCDX2 complex exists as 33
    # copies, one per variant reaction). Under DS_COMPOSITION_GROUP they are
    # aggregated as alternatives (max within the group) before the
    # limiting-reactant min across distinct components. Off, they are min'd
    # like any assembly input, so the container is capped by the LEAST copy.
    activator_group::Vector{Int}
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
        activator_is_and, inhibitor_is_and,
        fill(false, length(activator_uuids)),   # activator_is_catalyst
        fill(false, length(activator_uuids)),   # activator_is_assembly
        fill(0, length(activator_uuids)),       # activator_group
    )
end

"""
Convert logic network to reaction-based representation.
Groups edges by target node and creates reactions.
"""
function convert_to_reaction_network(network::ReactionNetwork)::Vector{Reaction}

    # Pass-through suppression (DS_DROP_PASSTHROUGH): drop an `output` edge R→E
    # when the reverse input edge E→R also exists — i.e. reaction R both consumes
    # and re-emits the same entity E (a stable participant / scaffold, e.g. active
    # p53 threaded through ~90 reactions). Such a "producer" only recycles E and
    # otherwise pins it at baseline via OR-max, masking upstream knockouts (the
    # dominant propagator_missed failure). Dropping it as a producer lets E track
    # its NET producers so a knockout can lower it. E is still CONSUMED (its input
    # edge stays) and R's other outputs are untouched. Faithful: R does not
    # produce E de novo. See memory project_loop_taxonomy_finding (Type II).
    drop_pt = _bool_env("DS_DROP_PASSTHROUGH", false)
    edge_pairs = drop_pt ?
        Set{Tuple{String,String}}((e.parent_uuid, e.child_uuid) for e in network.edges) :
        Set{Tuple{String,String}}()

    # Group edges by target (child) node
    target_groups = Dict{String, Vector{LogicNetworkEdge}}()

    for edge in network.edges
        if drop_pt && edge.edge_type == "output" &&
           (edge.child_uuid, edge.parent_uuid) in edge_pairs
            continue  # R→E is a pass-through recycle of E; not a net producer
        end
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
    activator_is_assembly = Bool[]
    activator_group = Int[]
    group_ids = Dict{String, Int}()   # base stable id -> dense group id
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
            # `composition` (complex -> complex that contains it, emitted by LNG
            # along Reactome's hasComponent hierarchy) shares assembly's
            # limiting-reactant semantics: a container cannot exceed the
            # component it is built from. specs/016. This is a DELIBERATE default
            # change for networks that carry the (new, LNG_COMPOSITION_EDGES=1)
            # edge type: before, an unknown edge_type was a plain AND input. No
            # shipped catalog contains it, so shipped results are unchanged.
            push!(activator_is_assembly, edge.edge_type in ("assembly", "composition"))
            # Group composition inputs by the source's base stable id, so node
            # copies of one entity can be aggregated as alternatives.
            g = 0
            if edge.edge_type == "composition"
                node = get(network.nodes, edge.parent_uuid, nothing)
                rid = node === nothing ? nothing : node.reactome_id
                # Key by base stable id so node copies of one entity share a
                # group; fall back to the uuid so EVERY composition edge has
                # g > 0 -- the propagator reads "g > 0" as "is composition".
                base = rid === nothing ? edge.parent_uuid : first(split(rid, "::variant::"))
                g = get!(group_ids, base, length(group_ids) + 1)
            end
            push!(activator_group, g)
        elseif edge.edge_type == "depletion"
            push!(depletions, edge.parent_uuid)
        else
            push!(inhibitors, edge.parent_uuid)
            push!(inhibitor_is_and, edge.is_and)
        end
    end

    # Collapse duplicate parallel activator edges from the SAME source
    # (DS_DEDUP_ACTIVATORS). Reactome routinely curates one entity as both the
    # `input` and the `catalyst` of a reaction, which arrives here as two
    # separate pos/and edges with the same parent_uuid. Each edge takes its own
    # slot in the activator vectors, so the aggregators count that entity twice
    # — under hill_log its log-fold is summed twice, i.e. its fold-change is
    # SQUARED (an entity at UI 10 drives the target 7.7x higher than a single
    # edge would). 4.84% of catalog edges sit in such duplicate groups, across
    # 78 of 92 pathways, so this silently skews every propagated result.
    #
    # One entity should contribute one factor. Keep the first slot and merge the
    # role flags (catalyst/assembly are unioned so the SCC catalyst-break and
    # assembly-limiting layers still see the role). Params are built after this,
    # so per-edge weights stay aligned with the deduplicated vectors.
    if _bool_env("DS_DEDUP_ACTIVATORS", false) && length(activators) > 1
        seen = Dict{String, Int}()
        d_act = String[]; d_and = Bool[]; d_cat = Bool[]; d_asm = Bool[]; d_grp = Int[]
        for k in eachindex(activators)
            src = activators[k]
            j = get(seen, src, 0)
            if j == 0
                push!(d_act, src)
                push!(d_and, activator_is_and[k])
                push!(d_cat, activator_is_catalyst[k])
                push!(d_asm, activator_is_assembly[k])
                push!(d_grp, activator_group[k])
                seen[src] = length(d_act)
            else
                # Merge deterministically so the result cannot depend on edge
                # order in the input file (the repo asserts edge-order
                # invariance). AND wins for is_and: the slot's cluster
                # membership decides whether it lands in and_vals or or_vals,
                # which DS_OR_COMBINE=gate makes load-bearing, so a
                # first-occurrence tie-break would let row order change the
                # propagated value. Catalyst unions (it only marks SCC
                # break-eligibility). Assembly requires ALL duplicates to be
                # assembly, so a mixed input+assembly pair stays a normal
                # activator contributing its own fold factor rather than being
                # reclassified into the assembly_min limiting rule.
                d_and[j] |= activator_is_and[k]
                d_cat[j] |= activator_is_catalyst[k]
                d_asm[j] &= activator_is_assembly[k]
                # A composition edge parallel to a plain activator from the same
                # source repeats a fold the plain edge already carries, so the
                # plain role wins: the group survives only if EVERY duplicate is
                # composition. `max` let the group id (read as "is composition"
                # under DS_COMPOSITION_MODE=limit*) delete the input role and turn
                # a 4x input into a <=1 limiter.
                d_grp[j] = (d_grp[j] > 0 && activator_group[k] > 0) ? max(d_grp[j], activator_group[k]) : 0
            end
        end
        activators, activator_is_and = d_act, d_and
        activator_is_catalyst, activator_is_assembly = d_cat, d_asm
        activator_group = d_grp
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
        activator_is_assembly,
        activator_group,
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

# --- Vector-form model (autodiff-friendly) ---
# compute_reaction_output_vec below is the ONLY propagator the solver drives.
# The dict-keyed compute_reaction_output / forward_model variants were removed
# (2026-07): they read none of the DS_* config and had diverged from this one.

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
    activator_in_loop::Vector{Bool}  # per-activator: does this edge CLOSE a
                                    # cycle -- source in the same SCC as the
                                    # target? Under DS_LOOP_ELASTICITY < 1 the
                                    # input fold is read through fold^eps on
                                    # exactly these edges, so a positive loop's
                                    # gain at baseline drops below 1 and
                                    # baseline becomes a stable state instead
                                    # of a knife-edge. Acyclic edges untouched.
    activator_is_assembly::Vector{Bool}  # per-activator: edge_type=="assembly"
    activator_comp_redundant::Vector{Bool} # per-activator: a composition edge
                                    # whose source ALREADY feeds a producing
                                    # reaction of this target (S -> R -> T).
                                    # 54% of the catalog's composition edges.
                                    # Under DS_COMPOSITION_MODE=limit_novel
                                    # these are skipped: the producing reaction
                                    # already carries the component's fold, and
                                    # multiplying it in again squares it at
                                    # every level of a nested hierarchy.
    activator_group::Vector{Int}     # per-activator composition group; > 0 iff
                                    # the edge is `composition`, copies of one
                                    # entity share an id. See DS_COMPOSITION_GROUP
                                    # and DS_COMPOSITION_MODE.
                                    # (member→complex). Under DS_ASSEMBLY_LIMITING
                                    # these inputs are aggregated with a
                                    # limiting-reactant (min) rule — a complex
                                    # cannot exceed its scarcest subunit — rather
                                    # than the permissive DS_AND_MODE.
    inhibitor_indices::Vector{Int}
    inhibitor_is_and::Vector{Bool}  # per-inhibitor AND/OR flag, mirroring
                                    # activator_is_and. AND-clustered inhibitors
                                    # repress cooperatively (their suppression
                                    # factors multiply); OR-clustered ones are
                                    # alternative repressors where any one
                                    # suffices, so the dominant (most
                                    # suppressing) one applies instead of the
                                    # product. Honored under DS_INHIBITOR_OR.
    inhibitor_in_short_loop::Vector{Bool}  # per-inhibitor: is the inhibitor a
                                            # node that the reaction's target
                                            # can reach in ≤ DS_LOOP_DEPTH
                                            # forward hops (or shares its SCC)?
                                            # If so, it's part of a negative-
                                            # feedback loop and DS_INHIBITOR_FLOOR
                                            # (scope=loops) applies.
    inhibitor_transcriptional::Vector{Bool} # per-inhibitor: is this a
                                            # TRANSCRIPTIONAL autoregulation edge
                                            # — the reaction has a gene input
                                            # (it's a transcription/expression
                                            # step) AND the inhibitor is in the
                                            # target's SCC (self-regulating
                                            # feedback). These get weakened
                                            # (DS_INHIBITOR_FLOOR scope=transcription)
                                            # because continuous gene dosage +
                                            # strong proportional repression is
                                            # unphysical; protein-level feedback
                                            # (no gene input) keeps full strength.
    inhibitor_shared::Vector{Vector{Int}}   # per-inhibitor: node indices of this
                                            # reaction's activators that the
                                            # inhibitor CONTAINS (specs/022).
                                            # Empty unless DS_SELF_INHIBITOR_WEIGHT
                                            # is set, so the default never reads it.
    depletion_indices::Vector{Int}          # catalyst→substrate "consumption"
    depletion_break::Vector{Bool}           # per-depletion: a recycling CLOSURE
                                            # (source and target in one SCC) under
                                            # DS_SCC_BREAK_ROLES=depletion; read at
                                            # the component-entry value (`supply`).
    depletion_own_product::Vector{Bool}     # per-depletion: is the depleter a
                                            # direct product of a reaction that
                                            # consumes this target (X -> R -> P,
                                            # P -| X), or a direct successor?
                                            # Under DS_DEPLETION_OWN_PRODUCT=
                                            # suppress_only such an edge may
                                            # suppress but not de-repress.
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

# Gene stable_ids (transcription-reaction inputs), loaded once from the file at
# DS_GENE_STIDS_FILE (one R-HSA stable_id per line). Empty when unset — then no
# edge is ever classified transcriptional and scope=transcription is inert.
const _GENE_STIDS = Ref{Union{Nothing, Set{String}}}(nothing)
function gene_stid_set()::Set{String}
    if _GENE_STIDS[] === nothing
        s = Set{String}()
        path = get(ENV, "DS_GENE_STIDS_FILE", "")
        if !isempty(path) && isfile(path)
            for line in eachline(path)
                t = strip(line)
                isempty(t) || push!(s, String(t))
            end
        end
        _GENE_STIDS[] = s
    end
    return _GENE_STIDS[]::Set{String}
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
    gene_uuids::Set{String} = Set{String}();
    stats::Union{Nothing, Dict{String, Any}} = nothing,
    self_shared::Dict{Tuple{String, String}, Vector{String}} =
        Dict{Tuple{String, String}, Vector{String}}(),
)::Tuple{Vector{IndexedReaction}, Vector{Int}, Int}
    # Node indices that are gene entities (inputs to transcription/expression
    # reactions). Used to detect transcriptional autoregulation loops.
    gene_indices = Set{Int}(uuid_to_idx[u] for u in gene_uuids if haskey(uuid_to_idx, u))

    # First pass: build the bare IndexedReactions (without loop flags) and the
    # forward adjacency for short-loop detection.
    raw = NamedTuple[]
    n_nodes = length(uuid_to_idx)
    fwd_adj = [Set{Int}() for _ in 1:n_nodes]
    # Substrate/assembly chains only: X -> R -> P through input, assembly and
    # output edges. Catalyst edges are excluded (an enzyme is not consumed into
    # its product, so the product is not the enzyme's "own product"), and so
    # are composition edges (a hierarchy hop is not a producing reaction, so it
    # cannot make another hierarchy edge "redundant").
    chain_adj = [Set{Int}() for _ in 1:n_nodes]

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
        act_is_assembly = Bool[]
        act_group = Int[]
        act_orig = Int[]
        for (k, uuid) in enumerate(r.activator_uuids)
            haskey(uuid_to_idx, uuid) || continue
            push!(act_indices, uuid_to_idx[uuid])
            push!(act_is_and, k <= length(r.activator_is_and) ? r.activator_is_and[k] : r.is_and_gate)
            push!(act_is_catalyst, k <= length(r.activator_is_catalyst) ? r.activator_is_catalyst[k] : false)
            push!(act_is_assembly, k <= length(r.activator_is_assembly) ? r.activator_is_assembly[k] : false)
            push!(act_group, k <= length(r.activator_group) ? r.activator_group[k] : 0)
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

        # Forward adjacency: every input the target's value actually depends on
        # points at the target. This graph is what Tarjan runs on, so an edge
        # missing here is a cycle the SCC solver cannot see and a topological
        # order that can be wrong — a singleton component is evaluated exactly
        # once, with no iteration to correct a stale input.
        #
        # Substrates are included for that reason: compute_reaction_output_vec
        # reads x[substrate_indices] into the availability factor L, so they are
        # a real dependency. Nothing populates `substrate_uuids` today, which is
        # the only reason their absence was harmless.
        for (kk, a) in enumerate(act_indices)
            push!(fwd_adj[a], target_idx)
            is_cat = kk <= length(act_is_catalyst) && act_is_catalyst[kk]
            is_comp = kk <= length(act_group) && act_group[kk] > 0
            (is_cat || is_comp) || push!(chain_adj[a], target_idx)
        end
        for i in inh_indices
            push!(fwd_adj[i], target_idx)
        end
        for d in dep_indices
            push!(fwd_adj[d], target_idx)
        end
        for sb in sub_indices
            push!(fwd_adj[sb], target_idx)
        end

        # Per-inhibitor AND/OR flag, compacted in lockstep with inh_indices via
        # the same surviving-position list used for the inhibitor params. Falls
        # back to the reaction-level flag for any inhibitor whose per-edge flag
        # is missing (legacy Reaction constructors).
        inh_is_and = [k <= length(r.inhibitor_is_and) ? r.inhibitor_is_and[k] : r.is_and_gate
                      for k in inh_orig]
        # specs/022: which of this reaction's activators each inhibitor contains,
        # keyed by (inhibitor uuid, target uuid) and aligned with inh_indices.
        inh_shared = [Int[uuid_to_idx[a]
                          for a in get(self_shared, (r.inhibitor_uuids[k], r.target_uuid), String[])
                          if haskey(uuid_to_idx, a)]
                      for k in inh_orig]

        push!(raw, (
            target_idx=target_idx, baseline=get(baselines, r.target_uuid, 0.01),
            act_indices=act_indices, act_is_and=act_is_and,
            act_is_catalyst=act_is_catalyst, act_is_assembly=act_is_assembly,
            act_group=act_group,
            inh_indices=inh_indices, inh_is_and=inh_is_and, inh_shared=inh_shared,
            dep_indices=dep_indices,
            sub_indices=sub_indices, params=params,
        ))
    end

    # SCC component ids. Computed when either SCC feature is enabled (keeps the
    # default path byte-for-byte unchanged otherwise). comp_id[i] = SCC id of
    # node i, numbered in reverse-topological order. Used for (a) the SCC-
    # condensation solver's processing order [DS_SCC_SOLVE] and (b) marking
    # recycling catalyst back-edges to relax [DS_SCC_BREAK_CATALYST].
    scc_solve = _bool_env("DS_SCC_SOLVE", true)  # SCC-condensation solve is the default; set DS_SCC_SOLVE=0 for legacy flat iteration
    break_catalyst = _bool_env("DS_SCC_BREAK_CATALYST", false)
    # Loop elasticity reads `activator_in_loop`, i.e. SCC membership, so the
    # components are needed when it is on even under the legacy flat solve --
    # otherwise DS_LOOP_ELASTICITY<1 with DS_SCC_SOLVE=0 was a silent no-op.
    elastic = _float_env("DS_LOOP_ELASTICITY", 1.0) < 1.0
    # DS_SCC_BREAK_ROLES (specs/018): comma list from {catalyst, assembly,
    # depletion}. An edge of a listed role whose source and target share a
    # first-pass SCC is a recycling CLOSURE: it keeps its feed-forward meaning
    # (read at the component-entry value through `supply`) but does not feed a
    # value back around the cycle, and components are recomputed without those
    # edges so a component welded together by our own derived edges falls apart
    # into its reaction-level cycles. Measured motivation: removing assembly +
    # depletion edges takes TP53's component from 836 nodes to 34 and DSB's
    # from 1,127 to 126; MP-BioPath's hand-curated networks keep the reaction
    # backbone inside our cycles and essentially never carry these classes.
    break_roles = _break_roles_env()
    roles_on = !isempty(break_roles)
    comp_id, n_comp = (scc_solve || break_catalyst || elastic || roles_on) ?
        tarjan_scc_components(fwd_adj) : (Int[], 0)
    n_comp_pass1 = n_comp
    closure_act = [falses(length(rec.act_indices)) for rec in raw]
    closure_dep = [falses(length(rec.dep_indices)) for rec in raw]
    n_closure = Dict("catalyst" => 0, "assembly" => 0, "depletion" => 0)
    if roles_on && !isempty(comp_id)
        cyclic_before = count(c -> c > 0, let sz = zeros(Int, n_comp); (for c in comp_id; c >= 1 && (sz[c] += 1); end); [s > 1 ? 1 : 0 for s in sz] end)
        fwd_adj2 = [copy(s) for s in fwd_adj]
        @inbounds for (ri, rec) in enumerate(raw)
            t = rec.target_idx
            for k in eachindex(rec.act_indices)
                s = rec.act_indices[k]
                comp_id[s] == comp_id[t] || continue
                s == t && continue      # a self-loop is iterated by comp_has_self_loop, not a closure
                # `assembly` covers composition edges too (activator_is_assembly is
                # set for both); they are counted under closures_assembly.
                role = rec.act_is_catalyst[k] ? "catalyst" : (rec.act_is_assembly[k] ? "assembly" : "")
                if role in break_roles
                    closure_act[ri][k] = true; n_closure[role] += 1
                    # drop only if no other edge from s to t survives (Set semantics: one entry per pair)
                    delete!(fwd_adj2[s], t)
                end
            end
            for k in eachindex(rec.dep_indices)
                s = rec.dep_indices[k]
                if "depletion" in break_roles && comp_id[s] == comp_id[t] && s != t
                    closure_dep[ri][k] = true; n_closure["depletion"] += 1
                    delete!(fwd_adj2[s], t)
                end
            end
        end
        # an s->t pair may carry several edges; restore the adjacency for any pair that
        # still has a non-closure edge
        @inbounds for (ri, rec) in enumerate(raw)
            t = rec.target_idx
            for k in eachindex(rec.act_indices)
                closure_act[ri][k] || push!(fwd_adj2[rec.act_indices[k]], t)
            end
            for i in rec.inh_indices; push!(fwd_adj2[i], t); end
            for k in eachindex(rec.dep_indices)
                closure_dep[ri][k] || push!(fwd_adj2[rec.dep_indices[k]], t)
            end
            for sb in rec.sub_indices; push!(fwd_adj2[sb], t); end
        end
        comp_id, n_comp = tarjan_scc_components(fwd_adj2)
        if stats !== nothing
            sz = zeros(Int, n_comp); for c in comp_id; c >= 1 && (sz[c] += 1); end
            stats["scc_closures_catalyst"] = n_closure["catalyst"]
            stats["scc_closures_assembly"] = n_closure["assembly"]
            stats["scc_closures_depletion"] = n_closure["depletion"]
            stats["scc_cyclic_before"] = cyclic_before
            stats["scc_cyclic_after"] = count(>(1), sz)
            stats["scc_largest_after"] = isempty(sz) ? 0 : maximum(sz)
        end
    elseif stats !== nothing
        stats["scc_closures_catalyst"] = 0; stats["scc_closures_assembly"] = 0; stats["scc_closures_depletion"] = 0
        sz = zeros(Int, n_comp); for c in comp_id; c >= 1 && (sz[c] += 1); end
        stats["scc_cyclic_before"] = count(>(1), sz); stats["scc_cyclic_after"] = count(>(1), sz)
        stats["scc_largest_after"] = isempty(sz) ? 0 : maximum(sz)
    end

    # Second pass: for each (target, inhibitor) pair, BFS forward from the
    # target up to max_depth hops to see if the inhibitor is reachable. If so,
    # they sit in a small negative-feedback loop and the inhibitor edge is
    # eligible for DS_INHIBITOR_FLOOR dampening.
    max_depth = parse(Int, get(ENV, "DS_LOOP_DEPTH", "3"))
    indexed = IndexedReaction[]
    for (ri, rec) in enumerate(raw)
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
            # Also flag SCC membership: an inhibitor whose source sits in the
            # SAME strongly-connected component as the target is part of a
            # feedback cycle of ANY length — this catches the long
            # gene→mRNA→protein→transcription autoregulation loops (e.g. MDM2
            # repressing its own transcription) that the depth-bounded BFS
            # above misses. Such self-regulating negative feedback is where
            # continuous gene dosage + strong proportional repression is
            # biologically wrong, so these edges are the ones DS_INHIBITOR_FLOOR
            # should weaken (while feed-forward inhibition keeps full divide
            # strength). comp_id is empty when SCC detection is off → skipped.
            if !found && !isempty(comp_id) && comp_id[inh] == comp_id[rec.target_idx]
                found = true
            end
            push!(in_loop, found)
        end

        # Transcriptional repression: this reaction is a transcription/
        # expression step (it has a gene entity among its inputs), so any
        # negative regulator of it is repressing gene transcription. These get
        # weakened under scope=transcription because genes are ~on/off — strong
        # proportional repression of a continuous gene-dosage node is unphysical.
        # NOTE: the SCC/loop requirement is intentionally NOT applied: positional
        # decomposition severs the gene→protein→own-transcription autoregulation
        # loop into separate UUIDs, so these edges are feed-forward in the graph
        # (in_loop≈0) even though biologically they are self-regulation.
        rxn_has_gene = any(a -> a in gene_indices, rec.act_indices)
        in_transcription = fill(rxn_has_gene, length(rec.inh_indices))

        # Recycling-catalyst back-edges: a catalyst activator whose source is in
        # the same SCC as this reaction's target closes a recycling cycle. Mark
        # it for relaxation under DS_SCC_SOLVE. All-false when SCC solve is off.
        act_break = Bool[]
        for k in eachindex(rec.act_indices)
            brk = closure_act[ri][k]
            if break_catalyst && rec.act_is_catalyst[k]
                src = rec.act_indices[k]
                if comp_id[src] == comp_id[rec.target_idx]
                    brk = true
                end
            end
            push!(act_break, brk)
        end

        # Loop-closing activators: source and target share an SCC. This is the
        # set of edges whose gain product decides whether a cycle's baseline is
        # stable. Empty comp_id (SCC detection off) => all false.
        act_in_loop = Bool[]
        for k in eachindex(rec.act_indices)
            src = rec.act_indices[k]
            push!(act_in_loop,
                  !isempty(comp_id) && comp_id[src] == comp_id[rec.target_idx])
        end

        # Own-product depleters: P -| X where P is produced from X within two
        # activator hops (X -> reaction -> P, or X -> P directly). Structural,
        # order-free, independent of any DS_* setting; the propagator decides
        # what to do with it (DS_DEPLETION_OWN_PRODUCT).
        dep_own = Bool[]
        for d in rec.dep_indices
            own = d in chain_adj[rec.target_idx]
            if !own
                for mid in chain_adj[rec.target_idx]
                    if d in chain_adj[mid]
                        own = true
                        break
                    end
                end
            end
            push!(dep_own, own)
        end

        # Redundant composition inputs: source -> some reaction -> this target
        # already exists through activator edges, so the hierarchy edge repeats
        # a fold the producing route carries. Structural, order-free.
        comp_redundant = Bool[]
        for k in eachindex(rec.act_indices)
            s = rec.act_indices[k]
            red = false
            if k <= length(rec.act_group) && rec.act_group[k] > 0
                for mid in chain_adj[s]
                    if mid != rec.target_idx && rec.target_idx in chain_adj[mid]
                        red = true
                        break
                    end
                end
            end
            push!(comp_redundant, red)
        end

        push!(indexed, IndexedReaction(
            rec.target_idx, rec.baseline,
            rec.act_indices, rec.act_is_and, act_break, act_in_loop,
            rec.act_is_assembly, comp_redundant, rec.act_group,
            rec.inh_indices, rec.inh_is_and, in_loop, in_transcription, rec.inh_shared,
            rec.dep_indices, closure_dep[ri], dep_own,
            rec.sub_indices, rec.params,
        ))
    end
    if stats !== nothing
        stats["self_inhibitors"] = sum((count(!isempty, r.inhibitor_shared) for r in indexed); init = 0)
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
Solver knobs read on every reaction evaluation, resolved ONCE per solve from
the DS_* environment (with the validated winning-config defaults, commit
c9805b2) instead of re-reading + re-parsing ENV inside the per-reaction /
per-iteration hot loop. Concrete field types keep it type-stable.

The two DS_INHIBITOR_BETA fields intentionally carry the raw string with the
per-mode default ("1.0" for devspec, "" for spec) and are parsed lazily inside
the branch that uses them — this preserves the existing behavior exactly,
including that an explicitly-empty DS_INHIBITOR_BETA errors only if the devspec
branch actually runs.
"""
struct ReactionEvalConfig
    inhibition_mode::String
    and_mode::String
    or_mode::String
    assembly_limiting::Bool
    hill_sat_eps::Float64
    hill_sat_h_max::Float64
    hill_log_zmax::Float64
    inhibitor_k::Float64
    inhibitor_eps::Float64
    inhibitor_floor::Float64
    floor_scope::String
    devspec_beta_raw::String
    spec_beta_raw::String
    depletion_h_max::Float64
    depletion_h_min::Float64
    inhibitor_or::Bool
    or_redundancy::Float64
    or_combine::String
    loop_elasticity::Float64
    loop_elasticity_width::Float64
    loop_elasticity_hi::Float64
    composition_group::Bool
    composition_mode::String
    depletion_own_product::String
    self_inhibitor_weight::Float64   # < 0 = off (default). specs/022.
end

"""
Valid values for each `DS_*` mode variable.

Every mode dispatch below is an `if/elseif/else` whose `else` is a real model,
so an unrecognised value used to select a *different model* silently — and for
`DS_INHIBITION_MODE` that fallback is `"spec"`, whose default β is 0, i.e. no
inhibition at all. `DS_INHIBITION_MODE="Divide"` (or a trailing space, easily
produced by a shell or compose file) therefore turned inhibition off and still
reported a number. These sets make that a startup error instead.

Each set below is the exhaustive list of values the corresponding dispatch
actually implements, INCLUDING the value the `else` branch represents (named
explicitly so that selecting it stays deliberate):

- `DS_INHIBITION_MODE`: `spec` is the `else`.
- `DS_AND_MODE`: `geomean` is the `else`.
- `DS_OR_MODE`: `max` is the `else`.
- `DS_OR_COMBINE`: `max` is the `else`.
- `DS_DEPLETION_OWN_PRODUCT`: `full` is the `else` (own-product depleters
  both suppress and de-repress, as every depleter does today).
- `DS_COMPOSITION_MODE`: `assembly` is the `else` (composition edges join the
  assembly-limiting AND cluster); `limit` makes them a pure limiter.
- `DS_INHIBITOR_FLOOR_SCOPE`: `none` is the `else` (no floor applied).
"""
const DS_VALID_MODES = Dict(
    "DS_INHIBITION_MODE" => Set([
        "divide", "devspec", "spec", "krep", "hill_sat", "inversion",
    ]),
    "DS_AND_MODE" => Set([
        "hill_log", "hill_log_asym", "hill_sat", "multiplicative", "signed",
        "signed_gated",
        "min", "geomean",
    ]),
    "DS_OR_MODE"               => Set(["mean", "median", "capacity", "max"]),
    "DS_OR_COMBINE"            => Set(["max", "gate"]),
    "DS_COMPOSITION_MODE"      => Set(["assembly", "limit", "limit_novel"]),
    "DS_DEPLETION_OWN_PRODUCT" => Set(["full", "suppress_only"]),
    "DS_INHIBITOR_FLOOR_SCOPE" => Set(["loops", "all", "transcription", "none"]),
)

const DS_BREAK_ROLES = Set(["catalyst", "assembly", "depletion"])

"""
`DS_SCC_BREAK_ROLES`: comma-separated roles whose cycle-closing edges are read
at the component-entry value and excluded from component detection
(specs/018). Empty (default) = off. A misspelt role is a startup error.
`assembly` includes `composition` edges (both are assembly-class). Self-loops
are never closures. Under a role list, `DS_LOOP_ELASTICITY` and the inhibitor
loop floor see the recomputed (smaller) components.
"""
function _break_roles_env()::Set{String}
    raw = strip(get(ENV, "DS_SCC_BREAK_ROLES", ""))
    isempty(raw) && return Set{String}()
    roles = Set{String}()
    for tok in split(raw, ",")
        r = String(strip(tok))
        isempty(r) && continue
        r in DS_BREAK_ROLES || throw(ArgumentError(
            "DS_SCC_BREAK_ROLES=$raw: \"$r\" is not a role; expected a comma list from " *
            "\"catalyst\", \"assembly\", \"depletion\"."))
        push!(roles, r)
    end
    return roles
end

"""
Read a `DS_*` mode variable, rejecting anything not in its allowlist.

Names are matched exactly: silently lowercasing or trimming would just move the
guess one level up. A typo is a configuration error, not a model choice.
"""
function _mode_env(name::String, default::String)::String
    value = get(ENV, name, default)
    valid = DS_VALID_MODES[name]
    if !(value in valid)
        throw(ArgumentError(
            "$name=\"$value\" is not a recognised mode. Valid values: " *
            join(sort(collect(valid)), ", ") * ". " *
            "(Unrecognised values used to fall through to a different model silently.)"
        ))
    end
    return value
end

"""
Read a boolean `DS_*` flag.

Flags were previously tested two incompatible ways — some `== "1"`, others
`!= "0"` — so `DS_ASSEMBLY_LIMITING="true"` turned the feature OFF (a 74x
change in output) while `DS_SCC_SOLVE="false"` left it ON. Both spellings are
now accepted, and anything else is an error rather than a silent default.
"""
function _bool_env(name::String, default::Bool)::Bool
    raw = get(ENV, name, nothing)
    raw === nothing && return default
    value = lowercase(strip(raw))
    value in ("1", "true", "yes", "on") && return true
    value in ("0", "false", "no", "off") && return false
    throw(ArgumentError(
        "$name=\"$raw\" is not a boolean. Use one of 1/0, true/false, yes/no, on/off."
    ))
end

"""
Read a numeric `DS_*` variable, reporting the variable name on a bad value.

A bare `parse(Float64, ...)` raises `ArgumentError: cannot parse "" as Float64`,
which does not say which of the two dozen DS_* knobs was wrong.
"""
function _float_env(name::String, default::Float64)::Float64
    raw = get(ENV, name, nothing)
    raw === nothing && return default
    value = tryparse(Float64, strip(raw))
    if value === nothing || !isfinite(value)
        throw(ArgumentError("$name=\"$raw\" is not a finite number."))
    end
    return value
end

"""
Resolve a ReactionEvalConfig from the DS_* environment. Defaults are the
validated winning config (commit c9805b2). Call once per solve and thread the
result into compute_reaction_output_vec / forward_model_vec.

Invalid values raise `ArgumentError` here, at the single resolution point,
rather than silently selecting a different model deeper in the dispatch.
"""
function resolve_reaction_eval_config()::ReactionEvalConfig
    return ReactionEvalConfig(
        _mode_env("DS_INHIBITION_MODE", "divide"),
        # hill_sat. AND is multiplication of fold-changes capped at 100 --
        # 0.5*0.5 = 0.25, 0.1*0.1 = 0.01, 0*x = 0, 10*10 = 100, 100*100 = 100 --
        # and hill_sat now reproduces that exactly. See
        # specs/010-and-multiplication-fidelity/research.md.
        #
        # This costs 16 held-out cases (0.085pp, macro-F1 -0.0014) against
        # hill_log, measured on 23,908 conditioned cases. It buys correct
        # arithmetic below baseline, where hill_log read +35% at 0.1*0.1 and
        # +168% at 0.001, systematically UPWARD, and where `0 x anything` never
        # reached 0 and rose with the co-input.
        #
        # It also buys DEPTH-INVARIANT magnitudes. hill_log's tanh compression
        # compounds along a cascade: through single-input reactions -- which
        # have nothing to combine and should be the identity -- a 100x source
        # arrives as 74x at one hop and 19x at ten, so two readouts with
        # identical biology and different path lengths get different predicted
        # folds. That is what made hill_log's continuous outputs unusable as a
        # scale, which was the reason specs/002 adopted it.
        #
        # The prior record of "hill_sat is 90 cases worse" was measured against
        # a hill_sat that could not represent a knockout AND inverted its own
        # saturation on wide reactions. Corrected, the gap is 16.
        _mode_env("DS_AND_MODE", "hill_sat"),
        _mode_env("DS_OR_MODE", "mean"),
        # OFF. The limiting-reactant rule aggregated a complex's subunits with
        # min, and min(elevated, baseline) is EXACTLY baseline — so a complex
        # transmitted scarcity perfectly and blocked abundance completely.
        # Every complex on a path hard-clamped any increase, which is why
        # upregulation collapsed over distance while knockouts propagated.
        #
        # It was adopted for a +0.9pp curator win measured when perturbations
        # were pinned across a median of 17 nodes, some adjacent to the
        # readout — conditions where an increase barely had to cross a complex
        # and the clamp therefore cost almost nothing.
        _bool_env("DS_ASSEMBLY_LIMITING", true),
        # 1e-5, not 1e-3. The epsilon smooths the saturation corners, but at
        # 1e-3 it EXCEEDS the internal values where knockouts live (~0.0006)
        # and acts as a floor on the whole network: 0.25x0.25 read 0.09
        # instead of 0.0625, and the model predicted UP on 317 of 564
        # benchmark cases against 246 actually UP, median output 1.587 rather
        # than 1.000. At 1e-5 the curve reproduces pure multiplication across
        # the range while keeping the boundary smooth.
        # 1e-9, not 1e-5. This eps is the smooth-max floor against zero, so
        # it bounds how far below baseline a value can go. At 1e-5 against a
        # baseline of 0.01 it distorts the low end badly: a fold of 0.001 read
        # 20.7% high and a fold of 0.0001 read 452% high, while contributing
        # nothing above baseline. At 1e-9 every product is exact to 0.0% across
        # the whole range and 10x10 still caps at 100. Affects hill_sat only.
        _float_env("DS_HILL_SAT_EPS", 1e-9),
        _float_env("DS_HILL_SAT_H_MAX", 10.0),
        _float_env("DS_HILL_LOG_ZMAX", 10.0),
        _float_env("DS_INHIBITOR_K", 0.1),
        # 1e-12, not 1e-3. This is a divide-by-zero guard and nothing else,
        # so it must be orders of magnitude below the scale it guards.
        # Baseline is 0.01, so the old default was TEN PERCENT of that, and at
        # that size it was not guarding — it was doing three jobs nobody wrote
        # down: setting the de-repression ceiling ((bl+eps)/eps = 11x from a
        # single inhibitor knockout), compressing the whole interior of the
        # response (x = bl/2 read 1.833 instead of 2.0), and shifting maximum
        # suppression (0.0110 instead of 0.0100).
        #
        # The blow-up as x -> 0 was never unhandled: clamp(result, 0, 10) below
        # is the real de-repression ceiling. Shrinking eps just stops it
        # distorting everything above zero.
        #
        # DS_HILL_SAT_EPS had the identical defect at the identical value and
        # was fixed in specs/002-upregulation-propagation without anyone
        # checking the sibling. Any epsilon here must stay << 0.01.
        #
        # Measured, 89 pathways / 23,788 curator cases: 19,431 -> 19,459
        # correct, macro-F1 0.7814 -> 0.7835, and the lowest
        # false_positive_change of any arm tested. See
        # specs/006-bounded-derepression.
        _float_env("DS_INHIBITOR_EPS", 1e-12),
        _float_env("DS_INHIBITOR_FLOOR", 0.0),
        _mode_env("DS_INHIBITOR_FLOOR_SCOPE", "loops"),
        get(ENV, "DS_INHIBITOR_BETA", "1.0"),  # devspec default
        get(ENV, "DS_INHIBITOR_BETA", ""),     # spec default (per-edge when empty)
        _float_env("DS_DEPLETION_H_MAX", 10.0),
        # Lower bound on depletion suppression. Defaults to 1/h_max, i.e.
        # symmetric with the de-repression cap in log space: depletion may
        # suppress at most as hard as it may de-repress. Before this existed
        # the bound was ZERO, so an abundant complex could deplete its free
        # subunit without limit and multiple depletion edges compounded
        # multiplicatively. Configurable so the value can be swept on the
        # TUNING split without editing code -- the symmetric default is
        # principled only relative to h_max, which is itself a tuned constant.
        _float_env("DS_DEPLETION_H_MIN", 1.0 / _float_env("DS_DEPLETION_H_MAX", 10.0)),
        _bool_env("DS_INHIBITOR_OR", false),
        # Clamped: w is documented as [0,1], and w > 1 would invert a
        # knockout (mean 0.96 + (1-w)*min with w=2 gives fold ~1.9, i.e. a KO
        # driving the target UP). Easy to hit with a sweep typo.
        clamp(_float_env("DS_OR_REDUNDANCY", 1.0), 0.0, 1.0),
        _mode_env("DS_OR_COMBINE", "max"),
        # Elasticity applied to LOOP-CLOSING activator edges only (source and
        # target in the same SCC): the input fold is read as fold^eps.
        #
        # Why: AND is multiplication of folds, so every edge has elasticity 1
        # and a positive cycle has loop gain exactly 1 at baseline. Baseline is
        # then a knife-edge, not a stable state. Measured on a two-reaction
        # AND-loop: U = 1.01 drives A to 1.09 / 2.30 / 100.0 at 50 / 500 / 5000
        # sweeps; U = 0.90 collapses it to 0.0026 and reports converged. Any
        # leak into a loop rails it, which way depends on sweep order and
        # budget -- the mechanism behind label-dependence, "more iterations
        # makes it worse", and false-change concentrating in cyclic pathways.
        #
        # With eps < 1 on the closing edges the loop's baseline root is unique
        # and stable and zero becomes unstable (A^eps >> A near zero pulls it
        # back up), so collapse needs an actual in-loop knockout. Acyclic edges
        # keep exact multiplication, so cascade depth-invariance is untouched.
        # This is the design doc's per-input sensitivity (Section 2.1) applied
        # where it matters, and it is the first learnable loop parameter.
        #
        # 1.0 = off (byte-identical default). Must be in (0, 1]: 0 would erase
        # the input, > 1 would sharpen the knife-edge.
        _loop_elasticity_env(),
        # Width (in |log fold|) of the band around baseline inside which the
        # loop elasticity applies. 0 = constant elasticity on every in-loop
        # edge (byte-identical to the flag above alone). With w > 0 the
        # elasticity is SIGMOIDAL in the input's fold:
        #
        #     eps(f) = eps_lo + (1 - eps_lo) * tanh(|log f| / w)
        #
        # so at f ~ 1 the edge reads fold^eps_lo (drift decays), and once
        # |log f| >> w it reads fold^1 (a real signal passes unattenuated).
        # Motivation, measured: constant eps = 0.5 cut false change 72% in
        # the loop-heavy pathways but traded it 1:1 for MISSED change -- 98%
        # of the cases it broke were genuine changes flattened to NORMAL, with
        # a median fixed-point value of 1.56x. In a giant SCC nearly every
        # edge is "in-loop", so constant eps compressed every path through it.
        # The leaks to damp sit at ~1.0x; the signals to keep sit at >= 1.5x.
        # This is the design doc's Section 2.1 sensitivity alpha(x) in the fold
        # domain, and the threshold is the first loop parameter to LEARN.
        _loop_elasticity_width_env(),
        # Ceiling the sigmoid rises to. 1.0 (default, byte-identical) means an
        # above-band signal is read at fold^1 -- exact multiplication -- which
        # restores loop gain 1 and therefore the knife-edge for every signal
        # outside the band. Measured: at w=0.15 that brought label-dependence
        # BACK and worse (73 predictions moved under UUID relabelling, vs 14
        # for the original solver and 0 for constant eps), a 0.19pp swing from
        # node names alone. With eps_hi < 1 the loop's gain stays below 1 at
        # every amplitude, so the root is unique everywhere, while a real
        # signal still clears the classifier: 1.5x through ten in-loop edges at
        # eps_hi=0.9 reads ~1.19, above the 1.15 UP cutoff.
        _loop_elasticity_hi_env(),
        # Aggregate `composition` inputs that are node copies of ONE entity as
        # alternatives (max within the group) before the limiting-reactant min
        # across distinct components. Motivation: HDR's BCDX2 complex exists as
        # 33 node copies from 5 variant reactions, all feeding the same container
        # -- min over copies caps the container by whichever copy a perturbation
        # happened to reach. Off by default; measured as its own arm.
        _bool_env("DS_COMPOSITION_GROUP", false),
        # How a `composition` edge (component complex -> containing complex, LNG
        # along hasComponent) enters the container's activity. specs/016.
        #
        #   "assembly" (default, byte-identical): the edge is one more
        #       assembly-limiting AND input, and the AND cluster then meets the
        #       container's OR cluster (its producing reactions) through
        #       DS_OR_COMBINE. Under `max` a baseline component MASKS a DOWN
        #       coming through the producing reaction (IFN-gamma 34 and DAP12 30
        #       DOWN->NORM); under `gate` an over-expressed component MULTIPLIES
        #       into every container above it (DSB Repair: 108 false UPs from
        #       RAD52/ERCC1/ERCC4/MUS81 OE; -151, identical with loops relaxed).
        #       Both measured on the deduplicated catalog.
        #
        #   "limit": a container cannot exceed its scarcest component, and a
        #       hierarchy edge is not a producing route. The composition inputs
        #       leave the AND cluster; their folds are capped at 1 and the
        #       smallest multiplies the reaction's result. A knocked-out or
        #       reduced component pulls the container down; an over-expressed
        #       one changes nothing (its partners still limit it); a container
        #       with no producing reaction in the network (the severed
        #       Interferon alpha/beta branch) reads baseline x that limiter.
        #
        #   "limit_novel": `limit`, but a composition edge whose source already
        #       feeds a producing reaction of the target (S -> R -> T; 54% of the
        #       catalog's composition edges) is skipped. Measured under `limit`:
        #       RAD52 OE on DSB Repair zeroes 281 of 365 composition sources
        #       through 35 nested containers, unmoved by an inhibitor floor or
        #       by loop elasticity -- the component's fold enters each container
        #       twice (producing reaction and hierarchy edge) and squares at
        #       every level. Only the hierarchy hops with NO reaction (the
        #       severed Interferon alpha/beta branch) carry new information.
        _mode_env("DS_COMPOSITION_MODE", "assembly"),
        # What a depletion edge may do when its source is the target's OWN
        # product (X -> R -> P, P -| X: a complex depleting the free subunit it
        # is built from). specs/016, traced on AKT1-KO -> TP53.
        #
        #   "full" (default, byte-identical): the depleter's fold acts both ways.
        #       When the substrate's supply halves, its products halve with it,
        #       and the edge then reads "fewer consumers -> less depletion" and
        #       de-represses the substrate back up: nuclear p-MDM2 arrives at
        #       0.82 from a 0.50 supply through three such edges, and TP53's
        #       de-repression never happens. A mass balance does not do this --
        #       a product that is low BECAUSE the substrate is low restores
        #       nothing.
        #
        #   "suppress_only": an own-product depleter may suppress (the abundant
        #       complex draining its free subunit -- specs/011's EGFR:CBL case)
        #       but its factor is capped at 1, so it cannot de-repress. What is
        #       lost: de-repression when the complex fell because of its OTHER
        #       partner (TP53 KO -> less MDM2:TP53 -> more free MDM2). The
        #       propagator cannot tell the two apart at the node; which the
        #       curators expect more often is the A/B.
        _mode_env("DS_DEPLETION_OWN_PRODUCT", "full"),
        # specs/022. An inhibitor that CONTAINS one of its reaction's inputs
        # (a sequestering complex such as WIF1:WNT) moves with that input, and
        # `divide` then counts the input twice: one such inhibitor cancels it,
        # two invert it. With w set, the part of the inhibitor's fold explained
        # by the shared input is kept only at power w (split across the
        # reaction's self-contained inhibitors), so the input sets the
        # direction and the inhibitor damps it: the reaction reads x^(1-w).
        # The inhibitor's independent part keeps full strength. Unset = off,
        # byte-identical. w is a single structural weight, fixed a priori and
        # meant to be learned later.
        _self_inhibitor_weight_env(),
    )
end

function _self_inhibitor_weight_env()::Float64
    raw = strip(get(ENV, "DS_SELF_INHIBITOR_WEIGHT", ""))
    isempty(raw) && return -1.0
    w = _float_env("DS_SELF_INHIBITOR_WEIGHT", -1.0)
    0.0 <= w <= 1.0 || throw(ArgumentError(
        "DS_SELF_INHIBITOR_WEIGHT=$raw must be in [0, 1] (unset = off); a weight " *
        "above 1 would invert the input, below 0 would amplify the double count."))
    return w
end

function _loop_elasticity_hi_env()::Float64
    hi = _float_env("DS_LOOP_ELASTICITY_HI", 1.0)
    if !(0.0 < hi <= 1.0) || !isfinite(hi)
        throw(ArgumentError(
            "DS_LOOP_ELASTICITY_HI=$hi is out of range; must be in (0, 1]. " *
            "1.0 = the sigmoid rises to exact multiplication (knife-edge above the band)."
        ))
    end
    hi
end

function _loop_elasticity_width_env()::Float64
    w = _float_env("DS_LOOP_ELASTICITY_WIDTH", 0.0)
    if !(w >= 0.0) || !isfinite(w)
        throw(ArgumentError(
            "DS_LOOP_ELASTICITY_WIDTH=$w must be finite and >= 0 " *
            "(0 = constant elasticity; > 0 = sigmoidal band in |log fold|)."
        ))
    end
    w
end

function _loop_elasticity_env()::Float64
    eps = _float_env("DS_LOOP_ELASTICITY", 1.0)
    if !(0.0 < eps <= 1.0) || !isfinite(eps)
        throw(ArgumentError(
            "DS_LOOP_ELASTICITY=$eps is out of range; must be in (0, 1]. " *
            "1.0 is off; below 1 stabilises positive loops at baseline."
        ))
    end
    eps
end

"""
Compute reaction output from a flat activity vector. Same math as the dict
version (sensitivity transform → geomean activators × Hill-suppression
inhibitors × geomean substrates → Hill output), but generic-typed over the
element type so ForwardDiff Duals propagate.

`config` is resolved once per solve (see resolve_reaction_eval_config). It
defaults to a fresh resolve so back-compat callers (tests, dead modules) keep
working; the solve hot path passes it explicitly so ENV is read once, not once
per reaction per iteration.
"""
function compute_reaction_output_vec(x::AbstractVector{T}, rxn::IndexedReaction;
                                     supply::Union{Nothing,AbstractVector}=nothing,
                                     config::ReactionEvalConfig=resolve_reaction_eval_config()) where {T<:Real}
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
    inhibition_mode = config.inhibition_mode
    bl = T(rxn.target_baseline)

    and_vals = T[]
    and_wts = Float64[]
    or_vals = T[]

    # Assembly (member→complex) inputs, when DS_ASSEMBLY_LIMITING is on, are
    # aggregated with a limiting-reactant (min) rule and injected as a SINGLE
    # AND input below: a complex cannot exceed the abundance of its scarcest
    # subunit, so overexpressing one member of a many-subunit complex must not
    # drive the complex up. When off, assembly inputs fall through to the normal
    # AND/OR handling (byte-for-byte unchanged default).
    assembly_limiting = config.assembly_limiting
    assembly_min = typemax(T)
    have_assembly = false
    # DS_COMPOSITION_GROUP: per-group running max for composition inputs that
    # are copies of one entity; folded into assembly_min after the loop.
    # Allocated only when the feature is on: this function is the hot path
    # (once per reaction per sweep) and an unconditional Dict here measured a
    # 3.7x slowdown / 3.5x allocation on a 3,000-node cyclic solve.
    use_groups = config.composition_group
    group_max = use_groups ? Dict{Int, T}() : nothing
    # DS_COMPOSITION_MODE=limit: composition inputs leave the AND cluster and
    # become a fold factor <= 1 applied to the whole reaction (see resolver).
    comp_limit_mode = config.composition_mode == "limit" || config.composition_mode == "limit_novel"
    comp_skip_redundant = config.composition_mode == "limit_novel"
    comp_limit = one(T)
    comp_group_max = (comp_limit_mode && use_groups) ? Dict{Int, T}() : nothing

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
        # Loop elasticity (DS_LOOP_ELASTICITY): on a cycle-closing edge, read
        # the input's FOLD through fold^eps so the loop's gain at baseline is
        # below 1. Fold 1 maps to fold 1, 0 to 0; only the slope changes.
        if config.loop_elasticity < 1.0 &&
           k <= length(rxn.activator_in_loop) && rxn.activator_in_loop[k] &&
           src_val > zero(T)
            fold = src_val / bl
            eps_lo = T(config.loop_elasticity)
            eps = if config.loop_elasticity_width > 0.0
                # Sigmoidal: eps_lo at baseline, -> eps_hi once |log fold| >> width.
                eps_hi = T(config.loop_elasticity_hi)
                eps_lo + (eps_hi - eps_lo) * tanh(abs(log(fold)) / T(config.loop_elasticity_width))
            else
                eps_lo
            end
            src_val = bl * fold^eps
        end
        transformed = apply_sensitivity_transform(
            src_val;
            s=p.activator_sensitivity_s[k],
            n=p.activator_sensitivity_n[k],
            K_α=p.activator_sensitivity_K[k],
        )
        if comp_limit_mode && k <= length(rxn.activator_group) && rxn.activator_group[k] > 0
            # Limiter: the container cannot exceed this component. Fold capped
            # at 1, so an over-expressed component never lifts the container.
            if comp_skip_redundant && k <= length(rxn.activator_comp_redundant) &&
               rxn.activator_comp_redundant[k]
                # the producing reaction already carries this component's fold
            elseif comp_group_max !== nothing
                g = rxn.activator_group[k]
                comp_group_max[g] = max(get(comp_group_max, g, zero(T)), transformed)
            else
                comp_limit = min(comp_limit, min(one(T), transformed / bl))
            end
        elseif assembly_limiting && k <= length(rxn.activator_is_assembly) &&
           rxn.activator_is_assembly[k]
            g = (group_max !== nothing && k <= length(rxn.activator_group)) ? rxn.activator_group[k] : 0
            if g > 0
                # Copies of one entity are alternatives: the strongest copy carries it.
                group_max[g] = max(get(group_max, g, zero(T)), transformed)
            else
                # Limiting-reactant: track the scarcest subunit; injected once below.
                assembly_min = min(assembly_min, transformed)
            end
            have_assembly = true
        elseif rxn.activator_is_and[k]
            push!(and_vals, transformed)
            push!(and_wts, p.activator_weights[k])
        else
            push!(or_vals, transformed)
        end
    end

    # Inject the assembled-complex limiting value as a single AND input, so it
    # combines with any catalytic/regulatory inputs via the normal DS_AND_MODE.
    # Each entity group contributes its strongest copy to the limiting-reactant min.
    if group_max !== nothing
        for v in values(group_max)
            assembly_min = min(assembly_min, v)
        end
    end
    # Copies of one entity are alternatives (strongest copy), components co-limit.
    if comp_group_max !== nothing
        for v in values(comp_group_max)
            comp_limit = min(comp_limit, min(one(T), v / bl))
        end
    end
    if have_assembly
        push!(and_vals, assembly_min)
        push!(and_wts, 1.0)
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
            mode = config.and_mode
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
                # A required input at ZERO means the product is zero. The
                # hardcoded eps below is 1e-6 against a baseline of 0.01, so
                # without this short-circuit a knockout contributes a fold of
                # 1e-6/0.010001 = 1e-4 rather than 0, and `0 x 100` came out
                # at 0.0100 instead of 0 -- rising with the co-input, so a
                # knockout could be "rescued" by an abundant partner. Same
                # epsilon sizing bug already fixed in DS_HILL_SAT_EPS and
                # DS_INHIBITOR_EPS; here it is not even configurable.
                if any(v -> v <= zero(T), and_vals)
                    return_zero = true
                else
                    return_zero = false
                end
                # No epsilon in the ratio. Zero is handled above, and node
                # baselines are validated in (0, 1], so v/bl is safe for every
                # value that reaches here. The epsilon was a floor: at 1e-6
                # against baseline 0.01 it put a hard lower bound of 1e-4 on
                # any fold, so 0.001 x 1.0 read 0.001293 -- 29% high -- and
                # everything further below baseline was compressed upward.
                sat_eps = T(config.hill_sat_eps)
                max_internal = one(T)
                prod_fold = one(T)
                for v in and_vals
                    prod_fold *= v / bl
                end
                raw = bl * prod_fold
                # Smooth-min with max_internal (caps at 1.0), in the form that
                # does NOT catastrophically cancel.
                #
                # The textbook form (raw + max - sqrt(diff^2 + eps^2))/2 loses
                # `max` entirely once raw is large: for 100 AND inputs at fold
                # 2, raw is 1.3e28, sqrt(diff^2) == diff to machine precision,
                # and (raw + 1 - (raw - 1))/2 evaluates to 0 rather than 1. The
                # smooth-max below then returned eps/2, so a strongly ELEVATED
                # wide reaction read as ~0 instead of saturating at 100x --
                # the saturation inverted. It is why Class_I_MHC, which has
                # reactions carrying hundreds of nodes, failed to solve at all.
                #
                # Algebraically identical, evaluated stably:
                #   (raw + max - sqrt(d^2+e^2))/2  where d = raw - max
                #     = max + (d - sqrt(d^2+e^2))/2
                #     = max - e^2 / (2*(d + sqrt(d^2+e^2)))
                # The last form has no subtraction of nearby large numbers.
                diff_hi = raw - max_internal
                root_hi = sqrt(diff_hi * diff_hi + sat_eps * sat_eps)
                denom_hi = diff_hi + root_hi
                hi = denom_hi > zero(T) ?
                    max_internal - (sat_eps * sat_eps) / (T(2.0) * denom_hi) :
                    (raw + max_internal - root_hi) / T(2.0)
                # Smooth-max with 0 (floors at 0). hi is bounded by
                # max_internal here, so this form cannot cancel.
                lo = (hi + sqrt(hi * hi + sat_eps * sat_eps)) / T(2.0)
                return_zero ? zero(T) : clamp(lo, zero(T), one(T))
            elseif mode == "hill_log_asym"
                # AND that multiplies fold-changes FAITHFULLY BELOW baseline
                # and keeps hill_log's compression above it.
                #
                # hill_log applies `z_max*tanh(log_fold/z_max)` symmetrically,
                # so it compresses downward exactly as hard as upward. Measured
                # against the product it is supposed to compute:
                #
                #     0.5  x 2.0   ->  1.00005    (  0.0%)
                #     2.0  x 2.0   ->  3.9645     ( -0.9%)
                #     0.5  x 0.5   ->  0.2523     ( +0.9%)
                #     0.25 x 0.25  ->  0.0670     ( +7.2%)
                #     0.1  x 0.1   ->  0.0135     (+35.2%)
                #     0.001 x 1.0  ->  0.0027     (+167.6%)
                #     0    x 100   ->  0.0135     (never zero)
                #
                # Above baseline it is accurate to ~1%. Below baseline the
                # error grows without bound and is systematically UPWARD, so
                # every down-regulated value is lifted toward baseline. That is
                # a directional bias against detecting DOWN, not a rail
                # artifact.
                #
                # There is nothing to saturate against downward. The internal
                # domain is [0, 1] with baseline 0.01, so an upward fold is
                # genuinely capped at 100x -- that ceiling is real and tanh
                # models it. A downward fold of 0.01 is perfectly
                # representable, so compressing it is not modelling a floor,
                # it is error. Hence: sigmoid above baseline, exact product
                # below.
                #
                # Zero is handled exactly rather than through an epsilon. In
                # hill_log a knockout contributes log(1e-6/0.010001) = -9.21
                # instead of -Inf, which is why 0 x anything came out at
                # 0.0007 and not 0. The epsilon was hardcoded and, at 1e-6
                # against a baseline of 0.01, violated the sizing rule that
                # already had to be applied to DS_HILL_SAT_EPS and
                # DS_INHIBITOR_EPS.
                z_max = T(config.hill_log_zmax)
                if any(v -> v <= zero(T), and_vals)
                    # A required input is absent: the product IS zero.
                    zero(T)
                else
                    log_fold = zero(T)
                    for v in and_vals
                        log_fold += log(v / bl)
                    end
                    log_out = log_fold > zero(T) ?
                        z_max * tanh(log_fold / z_max) :   # cap the real 100x ceiling
                        log_fold                            # exact product downward
                    clamp(bl * exp(log_out), zero(T), one(T))
                end
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
                # z_max controlled by DS_HILL_LOG_ZMAX. Default 10.0 is the
                # validated winning config (commit c9805b2) — higher than the
                # log(100)≈4.605 used in the illustrative numbers above, i.e.
                # closer to pure multiplication (softer saturation).
                eps_T = T(1e-6)
                z_max = T(config.hill_log_zmax)
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
            or_mode = config.or_mode
            if or_mode == "mean"
                sum(or_vals) / length(or_vals)
            elseif or_mode == "median"
                # Plain median; ties are fine because clamp keeps us in [0,1].
                sorted = sort(collect(or_vals))
                n = length(sorted)
                isodd(n) ? sorted[(n + 1) ÷ 2] : (sorted[n ÷ 2] + sorted[n ÷ 2 + 1]) / 2
            elseif or_mode == "capacity"
                # Redundant-alternatives capacity, in FOLD space so it composes
                # through a cascade the way the hill_log AND path does.
                #
                # Motivation: a set-valued catalyst is "any one of these plays
                # this role", so losing one alternative should cost a share of
                # the route — not nothing (max: KO of one of N leaves the max at
                # baseline, masking it) and not everything (AND: any KO kills the
                # reaction outright).
                #
                #   fold_i   = (vᵢ + ε) / (baseline + ε)
                #   capacity = mean(fold_i)          # share of surviving routes
                #   dominant = min(fold_i)           # any lost route is decisive
                #   fold     = w·capacity + (1-w)·dominant
                #
                # w = DS_OR_REDUNDANCY ∈ [0,1]. w=1 credits full redundancy and
                # is equivalent to `mean`; w=0 makes any lost alternative
                # decisive (AND-like for a knockout). Intermediate values credit
                # partial redundancy, which is the biologically interesting
                # regime: Reactome sets enumerate every paralog, but in a given
                # cell line only a subset is expressed, so a KO of the
                # functionally relevant member behaves closer to a single point
                # of failure than 1/N of a route. Pure capacity gives fold 0.96
                # for a KO of 1-of-26 — above the DOWN cutoff, hence invisible.
                eps_T = T(1e-6)
                w = T(config.or_redundancy)
                inv_n = one(T) / T(length(or_vals))
                cap = zero(T)
                dom = typemax(T)
                for v in or_vals
                    f = (v + eps_T) / (bl + eps_T)
                    cap += inv_n * f
                    f < dom && (dom = f)
                end
                clamp(bl * (w * cap + (one(T) - w) * dom), zero(T), one(T))
            else
                maximum(or_vals)
            end
        end
        if and_result === nothing
            or_result
        elseif or_result === nothing
            and_result
        elseif config.or_combine == "gate"
            # DS_OR_COMBINE=gate: treat the OR cluster as a REQUIRED route whose
            # remaining capacity throttles the reaction, rather than an
            # alternative to the AND inputs.
            #
            # The default `max(and, or)` silently discards any OR-cluster loss:
            # every reaction has AND-clustered input edges sitting at baseline,
            # so max() falls back to those and the reaction proceeds unchanged.
            # Measured: with one AND input present, knocking out an OR-catalyst
            # alternative leaves the target at fold 1.0 for EVERY redundancy
            # weight, including a total OR-cluster kill. So no OR aggregator can
            # gate a reaction under `max` — the combination rule has to change.
            #
            # Gating multiplies the AND result by the OR cluster's fold-change,
            # which is what "this reaction needs one of these catalysts" means.
            eps_T = T(1e-6)
            clamp(and_result * ((or_result + eps_T) / (bl + eps_T)), zero(T), one(T))
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
    #               controlled by DS_INHIBITOR_EPS (default 1e-12, tiny vs
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
        K = T(config.inhibitor_k)
        result = one(T)
        @inbounds for k in 1:length(rxn.inhibitor_indices)
            x_inh = clamp(x[rxn.inhibitor_indices[k]], zero(T), one(T))
            m = T(p.inhibitor_ms[k])
            Km = K^m
            result *= Km / (Km + x_inh^m)
        end
        clamp(result, zero(T), one(T))
    elseif inhibition_mode == "divide"
        eps_T = T(config.inhibitor_eps)
        # DS_INHIBITOR_FLOOR caps how strongly an inhibitor edge can suppress
        # its target. DS_INHIBITOR_FLOOR_SCOPE selects which edges it applies to:
        #   "all"           — every inhibitor (global variant, for comparison)
        #   "loops"         — inhibitors in any negative-feedback loop (SCC or
        #                     short BFS)
        #   "transcription" — ONLY transcriptional autoregulation loops (the
        #                     reaction has a gene input and the inhibitor is in
        #                     its SCC). This is the biologically-typed rule:
        #                     weaken self-regulating gene transcription, leave
        #                     protein-level feedback at full divide strength.
        h_floor = T(config.inhibitor_floor)
        floor_scope = config.floor_scope
        # DS_INHIBITOR_OR: honor the per-edge inhibitor AND/OR flag. AND-clustered
        # repressors act cooperatively, so their suppression factors multiply
        # (the existing behaviour, applied to every inhibitor). OR-clustered ones
        # are ALTERNATIVE repressors — "any one of these represses" — so the
        # dominant (most suppressing, i.e. smallest H) one applies instead of the
        # product. Multiplying alternatives both over-represses when they are all
        # present and, worse, manufactures de-repression when only one is knocked
        # out even though the redundant partners still repress.
        honor_or = config.inhibitor_or
        result = one(T)
        or_h = one(T)
        have_or = false
        w_self = config.self_inhibitor_weight
        n_self = w_self >= 0 ? count(!isempty, rxn.inhibitor_shared) : 0
        @inbounds for k in 1:length(rxn.inhibitor_indices)
            x_inh = clamp(x[rxn.inhibitor_indices[k]], zero(T), one(T))
            if n_self > 0 && !isempty(rxn.inhibitor_shared[k])
                # specs/022: split the inhibitor's fold into the part its shared
                # inputs explain (kept at power w/n) and the independent rest
                # (kept whole). All folds are relative to the reaction baseline.
                f_s = one(T)
                for a in rxn.inhibitor_shared[k]
                    f_s *= clamp(x[a], zero(T), one(T)) / bl
                end
                f_i = x_inh / bl
                indep = f_s > T(1e-12) ? f_i / f_s : one(T)
                x_inh = clamp(bl * indep * f_s^(T(w_self) / n_self), zero(T), one(T))
            end
            h_k = (bl + eps_T) / (x_inh + eps_T)
            apply_floor = floor_scope == "all" ||
                          (floor_scope == "loops" && rxn.inhibitor_in_short_loop[k]) ||
                          (floor_scope == "transcription" && rxn.inhibitor_transcriptional[k])
            if apply_floor
                h_k = max(h_floor, h_k)
            end
            if honor_or && k <= length(rxn.inhibitor_is_and) && !rxn.inhibitor_is_and[k]
                or_h = have_or ? min(or_h, h_k) : h_k
                have_or = true
            else
                result *= h_k
            end
        end
        have_or && (result *= or_h)
        # THIS is the de-repression ceiling: the most that removing inhibitors
        # can raise a target at one reaction. It is also what makes the
        # division safe as x_inh -> 0, which is why DS_INHIBITOR_EPS can be
        # numerically tiny. Ten is an assumption, not a measurement — a
        # tighter value was tried and rejected (2x cost 39 cases on 23,788),
        # but it has never been justified from data either.
        clamp(result, zero(T), T(10.0))
    elseif inhibition_mode == "devspec"
        beta_env = config.devspec_beta_raw
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
        eps_T = T(config.inhibitor_eps)
        h_max = T(config.hill_sat_h_max)
        sat_eps = T(config.hill_sat_eps)
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
        beta_env = config.spec_beta_raw
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
        eps_dep = T(config.inhibitor_eps)
        h_max_dep = T(config.depletion_h_max)
        suppress_only = config.depletion_own_product == "suppress_only"
        @inbounds for k in 1:length(rxn.depletion_indices)
            i_dep = rxn.depletion_indices[k]
            raw_dep = (supply !== nothing && k <= length(rxn.depletion_break) && rxn.depletion_break[k]) ?
                      supply[i_dep] : x[i_dep]      # recycling closure: entry value, no feedback
            x_dep = clamp(raw_dep, zero(T), one(T))
            f_dep = (bl + eps_dep) / (x_dep + eps_dep)
            if suppress_only && k <= length(rxn.depletion_own_product) &&
               rxn.depletion_own_product[k]
                f_dep = min(f_dep, one(T))
            end
            H_dep *= f_dep
        end
        # Bound suppression by the SAME factor as de-repression. The old
        # lower bound was zero, so depletion could suppress a node without
        # limit while de-repression was capped at h_max -- capped above,
        # unbounded below, the same asymmetry the AND modes had.
        #
        # Traced from a real failure: knocking out EPS15 in Signaling_by_EGFR
        # de-represses two EGFR:CBL complexes to 100x baseline. Both carry a
        # negative depletion edge onto free GRB2-1, each contributing
        # bl/x = 1/100, compounding to 1e-4. GRB2-1:SOS1 followed it down, and
        # the AND against a genuine 70x EGFR signal produced 0.007 -- so the
        # readout read DOWN when the truth is UP. The model was asserting that
        # an abundant complex depletes its free subunit ten-thousand-fold.
        #
        # 1/h_max_dep makes the bound symmetric in log space: depletion may
        # suppress at most as hard as it may de-repress.
        H_dep = clamp(H_dep, T(config.depletion_h_min), h_max_dep)
    end

    # Same rationale as the dict-version compute_reaction_output above:
    # drop the spec's output Hill step in favor of clamped linear propagation
    # to keep baseline as a stable fixed point. Hill primitives remain
    # available for explicit per-reaction use after training.
    return clamp(A * H * H_dep * L * comp_limit, zero(T), one(T))
end

"""
Forward model F(x; θ) in vector form. Returns a new vector where each
target node carries the output of its driving reaction; nodes with no
incoming reaction keep their current value.
"""
function forward_model_vec(x::AbstractVector{T}, reactions::Vector{IndexedReaction};
                           supply::Union{Nothing,AbstractVector}=nothing,
                           config::ReactionEvalConfig=resolve_reaction_eval_config()) where {T<:Real}
    y = copy(x)  # preserves element type, including Dual
    for rxn in reactions
        y[rxn.target_idx] = compute_reaction_output_vec(x, rxn; supply=supply, config=config)
    end
    return y
end
