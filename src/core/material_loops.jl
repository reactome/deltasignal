"""
Material-closed positive loops (specs/107): two kinds of loop the generator's
pool and carrier-loop tables miss, handed to the rules that already exist for
them.

**Metabolite pools** (`DS_METABOLITE_POOLS=on`, rule 1). A reaction node with
exactly one `input` and one `output` that are not held (cofactor, drug,
pathogen), both small molecules (`node_kind` `simple_entity`), converts one
into the other (a small molecule: every leaf an `R-ALL-` entity). Small molecules in one strongly connected set of that
conversion graph (LACT ⇄ PYR, R5P ⇄ RU5P, CHOL cytosol ⇄ extracellular) are one
conserved quantity redistributed by enzymes, and are solved as a pool
(specs/039): each conversion is a one-step path, enzyme-driven if it has a
catalyst. There is no active form, so every state is a resting state with an
equal baseline share. Nodes a generator pool already uses are left alone.

**Multi-step carriers** (`DS_CARRIER_MULTISTEP=strict|wide`, rules 2 and 3). A
free form F is carried when a walk from F along material edges (`input` into a
reaction, `output` out of one, the identity edges `variant_pool`,
`variant_split`, `specialization`), through nodes that carry F's core, reaches
a reaction R with an `output` edge into F within CYCLE_MAX_STEPS reactions: F
goes round and comes back unchanged (BAP31 in MHC class I loading,
microtubules, TRAF3:MAVS). Each (F, R) is added to the carrier-loop table,
closed when every activating in-edge of F comes from its releases, so the
carrier rules of specs/049/077/088 apply unchanged.
- `strict`: the core is F's macromolecule leaves (not `R-ALL-`); F must have
  one, and set nodes (`set_variant`, `variant_pool`) are not walked.
- `wide`: also a small-molecule F (its small-molecule leaves are the core: PE
  through LC3:PE), set nodes carry by their chosen members, and leaves compare
  by reference entity when the bundle ships `leaf_refs.csv` (GRIP2 at the
  membrane and in the cytosol). Without the table that part is inert and the
  rule reports it.
Free forms the generator already lists, pool-managed nodes and held nodes are
not candidates.
"""

const DS_VALID_METABOLITE_POOLS = Set(["off", "on"])
const DS_VALID_CARRIER_MULTISTEP = Set(["off", "strict", "wide"])
const MATERIAL_IDENTITY_EDGES = Set(["variant_pool", "variant_split", "specialization"])
const SET_NODE_KINDS = Set(["set_variant", "variant_pool"])

function metabolite_pools_mode()::String
    value = get(ENV, "DS_METABOLITE_POOLS", "off")
    value in DS_VALID_METABOLITE_POOLS || throw(ArgumentError(
        "DS_METABOLITE_POOLS must be one of off, on; got $(repr(value))"))
    return value
end

function carrier_multistep_mode()::String
    value = get(ENV, "DS_CARRIER_MULTISTEP", "off")
    value in DS_VALID_CARRIER_MULTISTEP || throw(ArgumentError(
        "DS_CARRIER_MULTISTEP must be one of $(join(sort(collect(DS_VALID_CARRIER_MULTISTEP)), ", ")); got $(repr(value))"))
    return value
end

held_uuids(network::ReactionNetwork) =
    union(cofactor_uuids(network), drug_uuids(network), pathogen_uuids(network))

"""Every node a generator pool uses: states, intermediates, step reaction copies, carriers."""
function pool_claimed(network::ReactionNetwork)::Set{String}
    out = Set{String}()
    network.pools === nothing && return out
    for p in network.pools
        union!(out, p.forms, p.intermediates)
        for q in p.paths, (_, _, rxs) in q.steps; union!(out, rxs); end
        for (c, _) in p.carriers; push!(out, c); end
    end
    return out
end

# strongly connected sets (iterative Tarjan) of a graph given as adjacency lists
function _scc_sets(nodes::Vector{String}, adj::Dict{String, Vector{String}})
    index = Dict{String, Int}(); low = Dict{String, Int}(); onstack = Set{String}()
    stack = String[]; out = Vector{Vector{String}}(); i = 0
    for root in nodes
        haskey(index, root) && continue
        work = [(root, 1)]
        while !isempty(work)
            v, k = pop!(work)
            if k == 1
                i += 1; index[v] = i; low[v] = i; push!(stack, v); push!(onstack, v)
            end
            nbrs = get(adj, v, String[])
            recurse = false
            for j in k:length(nbrs)
                w = nbrs[j]
                if !haskey(index, w)
                    push!(work, (v, j + 1)); push!(work, (w, 1)); recurse = true; break
                elseif w in onstack
                    low[v] = min(low[v], index[w])
                end
            end
            recurse && continue
            if low[v] == index[v]
                comp = String[]
                while true
                    w = pop!(stack); delete!(onstack, w); push!(comp, w)
                    w == v && break
                end
                push!(out, sort!(comp))
            end
            isempty(work) || (low[work[end][1]] = min(low[work[end][1]], low[v]))
        end
    end
    return out
end

"""
Rule 1: the metabolite pools of `network` (see the module docstring), in a
stable order. Empty without a node table.
"""
function metabolite_pools(network::ReactionNetwork)::Vector{CyclePool}
    leaves = network.node_leaves
    leaves === nothing && return CyclePool[]
    held = held_uuids(network)
    claimed = pool_claimed(network)
    # a protein node is a `simple_entity` too: small means every leaf is a small molecule
    small(u) = haskey(leaves, u) && leaves[u][1] == "simple_entity" && !isempty(leaves[u][2]) &&
               all(l -> startswith(l, "R-ALL-"), leaves[u][2]) && !(u in held) && !(u in claimed)
    ins = Dict{String, Vector{String}}(); outs = Dict{String, Vector{String}}()
    catalysed = Set{String}()
    for e in network.edges
        if e.edge_type == "input" && !(e.parent_uuid in held)
            push!(get!(ins, e.child_uuid, String[]), e.parent_uuid)
        elseif e.edge_type == "output" && !(e.child_uuid in held)
            push!(get!(outs, e.parent_uuid, String[]), e.child_uuid)
        elseif e.edge_type == "catalyst" && e.is_positive
            push!(catalysed, e.child_uuid)
        end
    end
    # conversions a -> b, grouped by the reaction's stable id: the copies of one
    # curated reaction are one step's parallel copies, as the generator's pools
    conv = Dict{Tuple{String, String, String}, Vector{String}}()
    for (rx, i) in ins
        (rx in claimed || !haskey(outs, rx)) && continue
        a = unique(i); b = unique(outs[rx])
        (length(a) == 1 && length(b) == 1 && a[1] != b[1] && small(a[1]) && small(b[1])) || continue
        nd = get(network.nodes, rx, nothing)
        key = nd === nothing || nd.reactome_id === nothing ? rx : String(first(split(nd.reactome_id, "::")))
        push!(get!(conv, (a[1], b[1], key), String[]), rx)
    end
    adj = Dict{String, Vector{String}}()
    for (a, b, _) in keys(conv); push!(get!(adj, a, String[]), b); end
    for v in values(adj); sort!(unique!(v)); end
    pools = CyclePool[]
    for comp in _scc_sets(sort!(collect(keys(adj))), adj)
        length(comp) >= 2 || continue
        cs = Set(comp)
        paths = PoolPath[]
        for k in sort!([k for k in keys(conv) if k[1] in cs && k[2] in cs])
            rxs = sort(conv[k])
            push!(paths, PoolPath(k[1], k[2], [(k[1], k[2], rxs)], any(in(catalysed), rxs)))
        end
        push!(pools, CyclePool("metabolite:" * comp[1], comp, comp[1], paths, String[],
                               Tuple{String, String}[], comp))
    end
    return pools
end

"""
Rules 2 and 3: the (free form, release, closed) rows multi-step carriers add
to the carrier-loop table under `mode` (see the module docstring), and whether
the reference-entity comparison was available.
"""
function multistep_carriers(network::ReactionNetwork, mode::String)
    leaves = network.node_leaves
    rows = CarrierLoopEdge[]
    (mode == "off" || leaves === nothing) && return rows, false
    wide = mode == "wide"
    refs = wide ? network.leaf_refs : nothing
    ref(l) = refs === nothing ? l : get(refs, l, l)
    macro_(l) = !startswith(l, "R-ALL-")
    held = held_uuids(network)
    claimed = pool_claimed(network)
    listed = network.carrier_loops === nothing ? Set{String}() : Set(f for (f, _, _) in network.carrier_loops)
    kind(u) = haskey(leaves, u) ? leaves[u][1] : ""
    refset(u) = Set(ref(l) for l in (haskey(leaves, u) ? leaves[u][2] : String[]))
    out_adj = Dict{String, Vector{Tuple{String, String}}}()   # node -> (child, edge type)
    producers = Dict{String, Vector{String}}()                # activating in-edges of entity nodes
    for e in network.edges
        push!(get!(out_adj, e.parent_uuid, Tuple{String, String}[]), (e.child_uuid, e.edge_type))
        e.is_positive && !(e.edge_type in ("depletion", "catalyst", "regulator")) && kind(e.child_uuid) != "reaction" &&
            push!(get!(producers, e.child_uuid, String[]), e.parent_uuid)
    end
    for v in values(out_adj); sort!(v); end
    for f in sort!(collect(keys(producers)))
        k = kind(f)
        (k == "" || k == "reaction" || f in held || f in claimed || f in listed) && continue
        k in SET_NODE_KINDS && !wide && continue
        fl = haskey(leaves, f) ? leaves[f][2] : String[]
        core = Set(ref(l) for l in fl if macro_(l))
        if isempty(core)
            (wide && k == "simple_entity") || continue
            core = Set(ref(l) for l in fl)
            isempty(core) && continue
        end
        carries(u) = !(u in held) && (wide || !(kind(u) in SET_NODE_KINDS)) && issubset(core, refset(u))
        # breadth-first by reactions crossed
        releases = Set{String}()
        depth = Dict{String, Int}(f => 0)
        frontier = [f]
        while !isempty(frontier)
            nxt = String[]
            for u in frontier, (w, t) in get(out_adj, u, Tuple{String, String}[])
                if kind(u) == "reaction"
                    t == "output" || continue
                    if w == f
                        push!(releases, u); continue
                    end
                    (haskey(depth, w) || !carries(w)) && continue
                    depth[w] = depth[u]; push!(nxt, w)
                elseif t == "input" && kind(w) == "reaction"
                    (haskey(depth, w) || w in claimed || depth[u] >= CYCLE_MAX_STEPS) && continue
                    depth[w] = depth[u] + 1; push!(nxt, w)
                elseif t in MATERIAL_IDENTITY_EDGES
                    w == f && continue
                    (haskey(depth, w) || !carries(w)) && continue
                    depth[w] = depth[u]; push!(nxt, w)
                end
            end
            frontier = nxt
        end
        isempty(releases) && continue
        closed = all(in(releases), producers[f])
        for r in sort!(collect(releases))
            push!(rows, (f, r, closed))
        end
    end
    return rows, refs !== nothing
end

"""
`(network, report)`: `network` with rule 1's pools added to its pool table and
rules 2–3's rows added to its carrier-loop table, and what was added. Applied to
the network as given (before recycle edges are stripped), since rebuilding a
network drops `leaf_refs`.
"""
function apply_material_loops(network::ReactionNetwork;
                              mpools::String = metabolite_pools_mode(),
                              multistep::String = carrier_multistep_mode())
    report = Dict{String, Any}("metabolite_pools_rule" => mpools, "metabolite_pools" => 0,
                               "carrier_multistep_rule" => multistep, "carrier_multistep_forms" => 0,
                               "carrier_multistep_rows" => 0)
    pools = network.pools; loops = network.carrier_loops
    if mpools == "on"
        if network.node_leaves === nothing
            report["metabolite_pools_rule"] = "on: no node table"
        else
            mp = metabolite_pools(network)
            report["metabolite_pools"] = length(mp)
            isempty(mp) || (pools = vcat(something(pools, CyclePool[]), mp))
        end
    end
    if multistep != "off"
        if network.node_leaves === nothing
            report["carrier_multistep_rule"] = "$multistep: no node table"
        else
            # rule 1's pools are pool-managed for rule 2
            net = pools === network.pools ? network :
                  ReactionNetwork(network.nodes, network.edges, network.set_mappings, network.cofactor_stids,
                                  network.containment, network.drug_stids, pools, network.pathogen_stids,
                                  network.carrier_loops, network.node_leaves, network.rna_stids, network.leaf_refs)
            rows, with_refs = multistep_carriers(net, multistep)
            multistep == "wide" && !with_refs && (report["carrier_multistep_rule"] = "wide: no reference table")
            report["carrier_multistep_rows"] = length(rows)
            report["carrier_multistep_forms"] = length(unique(first.(rows)))
            isempty(rows) || (loops = vcat(something(loops, CarrierLoopEdge[]), rows))
        end
    end
    (pools === network.pools && loops === network.carrier_loops) && return network, report
    out = ReactionNetwork(network.nodes, network.edges, network.set_mappings, network.cofactor_stids,
                          network.containment, network.drug_stids, pools, network.pathogen_stids,
                          loops, network.node_leaves, network.rna_stids, network.leaf_refs)
    return out, report
end
