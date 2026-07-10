#!/usr/bin/env julia
# Loop / cycle diagnostics for the DeltaSignal penalty solver.
#
# WHY THIS EXISTS
# ---------------
# The penalty solver (src/solvers/steady_state.jl) is a fixed-point iteration
# x ← F(x) (optionally damped x ← (1-λ)x + λF(x) via DS_DAMPING). For ACYCLIC
# networks plain Jacobi (λ=0) converges in `depth` iterations and is exact.
# For CYCLIC sub-networks (biological feedback loops) it can oscillate and just
# reports the last iterate.
#
# This tool measures the two numbers that decide how to handle loops:
#   1. SCC structure   — how much of each network is actually inside a cycle
#                        (Tarjan on the activator+inhibitor+depletion graph).
#   2. Condensation depth — longest path through the DAG of SCCs. This is the
#                        number of Jacobi sweeps an UNDAMPED solve needs to
#                        converge. Under global damping λ it needs ~depth/λ
#                        sweeps, so if depth/λ > max_iters the acyclic majority
#                        of the network never converges — which is the hidden
#                        cost of global DS_DAMPING.
#
# HOW TO RUN (inside the julia-api container; ./bench is mounted at /app/bench)
# ---------------------------------------------------------------------------
#   docker exec deltasignal-julia-api-1 \
#       julia /app/bench/analysis/loop_diagnostics.jl
#
# Optional: restrict to specific catalog dirs by passing them as args:
#   ... loop_diagnostics.jl Signaling_by_WNT_R-HSA-195721
#
# Reads DS_* env from the container, so rerun under different DS_DAMPING etc.
# to compare. Catalog is /app/pathway_catalog (the PATHWAY_CATALOG mount).

using Pkg
Pkg.activate("/app")
using DeltaSignal
const DS = DeltaSignal

const CATALOG = get(ENV, "DS_DIAG_CATALOG", "/app/pathway_catalog")

# Default = the 9-pathway experimental benchmark set.
const DEFAULT_PATHWAYS = [
    "Transcriptional_Regulation_by_TP53_R-HSA-3700989",
    "Mitotic_G1-G1_S_phases_R-HSA-453279",
    "S_Phase_R-HSA-69242",
    "Signaling_by_ERBB2_R-HSA-1227986",
    "Cell_Cycle_Checkpoints_R-HSA-69620",
    "Signaling_by_WNT_R-HSA-195721",
    "PIP3_activates_AKT_signaling_R-HSA-1257604",
    "Mitotic_Prophase_R-HSA-68875",
    "HDR_through_Homologous_Recombination_HRR_or_Single_Strand_Annealing_SSA__R-HSA-5693567",
]

"""Directed dependency graph: edge u → t whenever node u is an input
(activator, inhibitor, or depletion) to a reaction producing t."""
function build_dep_graph(reactions, uuid_to_idx)
    n = length(uuid_to_idx)
    adj = [Int[] for _ in 1:n]
    for r in reactions
        haskey(uuid_to_idx, r.target_uuid) || continue
        t = uuid_to_idx[r.target_uuid]
        for u in Iterators.flatten((r.activator_uuids, r.inhibitor_uuids, r.depletion_uuids))
            haskey(uuid_to_idx, u) || continue
            push!(adj[uuid_to_idx[u]], t)
        end
    end
    return adj
end

"""Iterative Tarjan SCC. Returns Vector of components (each a Vector{Int})."""
function tarjan_scc(adj::Vector{Vector{Int}})
    n = length(adj)
    index = fill(0, n); low = fill(0, n); onstack = falses(n)
    comp_id = fill(0, n)
    stack = Int[]; comps = Vector{Int}[]; counter = 0; ncomp = 0
    for s in 1:n
        index[s] != 0 && continue
        work = [(s, 1)]
        while !isempty(work)
            v, pi = work[end]
            if pi == 1
                counter += 1; index[v] = counter; low[v] = counter
                push!(stack, v); onstack[v] = true
            end
            advanced = false
            i = pi
            while i <= length(adj[v])
                w = adj[v][i]
                if index[w] == 0
                    work[end] = (v, i + 1)
                    push!(work, (w, 1))
                    advanced = true
                    break
                elseif onstack[w]
                    low[v] = min(low[v], index[w])
                end
                i += 1
            end
            advanced && continue
            # all neighbours processed
            pop!(work)
            if !isempty(work)
                p = work[end][1]
                low[p] = min(low[p], low[v])
            end
            if low[v] == index[v]
                ncomp += 1
                comp = Int[]
                while true
                    w = pop!(stack); onstack[w] = false; comp_id[w] = ncomp; push!(comp, w)
                    w == v && break
                end
                push!(comps, comp)
            end
        end
    end
    return comps, comp_id
end

"""Longest path (in #edges) through the SCC condensation DAG = number of
Jacobi sweeps an undamped solve needs. Tarjan emits comps in reverse
topological order, so we can DP over that order."""
function condensation_depth(adj, comps, comp_id)
    ncomp = length(comps)
    # condensed adjacency (between distinct comps)
    cadj = [Set{Int}() for _ in 1:ncomp]
    for u in 1:length(adj)
        cu = comp_id[u]
        for v in adj[u]
            cv = comp_id[v]
            cu != cv && push!(cadj[cu], cv)
        end
    end
    # comps[] is in reverse-topo order (sinks first). Process in topo order:
    order = reverse(1:ncomp)
    dist = fill(0, ncomp)
    for c in order
        for d in cadj[c]
            dist[d] = max(dist[d], dist[c] + 1)
        end
    end
    return isempty(dist) ? 0 : maximum(dist)
end

"""Run the real penalty solver with `pin_idx` pinned to `pin_val` (a stand-in
perturbation) and report convergence. Also returns the per-node final state so
callers can diff damped vs undamped. Honours DS_DAMPING from env."""
function probe_convergence(net, reactions, all_nodes, uuid_to_idx, pin_idx, pin_val)
    pin_uuid = all_nodes[pin_idx]
    obs = Dict{String,Tuple{Float64,Float64}}(pin_uuid => (pin_val, 1.0))
    params = DS.default_steady_state_params()
    res = DS.solve_steady_state(net, obs, params)
    return res
end

"""Highest out-degree node = a hub whose perturbation propagates widely."""
function hub_node(adj)
    degs = [length(a) for a in adj]
    return argmax(degs)
end

pathways = length(ARGS) > 0 ? ARGS : DEFAULT_PATHWAYS

println("# DS_DAMPING=$(get(ENV,"DS_DAMPING","0.0"))  max_iters=500")
println(rpad("pathway", 52), "\t",
        join(["nodes","rxns","loop_sccs","largest","%in_loops","dag_depth","sweeps_needed_λ0.3"], "\t"))

for pw in pathways
    dir = joinpath(CATALOG, pw)
    if !isdir(dir)
        println(rpad(pw, 52), "\tMISSING"); continue
    end
    net = DS.parse_complete_network(
        joinpath(dir, "logic_network.csv"),
        joinpath(dir, "stid_to_uuid_mapping.csv"),
        nothing,
    )
    reactions = DS.convert_to_reaction_network(net)
    all_nodes = collect(keys(net.nodes))
    uuid_to_idx = Dict(u => i for (i, u) in enumerate(all_nodes))
    adj = build_dep_graph(reactions, uuid_to_idx)
    comps, comp_id = tarjan_scc(adj)
    nontrivial = filter(c -> length(c) > 1, comps)
    largest = isempty(nontrivial) ? 0 : maximum(length, nontrivial)
    in_loops = isempty(nontrivial) ? 0 : sum(length, nontrivial)
    pct = round(100 * in_loops / length(all_nodes), digits=1)
    depth = condensation_depth(adj, comps, comp_id)
    sweeps_damped = round(Int, depth / 0.3)
    println(rpad(pw, 52), "\t",
            join([length(all_nodes), length(reactions), length(nontrivial),
                  largest, pct, depth, sweeps_damped], "\t"))

    # --- Composition of the largest SCC -----------------------------------
    # If a handful of stable_ids account for most member nodes, the SCC is a
    # tangle of REUSED entities (transcription machinery, recycled cofactors)
    # = artifactual recycling that should have been broken. If member nodes
    # map to many distinct stable_ids, it's a genuine multi-entity feedback
    # structure. stid_to_uuid maps BOTH entity UUIDs and reaction (VR) UUIDs.
    if get(ENV, "DS_DIAG_SCC_COMPOSITION", "0") == "1" && !isempty(nontrivial)
        # load uuid -> stable_id
        stid = Dict{String,String}()
        for line in eachline(joinpath(dir, "stid_to_uuid_mapping.csv"))
            startswith(line, "uuid,") && continue
            parts = split(line, ',')
            length(parts) >= 2 && (stid[parts[1]] = parts[2])
        end
        idx_to_uuid = Dict(i => u for (u, i) in uuid_to_idx)
        big = nontrivial[argmax(length.(nontrivial))]
        # count nodes per stable_id within the SCC
        per_stid = Dict{String,Int}()
        for nidx in big
            u = idx_to_uuid[nidx]
            s = get(stid, u, "(no-stid)")
            per_stid[s] = get(per_stid, s, 0) + 1
        end
        ndistinct = length(per_stid)
        println("    largest SCC: $(length(big)) nodes, $ndistinct distinct stable_ids ",
                "(ratio $(round(length(big)/ndistinct, digits=1)) nodes/stid)")
        top = sort(collect(per_stid), by = x -> -x[2])[1:min(12, ndistinct)]
        for (s, c) in top
            println("      $(rpad(s,16)) $c nodes")
        end
    end
end
