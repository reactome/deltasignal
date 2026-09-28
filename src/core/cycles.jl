"""
Interconversion cycles as conserved pools, solved at steady state (specs/039;
explained for readers in docs/MODEL.md §4).

A curated cycle A ⇄ B (RAS-GDP ⇄ RAS-GTP, protein ⇄ phospho-protein, …) is one
protein moved between forms. Multiplying fold-changes around A → F → B → R → A
gives the loop a gain of exactly 1 at baseline, which is the knife-edge that rails
or collapses loops. The model treats "how fast B is made" as "how much B there
is", and nothing supplies a scale. Conservation supplies it.

Each pool is a Markov chain:
- the states are the protein's forms;
- a transition i → j has rate k_ij · u_ij, where u_ij is the DRIVE of its reaction
  (the reaction evaluated with the existing semantics, with its source form held
  at baseline, divided by the reaction's baseline);
- the stationary distribution π (πQ = 0, Σπ = 1) is the split between forms.

Form i then reads baseline × s × π_i / π0_i, where s is the supply fold of the
protein. Each transition's reaction node reads its flux relative to baseline:
s × u_ij × π_i / π0_i.

**Baseline.** φ₀ (DS_CYCLE_PHI, default 0.1; Adam 2026-09-27) is the share of the
modified form at rest. It uses detailed balance, with each step away from the base
(resting) form holding φ₀ / (1 − φ₀) of its predecessor's share. The baseline
rates are k_ij = sqrt(π0_j / π0_i), which satisfy detailed balance, so π = π0 when
nothing is perturbed.

  DS_CYCLE_MODE=off      (default) previous behaviour: loops are iterated.
  DS_CYCLE_MODE=balance  pools listed in the bundle's pools.csv are solved as above.
"""

const DS_VALID_CYCLE_MODES = Set(["off", "balance"])

function cycle_mode()::String
    value = get(ENV, "DS_CYCLE_MODE", "off")
    value in DS_VALID_CYCLE_MODES || throw(ArgumentError(
        "DS_CYCLE_MODE must be one of $(join(sort(collect(DS_VALID_CYCLE_MODES)), ", ")); got $(repr(value))"))
    return value
end

function cycle_phi()::Float64
    raw = strip(get(ENV, "DS_CYCLE_PHI", "0.1"))
    phi = tryparse(Float64, raw)
    (phi !== nothing && 0.0 < phi < 1.0) || throw(ArgumentError(
        "DS_CYCLE_PHI=$(repr(raw)) must be a number strictly between 0 and 1."))
    return phi
end

"""A pool resolved to solver indices."""
struct IndexedPool
    forms::Vector{Int}                  # node indices of the forms
    pi0::Vector{Float64}                # baseline split, aligned with forms
    trans_from::Vector{Int}             # position in `forms`
    trans_to::Vector{Int}
    trans_rxn::Vector{Int}              # index into rxns_idx of the transition reaction
    trans_node::Vector{Int}             # node index of the transition reaction
    trans_k::Vector{Float64}            # baseline rate
    supply_rxns::Vector{Vector{Int}}    # per form: rxns_idx of its OUTSIDE producers' reaction nodes
end

"""
Resolve `pools` to solver indices. A pool is kept only if every form and every
transition reaction has a node, and every transition reaction has an
IndexedReaction; otherwise it is skipped and counted (reported as a fallback).
"""
function index_pools(pools, uuid_to_idx::Dict{String, Int}, rxns_idx, phi::Float64)
    by_target = Dict{Int, Int}(r.target_idx => ri for (ri, r) in enumerate(rxns_idx))
    out = IndexedPool[]
    skipped = 0
    ratio = phi / (1.0 - phi)
    for p in pools
        if !all(f -> haskey(uuid_to_idx, f), p.forms) || isempty(p.transitions)
            skipped += 1; continue
        end
        pos = Dict(f => i for (i, f) in enumerate(p.forms))
        # distance from the base form over the (undirected) transition graph
        adj = Dict(f => Set{String}() for f in p.forms)
        for (a, b, _) in p.transitions
            (haskey(adj, a) && haskey(adj, b)) || continue
            push!(adj[a], b); push!(adj[b], a)
        end
        dist = Dict(p.base => 0); frontier = [p.base]
        while !isempty(frontier)
            nxt = String[]
            for f in frontier, g in sort!(collect(adj[f]))
                haskey(dist, g) && continue
                dist[g] = dist[f] + 1; push!(nxt, g)
            end
            frontier = nxt
        end
        if length(dist) != length(p.forms)
            skipped += 1; continue           # disconnected pool: not solvable as one chain
        end
        w = [ratio ^ dist[f] for f in p.forms]
        pi0 = w ./ sum(w)
        tf = Int[]; tt = Int[]; tr = Int[]; tn = Int[]; tk = Float64[]
        ok = true
        mult = Dict{Tuple{Int, Int}, Int}()
        for (a, b, rx) in p.transitions
            (haskey(pos, a) && haskey(pos, b) && haskey(uuid_to_idx, rx)) || (ok = false; break)
            node = uuid_to_idx[rx]
            haskey(by_target, node) || (ok = false; break)
            i, j = pos[a], pos[b]
            push!(tf, i); push!(tt, j); push!(tr, by_target[node]); push!(tn, node)
            mult[(i, j)] = get(mult, (i, j), 0) + 1
        end
        ok || (skipped += 1; continue)
        for q in eachindex(tf)
            i, j = tf[q], tt[q]
            # parallel copies share the pair's baseline rate
            push!(tk, sqrt(pi0[j] / pi0[i]) / mult[(i, j)])
        end
        trans_nodes = Set(tn)
        supply = Vector{Vector{Int}}()
        for f in p.forms
            fi = uuid_to_idx[f]
            prod = Int[]
            if haskey(by_target, fi)
                for a in rxns_idx[by_target[fi]].activator_indices
                    a in trans_nodes && continue
                    push!(prod, a)
                end
            end
            push!(supply, prod)
        end
        push!(out, IndexedPool([uuid_to_idx[f] for f in p.forms], pi0, tf, tt, tr, tn, tk, supply))
    end
    return out, skipped
end

"""
Stationary distribution of the pool's chain at the current drives, and the
drives themselves. The drives are evaluated with each transition's source form
held at baseline (the rate is k · u · π_i; its π_i is the unknown).
"""
function pool_stationary(x::AbstractVector{Float64}, p::IndexedPool, rxns_idx, baseline_vec, config)
    n = length(p.forms)
    u = Vector{Float64}(undef, length(p.trans_from))
    for q in eachindex(p.trans_from)
        src = p.forms[p.trans_from[q]]
        keep = x[src]
        x[src] = baseline_vec[src]
        out = compute_reaction_output_vec(x, rxns_idx[p.trans_rxn[q]]; config = config)
        x[src] = keep
        bl = baseline_vec[p.trans_node[q]]
        u[q] = bl > 0 ? out / bl : 0.0
    end
    Q = zeros(n, n)
    for q in eachindex(p.trans_from)
        i, j = p.trans_from[q], p.trans_to[q]
        rate = p.trans_k[q] * u[q] + 1e-12       # irreducible even with a drive at 0
        Q[i, j] += rate; Q[i, i] -= rate
    end
    A = Matrix(transpose(Q)); A[n, :] .= 1.0
    b = zeros(n); b[n] = 1.0
    pi = A \ b
    pi = max.(pi, 0.0); pi ./= sum(pi)
    return pi, u
end

"""
Write the pool's forms and transition-reaction nodes from its steady state.
Pinned forms keep their value and fix the supply; otherwise the supply is the
mean fold of the forms' producers outside the pool (1 if there are none).
Returns the largest absolute change written.
"""
function update_pool!(x::Vector{Float64}, p::IndexedPool, rxns_idx, baseline_vec, obs_set, config)
    pi, u = pool_stationary(x, p, rxns_idx, baseline_vec, config)
    s = NaN
    for (k, f) in enumerate(p.forms)
        if f in obs_set
            rel = pi[k] / p.pi0[k]
            s = rel > 0 ? (x[f] / baseline_vec[f]) / rel : 0.0
            break
        end
    end
    if isnan(s)
        tot = 0.0; cnt = 0
        for prods in p.supply_rxns, a in prods
            bl = baseline_vec[a]
            bl > 0 || continue
            tot += x[a] / bl; cnt += 1
        end
        s = cnt == 0 ? 1.0 : tot / cnt
    end
    maxch = 0.0
    for (k, f) in enumerate(p.forms)
        f in obs_set && continue
        v = clamp(baseline_vec[f] * s * pi[k] / p.pi0[k], 0.0, 1.0)
        maxch = max(maxch, abs(v - x[f])); x[f] = v
    end
    for q in eachindex(p.trans_from)
        t = p.trans_node[q]
        t in obs_set && continue
        i = p.trans_from[q]
        v = clamp(baseline_vec[t] * s * u[q] * pi[i] / p.pi0[i], 0.0, 1.0)
        maxch = max(maxch, abs(v - x[t])); x[t] = v
    end
    return maxch
end
