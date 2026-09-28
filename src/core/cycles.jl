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

**Multi-step cycles (amendment 2).** Reactome usually curates a modification in
steps (S + E → S:E → S*:E → S* + E). The chain's states are the least-bound form
of each modification signature; a transition is a PATH state → intermediates →
state, with drive = the product of its steps' drives, and the intermediates
(enzyme complexes) hold no baseline share and read the flux through them. This
reduces exactly to the two-form rule. Baseline rates come from the stationary
flow of a random walk on the states (exits weighted 1 if enzyme-driven, else
CYCLE_INTRINSIC_WEIGHT beside an enzyme-driven exit), so π0 is stationary on any
graph, reversible or not.

**Carriers.** The enzyme's free form is fed back by the release steps it
catalyses; that closes a gain-1 loop of its own. Under DS_CYCLE_CARRIERS=on
(default) the free form reads its OTHER producers (else baseline).

  DS_CYCLE_MODE=off      (default) previous behaviour: loops are iterated.
  DS_CYCLE_MODE=balance  pools listed in the bundle's pools.csv are solved as above.
  DS_CYCLE_CARRIERS=on|off  (balance only) the carrier rule above.
"""

const DS_VALID_CYCLE_MODES = Set(["off", "balance"])

function cycle_mode()::String
    value = get(ENV, "DS_CYCLE_MODE", "off")
    value in DS_VALID_CYCLE_MODES || throw(ArgumentError(
        "DS_CYCLE_MODE must be one of $(join(sort(collect(DS_VALID_CYCLE_MODES)), ", ")); got $(repr(value))"))
    return value
end

function cycle_carriers()::Bool
    value = get(ENV, "DS_CYCLE_CARRIERS", "on")
    value in ("on", "off") || throw(ArgumentError("DS_CYCLE_CARRIERS must be on or off; got $(repr(value))"))
    return value == "on"
end

function cycle_phi()::Float64
    raw = strip(get(ENV, "DS_CYCLE_PHI", "0.1"))
    phi = tryparse(Float64, raw)
    (phi !== nothing && 0.0 < phi < 1.0) || throw(ArgumentError(
        "DS_CYCLE_PHI=$(repr(raw)) must be a number strictly between 0 and 1."))
    return phi
end

"""A pool resolved to solver indices. Positions refer to `forms` (the states)."""
struct IndexedPool
    forms::Vector{Int}                  # node indices of the states
    pi0::Vector{Float64}                # baseline split over states
    path_from::Vector{Int}              # position in `forms`
    path_to::Vector{Int}
    path_k::Vector{Float64}             # baseline rate
    path_j0::Vector{Float64}            # baseline flow
    path_steps::Vector{Vector{Tuple{Int, Int}}}   # (rxns_idx of the step's reaction, its source node)
    flux_nodes::Vector{Int}             # step reaction nodes and intermediates
    flux_paths::Vector{Vector{Int}}     # the paths through each
    supply_rxns::Vector{Vector{Int}}    # per state: its producers outside the pool
    form_rxn::Vector{Int}               # per state: rxns_idx of the state's own update (0 if none)
    carriers::Vector{Int}               # enzyme free-form nodes
    carrier_rxn::Vector{Int}            # rxns_idx of each carrier's update
    carrier_release::Vector{Set{Int}}   # release step nodes excluded from each carrier's producers
end

"""Drive regulariser: rates are k·(u + ε), so with every drive of a pool at zero
the split tends to π0 rather than to whatever a constant added to Q selects
(review of specs/039, finding 5). ε is relative to a drive of 1 at baseline."""
const CYCLE_DRIVE_EPS = 1e-9

"""Baseline weight of a NON-enzyme exit from a state that also has an
enzyme-driven exit, relative to each enzyme-driven one (specs/039 amendment 2).
Intrinsic GTP hydrolysis and nucleotide exchange run orders of magnitude slower
than the GAP- or GEF-stimulated reactions, and basal (de)modification slower
than the enzyme's. Weighed per exit, so RAS's intrinsic hydrolysis competes
with the GAP path. Declared, not fitted."""
const CYCLE_INTRINSIC_WEIGHT = 1e-3

"""Paths longer than this are not resolved (the generator already caps at 6)."""
const CYCLE_MAX_STEPS = 6

"""
Stationary distribution of an n-state chain with rate matrix Q (rows sum to 0),
normalised; negative round-off clipped.
"""
function _stationary(Q::Matrix{Float64})
    n = size(Q, 1)
    A = Matrix(transpose(Q)); A[n, :] .= 1.0
    b = zeros(n); b[n] = 1.0
    pi = A \ b
    pi = max.(pi, 0.0)
    return pi ./ sum(pi)
end

"""
Resolve `pools` to solver indices. A pool is kept only if every state,
intermediate and step reaction has a node, every step reaction has an
IndexedReaction, every path is at most CYCLE_MAX_STEPS long, and the state graph
is connected; otherwise it is skipped and counted (reported as a fallback).
Carriers whose node or release node is missing are dropped from the pool.
"""
function index_pools(pools, uuid_to_idx::Dict{String, Int}, rxns_idx, phi::Float64)
    by_target = Dict{Int, Int}(r.target_idx => ri for (ri, r) in enumerate(rxns_idx))
    out = IndexedPool[]
    skipped = 0
    ratio = phi / (1.0 - phi)
    for p in pools
        nodes_ok = all(f -> haskey(uuid_to_idx, f), p.forms) && all(f -> haskey(uuid_to_idx, f), p.intermediates)
        if !nodes_ok || isempty(p.paths) || !(p.base in p.forms)
            skipped += 1; continue
        end
        pos = Dict(f => i for (i, f) in enumerate(p.forms))
        n = length(p.forms)
        # π0 over states: ρ^(distance from the base state) on the state graph
        adj = Dict(f => Set{String}() for f in p.forms)
        for q in p.paths
            (haskey(adj, q.from) && haskey(adj, q.to)) || continue
            push!(adj[q.from], q.to); push!(adj[q.to], q.from)
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
        if length(dist) != n
            skipped += 1; continue
        end
        w0 = [ratio ^ dist[f] for f in p.forms]
        pi0 = w0 ./ sum(w0)
        # paths
        pf = Int[]; pt = Int[]; psteps = Vector{Tuple{Int, Int}}[]; penz = Bool[]
        ok = true
        for q in p.paths
            (haskey(pos, q.from) && haskey(pos, q.to) && q.from != q.to &&
             1 <= length(q.steps) <= CYCLE_MAX_STEPS) || (ok = false; break)
            st = Tuple{Int, Int}[]
            for (a, b, rx) in q.steps
                (haskey(uuid_to_idx, a) && haskey(uuid_to_idx, rx)) || (ok = false; break)
                node = uuid_to_idx[rx]
                haskey(by_target, node) || (ok = false; break)
                push!(st, (by_target[node], uuid_to_idx[a]))
            end
            ok || break
            push!(pf, pos[q.from]); push!(pt, pos[q.to]); push!(psteps, st); push!(penz, q.enzyme)
        end
        (ok && !isempty(pf)) || (skipped += 1; continue)
        # baseline flow: the random walk on states, exits weighted by kind
        has_enz = falses(n)
        for q in eachindex(pf); penz[q] && (has_enz[pf[q]] = true); end
        w = [(!penz[q] && has_enz[pf[q]]) ? CYCLE_INTRINSIC_WEIGHT : 1.0 for q in eachindex(pf)]
        W = zeros(n)
        for q in eachindex(pf); W[pf[q]] += w[q]; end
        if any(==(0.0), W)
            skipped += 1; continue                 # a state with no exit: not one chain
        end
        P = zeros(n, n)
        for q in eachindex(pf); P[pf[q], pt[q]] += w[q] / W[pf[q]]; end
        Qw = P - Matrix{Float64}(I, n, n)
        mu = _stationary(Qw)
        if any(<=(0.0), mu)
            skipped += 1; continue                 # not strongly connected
        end
        j0 = [mu[pf[q]] * w[q] / W[pf[q]] for q in eachindex(pf)]
        k = [j0[q] / pi0[pf[q]] for q in eachindex(pf)]
        # flux nodes: step reaction nodes and intermediates
        fl = Dict{Int, Vector{Int}}()
        for (q, st) in enumerate(psteps)
            for (s_i, (ri, src)) in enumerate(st)
                push!(get!(fl, rxns_idx[ri].target_idx, Int[]), q)
                s_i > 1 && push!(get!(fl, src, Int[]), q)
            end
        end
        fnodes = sort!(collect(keys(fl)))
        fpaths = [unique(fl[v]) for v in fnodes]
        step_nodes = Set(rxns_idx[ri].target_idx for st in psteps for (ri, _) in st)
        supply = Vector{Vector{Int}}()
        for f in p.forms
            fi = uuid_to_idx[f]
            prod = Int[]
            if haskey(by_target, fi)
                for a in rxns_idx[by_target[fi]].activator_indices
                    a in step_nodes && continue
                    push!(prod, a)
                end
            end
            push!(supply, prod)
        end
        frx = [get(by_target, uuid_to_idx[f], 0) for f in p.forms]
        # carriers
        cs = Int[]; crx = Int[]; crel = Set{Int}[]
        cmap = Dict{Int, Set{Int}}()
        for (c, rel) in p.carriers
            (haskey(uuid_to_idx, c) && haskey(uuid_to_idx, rel)) || continue
            ci = uuid_to_idx[c]
            haskey(by_target, ci) || continue
            push!(get!(cmap, ci, Set{Int}()), uuid_to_idx[rel])
        end
        for ci in sort!(collect(keys(cmap)))
            push!(cs, ci); push!(crx, by_target[ci]); push!(crel, cmap[ci])
        end
        push!(out, IndexedPool([uuid_to_idx[f] for f in p.forms], pi0, pf, pt, k, j0, psteps,
                               fnodes, fpaths, supply, frx, cs, crx, crel))
    end
    return out, skipped
end

"""
Drive of each path: the product over its steps of the step's reaction, evaluated
with that step's source held at baseline, over the reaction's baseline.
"""
function path_drives(x::AbstractVector{Float64}, p::IndexedPool, rxns_idx, baseline_vec, config)
    u = Vector{Float64}(undef, length(p.path_from))
    for q in eachindex(p.path_from)
        d = 1.0
        for (ri, src) in p.path_steps[q]
            keep = x[src]
            x[src] = baseline_vec[src]
            out = compute_reaction_output_vec(x, rxns_idx[ri]; config = config)
            x[src] = keep
            bl = baseline_vec[rxns_idx[ri].target_idx]
            d *= bl > 0 ? out / bl : 0.0
            d == 0.0 && break
        end
        u[q] = d
    end
    return u
end

"""
Stationary distribution over the pool's states at the current drives, and the
drives themselves.
"""
function pool_stationary(x::AbstractVector{Float64}, p::IndexedPool, rxns_idx, baseline_vec, config)
    u = path_drives(x, p, rxns_idx, baseline_vec, config)
    n = length(p.forms)
    Q = zeros(n, n)
    for q in eachindex(p.path_from)
        i, j = p.path_from[q], p.path_to[q]
        rate = p.path_k[q] * (u[q] + CYCLE_DRIVE_EPS)
        Q[i, j] += rate; Q[i, i] -= rate
    end
    return _stationary(Q), u
end

"""
Fold applied to a node by its non-activator inputs (depletion and inhibitor
edges): its own update evaluated with every activator, and every node in
`extra_held`, at baseline, divided by its baseline. 1 when it has none.
"""
function _modifier(x::AbstractVector{Float64}, ri::Int, rxns_idx, baseline_vec, config, extra_held)
    ri == 0 && return 1.0
    r = rxns_idx[ri]
    (isempty(r.inhibitor_indices) && isempty(r.depletion_indices)) && return 1.0
    held = Set(r.activator_indices)
    union!(held, extra_held)
    held = sort!(collect(held))
    keep = [x[a] for a in held]
    for a in held; x[a] = baseline_vec[a]; end
    out = compute_reaction_output_vec(x, r; config = config)
    for (j, a) in enumerate(held); x[a] = keep[j]; end
    bl = baseline_vec[r.target_idx]
    return bl > 0 ? out / bl : 1.0
end

"""
Fold applied to a state by its non-activator inputs from outside the pool
(review of specs/039, finding 4). Sources that drive the pool's own steps (the
catalyst ⊣ source-form depletion edge) are also held at baseline: the balance
already contains them, as pre-registered. The depleted share is removed, not
redistributed to other states.
"""
function form_modifier(x::AbstractVector{Float64}, p::IndexedPool, k::Int, rxns_idx, baseline_vec, config)
    held = Set{Int}()
    for st in p.path_steps, (ri, _) in st
        union!(held, rxns_idx[ri].activator_indices)
    end
    return _modifier(x, p.form_rxn[k], rxns_idx, baseline_vec, config, held)
end

"""
Write the pool's states, step reaction nodes and intermediates from its steady
state, and (with `carriers`) its carriers. Pinned nodes keep their value; a
pinned state fixes the supply, otherwise the supply is the mean fold of the
states' producers outside the pool (1 if none). A path's flux fold is its drive
times its source state's fold; a flux node reads the baseline-flow-weighted mean
of the paths through it. A carrier reads the mean fold of its producers other
than the release steps (1 if none), times its modifier. Returns the largest
absolute change written.
"""
function update_pool!(x::Vector{Float64}, p::IndexedPool, rxns_idx, baseline_vec, obs_set, config;
                      carriers::Bool = true)
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
        m = form_modifier(x, p, k, rxns_idx, baseline_vec, config)
        v = clamp(baseline_vec[f] * s * m * pi[k] / p.pi0[k], 0.0, 1.0)
        maxch = max(maxch, abs(v - x[f])); x[f] = v
    end
    flux = [u[q] * x[p.forms[p.path_from[q]]] / baseline_vec[p.forms[p.path_from[q]]] for q in eachindex(u)]
    for (v, qs) in zip(p.flux_nodes, p.flux_paths)
        v in obs_set && continue
        num = 0.0; den = 0.0
        for q in qs
            num += p.path_j0[q] * flux[q]; den += p.path_j0[q]
        end
        val = clamp(baseline_vec[v] * (den > 0 ? num / den : 1.0), 0.0, 1.0)
        maxch = max(maxch, abs(val - x[v])); x[v] = val
    end
    if carriers
        for (c, ri, rel) in zip(p.carriers, p.carrier_rxn, p.carrier_release)
            c in obs_set && continue
            tot = 0.0; cnt = 0
            for a in rxns_idx[ri].activator_indices
                a in rel && continue
                bl = baseline_vec[a]
                bl > 0 || continue
                tot += x[a] / bl; cnt += 1
            end
            sc = cnt == 0 ? 1.0 : tot / cnt
            m = _modifier(x, ri, rxns_idx, baseline_vec, config, Set{Int}())
            val = clamp(baseline_vec[c] * sc * m, 0.0, 1.0)
            maxch = max(maxch, abs(val - x[c])); x[c] = val
        end
    end
    return maxch
end
