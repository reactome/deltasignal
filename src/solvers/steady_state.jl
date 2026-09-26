# Steady-State Solver Implementation

struct SolverResult
    node_activities::Dict{String, Float64}
    converged::Bool
    iterations::Int
    final_residual::Float64
    solve_time::Float64
    diagnostics::Dict{String, Any}
end

struct SteadyStateParams
    mu::Float64          # Model consistency penalty weight
    gamma::Float64       # Baseline prior weight  
    max_iters::Int       # Maximum optimization iterations
    tolerance::Float64   # Convergence tolerance
    method::String       # Solver method ("penalty" or "fixed_point")
end

function default_steady_state_params()
    # DS_MAX_ITERS overrides the per-SCC iteration budget. Added because the
    # budget was hardcoded with no way to distinguish "this component is
    # oscillating" from "this component is converging but needs more sweeps" —
    # a distinction that decides whether the fixed-point iteration can work at
    # all here. Weaker damping made non-convergence WORSE (122 -> 152 TP53
    # cases at lambda=0.1), which points at slow convergence rather than
    # oscillation, but only a larger budget settles it.
    max_iters = round(Int, _float_env("DS_MAX_ITERS", 500.0))
    if max_iters < 1
        throw(ArgumentError("DS_MAX_ITERS=$max_iters must be at least 1."))
    end
    # mu and gamma are read ONLY by DS_SCC_METHOD=minimize (specs/003). Under
    # the default fixed-point method they are still inert, which is why they
    # must not be reported as if they shaped the answer — see FR5.
    #
    # gamma's default is 1e-6, NOT the design doc's 0.1. Measured both ways
    # (specs/003 research.md section 6): gamma must clear ~5e-7 or the prior
    # cannot lift a collapsed loop out of the all-zero root (the escape
    # gradient 2*gamma*x0 falls below the optimiser's gradient tolerance), and
    # at the doc's 0.1 it competes with model consistency hard enough to read a
    # 100x perturbation as 51x. 1e-6 is three orders above the floor and five
    # below the doc.
    mu = _float_env("DS_MU", 1.0)
    gamma = _float_env("DS_GAMMA", 1e-6)
    if !(mu > 0.0) || !isfinite(mu)
        throw(ArgumentError("DS_MU=$mu must be finite and > 0."))
    end
    if gamma < 0.0 || !isfinite(gamma)
        throw(ArgumentError("DS_GAMMA=$gamma must be finite and >= 0."))
    end
    return SteadyStateParams(
        mu,
        gamma,
        max_iters,
        _float_env("DS_TOLERANCE", 1e-6),
        "penalty"
    )
end

"""
Confidence below this is treated as "no observation at all".

Shared by the initial guess and the pinning step. They were two separate
literals, and the initial guess did not apply the gate — which made it
inoperative for root nodes, the case the benchmark almost always exercises.
"""
const OBS_CONFIDENCE_TOL = 1e-6

"""
Solve steady-state network given observations.
Implements the penalty formulation:

minimize: Σᵢ ωᵢ(xᵢ - yᵢ)² + μ||x - F(x;θ)||² + γ||x - x₀||²

where:
- yᵢ are observed node activities
- F(x;θ) is the forward model
- x₀ are baseline activities
"""
function solve_steady_state(
    network::ReactionNetwork,
    observations::Dict{String, Tuple{Float64, Float64}},  # node_uuid -> (activity, confidence)
    params::SteadyStateParams = default_steady_state_params()
)::SolverResult
    
    start_time = time()
    
    # Cofactor handling is a MODELLING choice made here, not upstream: the
    # networks stay faithful to the curated data and this decides whether a
    # small molecule may carry a perturbation. See src/core/cofactors.jl.
    #
    # Deleting them instead was measured at -84 cases; pinning them gains +37.
    mode = cofactor_mode()
    cofactors = mode == "inert" ? cofactor_uuids(network) : Set{String}()

    # specs/032: drug-derived nodes held at baseline for a cell without the
    # drug. Same mechanism as cofactors, own switch and report.
    dmode = drug_mode()
    drugs = dmode == "inert" ? drug_uuids(network) : Set{String}()
    drug_rule = dmode == "propagate" ? "propagate" :
                network.drug_stids === nothing ? "inert: no drug table" : "inert"

    # Whether to traverse a positional-decomposition silo is a PROCESSING
    # decision, like the one above: the generator faithfully records that a
    # curated entity occurs in two places, and this decides whether a signal
    # reaching one occurrence is available at the other. Off by default; see
    # src/core/silo_bridges.jl for why it is capped.
    bridges = silo_bridge_edges(network)
    if !isempty(bridges)
        @info "Adding $(length(bridges)) silo bridge edge(s)" max_reach=silo_bridge_max_reach()
        network = ReactionNetwork(network.nodes, vcat(network.edges, bridges),
                                  network.set_mappings, network.cofactor_stids,
                                  network.containment, network.drug_stids)
    end

    # Convert network to reactions
    reactions = convert_to_reaction_network(network)

    # specs/022: inhibitors that contain their own reaction's input. Built unless
    # DS_SELF_INHIBITOR_WEIGHT=off.
    self_shared, self_rule = self_inhibitor_setup(network)
    
    # Initialize node activities
    all_nodes = collect(keys(network.nodes))
    n_nodes = length(all_nodes)
    node_to_idx = Dict(uuid => i for (i, uuid) in enumerate(all_nodes))
    
    if !isempty(cofactors)
        # Pinning at baseline is exactly "participant, not conduit": the node
        # still contributes its fold of 1.0 to every AND it belongs to, and
        # every reaction keeps the same input set, but its value never moves
        # so no perturbation can travel through it.
        #
        # An explicit observation always wins. Someone measuring an ATP
        # depletion is not making the modelling assumption this mode encodes,
        # and silently overwriting their input would be the same silent
        # substitution the config guard rails exist to prevent.
        # Annotate the type: when every cofactor is already observed this
        # comprehension is empty, which Julia infers as Dict{Any, Any}, and the
        # merge then widens `observations` away from the signature the solver
        # is dispatched on.
        pins = Dict{String, Tuple{Float64, Float64}}(
            u => (network.nodes[u].baseline * 100.0, 1.0)
            for u in cofactors if !haskey(observations, u))
        observations = merge(observations, pins)
    end
    drugs_held = 0
    if !isempty(drugs)
        # An observation (a measured drug, or a cofactor pin) wins; see above.
        dpins = Dict{String, Tuple{Float64, Float64}}(
            u => (network.nodes[u].baseline * 100.0, 1.0)
            for u in drugs if !haskey(observations, u))
        drugs_held = length(dpins)
        observations = merge(observations, dpins)
    end

    # Initial guess: use observations where available, baseline elsewhere
    x0 = Dict{String, Float64}()
    baseline_activities = Dict{String, Float64}()
    
    # Gene-entity UUIDs (nodes whose Reactome stable_id is a gene) — used to
    # detect transcriptional autoregulation loops. Empty unless DS_GENE_STIDS_FILE
    # is set, so this is a no-op by default.
    gene_stids = gene_stid_set()
    gene_uuids = Set{String}()
    for (uuid, node) in network.nodes
        # Seed from an observation ONLY if it clears the same confidence gate
        # the pinning step applies. Seeding an unpinned observation silently
        # defeated that gate: a ROOT node has no incoming reaction, so nothing
        # ever overwrites its initial value and a confidence-0 observation was
        # applied exactly as hard as a confidence-1 one. Non-root nodes hid the
        # bug, because the forward model overwrote them on the first sweep.
        if haskey(observations, uuid) && observations[uuid][2] > OBS_CONFIDENCE_TOL
            x0[uuid] = observations[uuid][1] / 100.0  # Convert from UI scale 0-100 to internal 0-1
        else
            x0[uuid] = node.baseline
        end
        baseline_activities[uuid] = node.baseline
        if !isempty(gene_stids) && node.reactome_id !== nothing && node.reactome_id in gene_stids
            push!(gene_uuids, uuid)
        end
    end

    # The penalty/SCC-condensation solver is the only path (see
    # solve_steady_state_penalty). `params.method` is retained for API
    # compatibility but only "penalty" is supported.
    params.method == "penalty" ||
        error("Unsupported solver method $(params.method); only \"penalty\" is available.")
    result = solve_steady_state_penalty(reactions, observations, x0, baseline_activities, params, start_time, gene_uuids;
                                        self_shared = self_shared)
    # Say which model ran: "on", "off", or inert because nothing told the solver
    # what contains what (a POSTed network or an older bundle).
    result.diagnostics["self_inhibitor_rule"] = self_rule
    result.diagnostics["drug_rule"] = drug_rule
    result.diagnostics["drugs_held"] = drugs_held
    return result
end

"""
(map, status) for DS_SELF_INHIBITOR_WEIGHT on this network. status is "off",
"on", or "inert: no containment table" -- the last means the default model is
NOT the one being run, which a caller must be able to see.
"""
function self_inhibitor_setup(network::ReactionNetwork)
    empty = Dict{Tuple{String, String}, Vector{String}}()
    _self_inhibitor_weight_env() < 0 && return empty, "off"
    isempty(network.containment) && return empty, "inert: no containment table"
    return self_contained_inhibitor_map(network), "on"
end

"""
(inhibitor uuid, target uuid) -> the target's activator uuids that the
inhibitor CONTAINS, according to the network's containment table (specs/022).
A node's identity is its own `reactome_id`; activators are collected per target
node, so an inhibitor is never matched against a sibling variant's input.
Depletion edges are not inhibitors here and are never included.
"""
function self_contained_inhibitor_map(network::ReactionNetwork)::Dict{Tuple{String, String}, Vector{String}}
    out = Dict{Tuple{String, String}, Vector{String}}()
    isempty(network.containment) && return out
    sid(u) = (n = get(network.nodes, u, nothing); n === nothing ? nothing : n.reactome_id)
    acts = Dict{String, Vector{String}}()
    for e in network.edges
        e.is_positive && push!(get!(acts, e.child_uuid, String[]), e.parent_uuid)
    end
    # Containment alone is not enough: the rule divides the shared input's fold
    # out of the inhibitor, which is only right if the inhibitor is actually
    # COMPUTED from that input here. 71 of 625 catalog pairs name an input that
    # cannot reach the inhibitor at all (all 48 in Signaling by ERBB2).
    fwd = Dict{String, Vector{String}}()
    for e in network.edges
        push!(get!(fwd, e.parent_uuid, String[]), e.child_uuid)
    end
    reach_cache = Dict{String, Set{String}}()
    function reaches(a::String)
        get!(reach_cache, a) do
            seen = Set{String}(); stack = copy(get(fwd, a, String[]))
            while !isempty(stack)
                n = pop!(stack)
                n in seen && continue
                push!(seen, n)
                append!(stack, get(fwd, n, String[]))
            end
            seen
        end
    end
    for e in network.edges
        (e.is_positive || e.edge_type == "depletion") && continue
        s = sid(e.parent_uuid)
        s === nothing && continue
        inside = get(network.containment, s, nothing)
        inside === nothing && continue
        shared = String[a for a in unique(get(acts, e.child_uuid, String[]))
                        if (as = sid(a)) !== nothing && as in inside && e.parent_uuid in reaches(a)]
        isempty(shared) || (out[(e.parent_uuid, e.child_uuid)] = sort(shared))
    end
    return out
end

"""
SCC-condensation solve. Processes strongly-connected components in topological
order (Tarjan numbers them in reverse-topo order, so DESCENDING comp_id is
upstream-first). Singleton (acyclic) components are evaluated exactly in one
pass — Gauss-Seidel reading the already-final upstream values — so the acyclic
majority of the network (≈90%+) is solved without any damping and matches the
plain feed-forward result. Non-trivial components (genuine loops) get a damped
fixed-point iteration CONFINED to that component, with external inputs held
fixed at their solved upstream values. This makes feedback loops converge to a
stable point instead of the flat solver's oscillating last-iterate.

Mutates `x` in place. Returns (total_inner_iters, last_max_change). Observations
in `obs_set` are never updated (hard-pinned). `DS_SCC_DAMPING` (default 0.5)
sets the per-component blend; `DS_SCC_BREAK_CATALYST` additionally freezes
recycling-catalyst edges at their component-entry value (conserved moiety).
"""

"""
A full-length state vector with a handful of entries overlaid by dual numbers.

`component_gradient!` differentiates one reaction at a time with respect to
that reaction's own inputs. Materialising a full-length `Vector{Dual}` for each
one allocates the whole network per reaction per gradient — ~2M allocations per
gradient on a 840-reaction component, which dominated everything else. The
reaction reads only its own few entries, so this overlays them on the Float64
state instead: O(1) construction, O(#overlaid) lookup, nothing allocated.

`idx` is short (a reaction's input count), so the linear scan beats a Dict.
"""
struct OverlayVector{T, S} <: AbstractVector{T}
    base::Vector{S}
    idx::Vector{Int}
    vals::Vector{T}
end

Base.size(v::OverlayVector) = size(v.base)
Base.length(v::OverlayVector) = length(v.base)
Base.IndexStyle(::Type{<:OverlayVector}) = IndexLinear()

Base.@propagate_inbounds function Base.getindex(v::OverlayVector{T, S}, i::Int) where {T, S}
    @inbounds for k in eachindex(v.idx)
        v.idx[k] == i && return v.vals[k]
    end
    return convert(T, @inbounds v.base[i])
end

"""
One cyclic component's minimisation problem, in the form the optimiser needs.

`free` are the component's unobserved target nodes — the variables. `xbuf` is a
full-length state vector whose entries OUTSIDE the component stay at their final
upstream values; the free entries are overwritten per evaluation. `deps[k]` is
the set of free variables reaction `rs[k]` actually reads, which is what makes
the gradient sparse.
"""
struct ComponentProblem
    rs::Vector{Int}
    rxns::Vector{IndexedReaction}
    free::Vector{Int}
    slot::Dict{Int, Int}
    deps::Vector{Vector{Int}}
    xbuf::Vector{Float64}
    baseline::Vector{Float64}
    mu::Float64
    gamma::Float64
    config::ReactionEvalConfig
    supply::Union{Nothing, Vector{Float64}}
end

function ComponentProblem(x, rs, rxns_idx, obs_set, baseline_vec,
                          mu, gamma, config, supply)
    free = Int[]
    @inbounds for ri in rs
        t = rxns_idx[ri].target_idx
        t in obs_set && continue
        push!(free, t)
    end
    slot = Dict{Int, Int}(t => k for (k, t) in enumerate(free))
    deps = Vector{Vector{Int}}(undef, length(rs))
    @inbounds for (k, ri) in enumerate(rs)
        r = rxns_idx[ri]
        d = Int[]
        append!(d, r.activator_indices); append!(d, r.inhibitor_indices)
        append!(d, r.depletion_indices); append!(d, r.substrate_indices)
        push!(d, r.target_idx)
        deps[k] = unique!(sort!(filter(i -> haskey(slot, i), d)))
    end
    ComponentProblem(collect(rs), rxns_idx, free, slot, deps, copy(x),
                     baseline_vec, mu, gamma, config, supply)
end

"""`L(z) = mu*||F(x) - x||^2 + gamma*||x - x0||^2` restricted to the component."""
function component_objective(p::ComponentProblem, z::Vector{Float64})::Float64
    @inbounds for (k, t) in enumerate(p.free)
        p.xbuf[t] = z[k]
    end
    acc = 0.0
    @inbounds for ri in p.rs
        r = p.rxns[ri]
        t = r.target_idx
        haskey(p.slot, t) || continue
        d = compute_reaction_output_vec(p.xbuf, r; supply=p.supply, config=p.config) - p.xbuf[t]
        acc += d * d
    end
    prior = 0.0
    @inbounds for (k, t) in enumerate(p.free)
        d = z[k] - p.baseline[t]
        prior += d * d
    end
    p.mu * acc + p.gamma * prior
end

"""
Gradient of `component_objective`, assembled from per-reaction pieces.

Differentiating the whole objective at once with forward-mode AD costs O(n)
evaluations per gradient, and a component here reaches 836 nodes: measured at
~400x the fixed point's cost by 80 nodes, and one catalog solve ran 80 minutes
without returning. But the objective is a sum of per-reaction terms and each
reaction reads only its own few inputs, so

    dL/dx_j = 2*mu * SUM_r (F_r - x_tr) * (dF_r/dx_j - delta_{j,tr})
              + 2*gamma * (x_j - x0_j)

is built by taking a SMALL ForwardDiff gradient per reaction over that
reaction's own inputs and scattering it. Cost is O(sum of in-degrees) =
O(edges), independent of component size, and adds no dependency.
"""
function component_gradient!(G::Vector{Float64}, p::ComponentProblem, z::Vector{Float64})
    @inbounds for i in eachindex(G); G[i] = 0.0; end
    @inbounds for (k, t) in enumerate(p.free); p.xbuf[t] = z[k]; end
    @inbounds for (k, ri) in enumerate(p.rs)
        r = p.rxns[ri]
        t = r.target_idx
        haskey(p.slot, t) || continue
        d = p.deps[k]
        fval = compute_reaction_output_vec(p.xbuf, r; supply=p.supply, config=p.config)
        resid = fval - p.xbuf[t]
        if !isempty(d)
            local_in = Float64[p.xbuf[i] for i in d]
            fr = function (v)
                yy = OverlayVector{eltype(v), Float64}(p.xbuf, d, v)
                compute_reaction_output_vec(yy, r; supply=p.supply, config=p.config)
            end
            gloc = ForwardDiff.gradient(fr, local_in)
            @inbounds for (m, i) in enumerate(d)
                G[p.slot[i]] += 2.0 * p.mu * resid * gloc[m]
            end
        end
        G[p.slot[t]] -= 2.0 * p.mu * resid
    end
    @inbounds for (k, t) in enumerate(p.free)
        G[k] += 2.0 * p.gamma * (z[k] - p.baseline[t])
    end
    G
end


"""
Residual vector and its sparse Jacobian for one component, in least-squares form.

    r = [ sqrt(mu) * (F_r(x) - x_tr)   for each reaction r in the component
          sqrt(gamma) * (x_j - x0_j)   for each free node j ]

so that `L = ||r||^2` is exactly `component_objective`. Returning `r` and `J`
separately is the whole point: the objective is a sum of squares, so its
Hessian is well approximated by `2*J'J` and a Gauss-Newton / Levenberg-Marquardt
step uses that directly. LBFGS only sees the scalar and has to rebuild curvature
from gradient history, which is why it needed hundreds of iterations and was
still drifting at 100 on an 840-node component.

`J` is very sparse: a reaction row has one entry per input it reads plus its own
target, and each prior row is a single diagonal entry.
"""
function component_residual_jacobian(p::ComponentProblem, z::Vector{Float64};
                                     gamma::Float64 = p.gamma)
    nfree = length(p.free)
    nrx = 0
    @inbounds for ri in p.rs
        haskey(p.slot, p.rxns[ri].target_idx) && (nrx += 1)
    end
    m = nrx + nfree
    r = zeros(Float64, m)
    I = Int[]; Jc = Int[]; V = Float64[]
    sm = sqrt(p.mu); sg = sqrt(gamma)

    @inbounds for (k, t) in enumerate(p.free); p.xbuf[t] = z[k]; end

    row = 0
    @inbounds for (k, ri) in enumerate(p.rs)
        rx = p.rxns[ri]
        t = rx.target_idx
        haskey(p.slot, t) || continue
        row += 1
        d = p.deps[k]
        fval = compute_reaction_output_vec(p.xbuf, rx; supply=p.supply, config=p.config)
        r[row] = sm * (fval - p.xbuf[t])
        if !isempty(d)
            local_in = Float64[p.xbuf[i] for i in d]
            fr = function (v)
                yy = OverlayVector{eltype(v), Float64}(p.xbuf, d, v)
                compute_reaction_output_vec(yy, rx; supply=p.supply, config=p.config)
            end
            gloc = ForwardDiff.gradient(fr, local_in)
            for (mm, i) in enumerate(d)
                push!(I, row); push!(Jc, p.slot[i]); push!(V, sm * gloc[mm])
            end
        end
        # d(F_r - x_t)/dx_t also carries the -1 from the subtracted target.
        push!(I, row); push!(Jc, p.slot[t]); push!(V, -sm)
    end
    @inbounds for k in 1:nfree
        row += 1
        r[row] = sg * (z[k] - p.baseline[p.free[k]])
        push!(I, row); push!(Jc, k); push!(V, sg)
    end
    J = sparse(I, Jc, V, m, nfree)
    return r, J
end

"""
Levenberg-Marquardt minimisation of `component_objective`, box-clamped to [0,1].

Solves `(J'J + lambda*diag(J'J)) delta = -J'r` each step, growing `lambda` when a
step fails and shrinking it when it succeeds. `J'J` is sparse and factorised
directly.

Returns `(converged, iterations)`.
"""
function lm_minimize!(z::Vector{Float64}, p::ComponentProblem;
                      max_iter::Int = 100, g_tol::Float64 = 1e-10,
                      step_tol::Float64 = 1e-12, gamma::Float64 = p.gamma)
    lambda = 1e-3
    r, J = component_residual_jacobian(p, z; gamma=gamma)
    cost = dot(r, r)
    local_iters = 0
    for it in 1:max_iter
        local_iters = it
        JtJ = Symmetric(Matrix(J' * J))
        g = J' * r
        norm(g, Inf) < g_tol && return (true, it)
        stepped = false
        for _ in 1:12          # lambda back-off within one iteration
            A = Matrix(JtJ) + lambda * Diagonal(max.(diag(JtJ), 1e-12))
            local delta
            try
                delta = -(A \ g)
            catch
                lambda *= 10.0
                continue
            end
            znew = clamp.(z .+ delta, 0.0, 1.0)
            rnew, Jnew = component_residual_jacobian(p, znew; gamma=gamma)
            cnew = dot(rnew, rnew)
            if cnew < cost
                if norm(znew .- z, Inf) < step_tol
                    z .= znew
                    return (true, it)
                end
                z .= znew; r = rnew; J = Jnew; cost = cnew
                lambda = max(lambda * 0.3, 1e-12)
                stepped = true
                break
            else
                lambda *= 10.0
                lambda > 1e12 && break
            end
        end
        stepped || return (false, it)  # no downhill step found (lambda exhausted / solve failed): NOT a certified minimum
    end
    return (false, local_iters)
end

"""
Minimise the specified objective over one cyclic component's free nodes,
holding everything outside the component fixed.

Returns the component's final model residual `max|F(x) - x|`, on the same scale
the fixed-point path reports, so the convergence verdict keeps its meaning.

Observed nodes are excluded from the free set rather than penalised, which is
the current hard-constraint semantics (FR4 makes them weighted later).
"""
function minimize_component!(
    x::Vector{Float64},
    rs::Vector{Int},
    rxns_idx::Vector{IndexedReaction},
    obs_set::Set{Int},
    baseline_vec::Vector{Float64},
    params::SteadyStateParams,
    config::ReactionEvalConfig,
    supply::Union{Nothing,Vector{Float64}},
)::Float64
    prob = ComponentProblem(x, rs, rxns_idx, obs_set, baseline_vec,
                            params.mu, params.gamma, config, supply)
    isempty(prob.free) && return 0.0

    z0 = Float64[x[t] for t in prob.free]

    # Optimiser. LM is the default because the objective is a sum of squares:
    # measured on an 840-node component, LBFGS was still drifting at 100
    # iterations (readout 1.01 -> 1.31 -> 1.12 -> 0.22 at 5/10/25/100) and took
    # 29s, because it rebuilds curvature from gradient history instead of using
    # the J'J that a least-squares problem hands you for free.
    opt = get(ENV, "DS_SCC_OPTIMIZER", "lm")
    if !(opt in ("lm", "lbfgs"))
        throw(ArgumentError(
            "DS_SCC_OPTIMIZER=$opt is not an optimiser; expected \"lm\" or \"lbfgs\"."
        ))
    end

    try
        if opt == "lm"
            zb = copy(z0)
            cap = max(1, min(params.max_iters, _lm_iter_cap()))
            gtol = params.tolerance * 1e-2
            # gamma-continuation (homotopy). A cyclic component can have several
            # self-consistent states -- on an 840-node ring, collapsed (~0.2 UI),
            # baseline (1.0) and saturated (100) are all roots -- and a plain
            # solve lands in whichever one the warm start sits nearest. That
            # makes the answer a property of where we happened to start.
            #
            # Instead, start at a LARGE gamma, where the prior dominates and the
            # minimiser is unambiguous (near baseline), then step gamma down,
            # warm-starting each solve from the last. This tracks the branch
            # CONNECTED TO BASELINE down to the target gamma rather than picking
            # a root by accident. Off by default; DS_GAMMA_ANNEAL=n sets steps.
            for g in _gamma_schedule(params.gamma)
                lm_minimize!(zb, prob; max_iter=cap, g_tol=gtol, gamma=g)
            end
            @inbounds for (k, t) in enumerate(prob.free)
                x[t] = clamp(zb[k], 0.0, 1.0)
            end
            return _component_residual(x, rs, rxns_idx, obs_set, supply, config)
        end
        lower = zeros(Float64, length(prob.free))
        upper = ones(Float64, length(prob.free))
        res = Optim.optimize(
            z -> component_objective(prob, z),
            (G, z) -> component_gradient!(G, prob, z),
            lower, upper, z0,
            Optim.Fminbox(Optim.LBFGS()),
            Optim.Options(iterations = params.max_iters,
                          g_tol = params.tolerance),
        )
        zbest = Optim.minimizer(res)
        @inbounds for (k, t) in enumerate(prob.free)
            x[t] = clamp(zbest[k], 0.0, 1.0)
        end
    catch err
        # A failed minimisation must not hand back the warm start as though it
        # were a minimiser. Keep it (it is a valid fixed-point answer) but say so.
        @warn "component minimisation failed; keeping the fixed-point warm start" exception=(err, catch_backtrace()) n_free=length(prob.free)
    end

    return _component_residual(x, rs, rxns_idx, obs_set, supply, config)
end

"""Worst `|F(x) - x|` over a component's unobserved targets."""
function _component_residual(x, rs, rxns_idx, obs_set, supply, config)::Float64
    resid = 0.0
    @inbounds for ri in rs
        r = rxns_idx[ri]
        t = r.target_idx
        t in obs_set && continue
        d = abs(compute_reaction_output_vec(x, r; supply=supply, config=config) - x[t])
        d > resid && (resid = d)
    end
    return resid
end


"""
Decreasing gamma schedule ending at the target, for continuation.

`DS_GAMMA_ANNEAL=0` (the default) returns just the target, i.e. no continuation.
`n > 0` prepends `n` geometrically-spaced larger values starting from
`DS_GAMMA_ANNEAL_START` (default 1.0).
"""
function _gamma_schedule(target::Float64)
    n = round(Int, _float_env("DS_GAMMA_ANNEAL", 0.0))
    n < 0 && throw(ArgumentError("DS_GAMMA_ANNEAL=$n must be >= 0."))
    (n == 0 || target <= 0.0) && return [target]
    start = _float_env("DS_GAMMA_ANNEAL_START", 1.0)
    start <= target && return [target]
    ratio = (target / start)^(1.0 / n)
    sched = [start * ratio^(k - 1) for k in 1:n]
    push!(sched, target)
    sched
end

"""LM iteration cap. Far smaller than the fixed-point budget by design."""
function _lm_iter_cap()
    n = round(Int, _float_env("DS_LM_ITERS", 60.0))
    n < 1 && throw(ArgumentError("DS_LM_ITERS=$n must be at least 1."))
    n
end

function solve_scc_ordered!(
    x::Vector{Float64},
    rxns_idx::Vector{IndexedReaction},
    comp_id::Vector{Int},
    n_comp::Int,
    obs_set::Set{Int},
    params::SteadyStateParams,
    config::ReactionEvalConfig = resolve_reaction_eval_config(),
    baseline_vec::Union{Nothing,Vector{Float64}} = nothing,
)
    # Damping must lie in (0, 1]. The update nv = (1-λ)·x + λ·F(x) is written
    # back unclamped, so λ outside [0,1] extrapolates past the model's [0,1]
    # domain and reaches log((v+ε)/(bl+ε)) as a negative argument — an uncaught
    # DomainError out of the solver. λ = 0 never updates the component at all
    # and silently burns every iteration. Both are sweep typos, not choices.
    λ = _float_env("DS_SCC_DAMPING", 0.5)
    if !(0.0 < λ <= 1.0)
        throw(ArgumentError(
            "DS_SCC_DAMPING=$λ is out of range; must be in (0, 1]. " *
            "λ ≤ 0 never updates the component; λ > 1 extrapolates outside " *
            "the model domain and throws from inside the propagator."
        ))
    end
    # `supply` (component-entry state) is read by recycling-closure edges: the
    # legacy catalyst knob and the role list of specs/018.
    use_supply = _bool_env("DS_SCC_BREAK_CATALYST", false) || !isempty(_break_roles_env())
    max_inner = params.max_iters
    tol = params.tolerance

    # Sweep scheme inside a cyclic component.
    #
    #   "gauss_seidel" (default, historical): x[t] is written back inside the
    #       sweep, so each reaction reads values some of its neighbours have
    #       already updated this pass. Converges faster, but the result depends
    #       on the ORDER reactions are visited in. That order comes from Julia
    #       `Dict` iteration in `convert_to_reaction_network` and
    #       `solve_steady_state`, i.e. from UUID hashing — so relabelling a
    #       network's nodes, with the graph otherwise identical, can land the
    #       component in a different fixed point. Measured: renaming every UUID
    #       in the 92-pathway catalog (verified isomorphic — identical
    #       stable-id edge multiset for all 92) moved 14 of 23,908 curator
    #       predictions, every one of them inside a single cyclic pathway.
    #       No acyclic pathway moved a case.
    #
    #   "jacobi": every reaction is evaluated against the state at the START of
    #       the sweep and the new values are committed together, so the sweep is
    #       invariant to the order of `rs` and the solve becomes a function of
    #       the graph rather than of its node names.
    #
    # This does not resolve WHY a cyclic component has more than one fixed point
    # (that is the all-zero root, specs/004); it removes node labelling as the
    # thing that picks between them.
    sweep = get(ENV, "DS_SCC_SWEEP", "gauss_seidel")
    if !(sweep in ("gauss_seidel", "jacobi"))
        throw(ArgumentError(
            "DS_SCC_SWEEP=$sweep is not a sweep scheme; " *
            "expected \"gauss_seidel\" or \"jacobi\"."
        ))
    end
    jacobi = sweep == "jacobi"

    # How a CYCLIC component is resolved. specs/003-solver-objective.
    #
    #   "fixed_point" (default, historical): sweep x <- F(x) until the residual
    #       stops moving. On a component with several self-consistent states
    #       this finds *a* root, not a defined one, and which root depends on
    #       the sweep order. Raising the iteration budget does not help: at 20x
    #       the budget two labellings of one network diverged FURTHER (14 -> 25
    #       differing predictions). See specs/003 research.md section 3.
    #
    #   "minimize": minimise the objective the design specifies over the
    #       component's free nodes,
    #
    #           L = mu*||F(x) - x||^2 + gamma*||x - x0||^2
    #
    #       box-constrained to [0,1], LBFGS with ForwardDiff gradients, warm
    #       started from the fixed-point result. The gamma term is the point:
    #       with gamma = 0 the all-zero state is a stationary point AND a
    #       global minimum (gradient exactly 0, L = 0), tied with the correct
    #       answer, so nothing prefers the right root -- which is the defect
    #       described in review. With gamma = 0.1 the correct state scores ~1e28
    #       better. See specs/003 research.md section 2.
    #
    # Observations stay pinned here; making them weighted is FR4, Stage 3.
    #   "pool" / "pool_parity" (specs/017): the component is a conserved POOL.
    #       Every member reaction is evaluated with its in-component inputs held
    #       at baseline (fold 1), giving an EXTERNAL fold; the entry folds
    #       multiply into one pool fold (0 absorbs, capped at 100x like hill_sat)
    #       and every non-pinned member reads its own baseline x that fold. No
    #       iteration inside the component, so the gain-1 knife-edge cannot rail
    #       or collapse it and the answer does not depend on node labels.
    #       "pool" leaves any component with an internal inhibitor/depletion
    #       edge to the damped fixed point; "pool_parity" reads each entry's
    #       fold to the power +-1 by the parity of negative internal edges on
    #       the path from that entry, and falls back only when the signs are
    #       inconsistent (an odd negative cycle: genuine negative feedback).
    #       "pool_all" is the literal rule: pool every component, and an
    #       internal negative edge contributes nothing (its source reads
    #       baseline, factor 1). A member downstream of an internal inhibitor
    #       therefore moves WITH the entry, not against it -- known wrong for
    #       MDM2 -| TP53, kept because it is the only variant that pools the
    #       giant components at all (every one of them has a negative edge).
    scc_method = get(ENV, "DS_SCC_METHOD", "fixed_point")
    if !(scc_method in ("fixed_point", "minimize", "pool", "pool_parity", "pool_all"))
        throw(ArgumentError(
            "DS_SCC_METHOD=$scc_method is not a method; " *
            "expected \"fixed_point\", \"minimize\", \"pool\", \"pool_parity\" or \"pool_all\"."
        ))
    end
    minimize = scc_method == "minimize"
    pooling = scc_method in ("pool", "pool_parity", "pool_all")
    parity = scc_method == "pool_parity"
    pool_all = scc_method == "pool_all"
    if pooling && baseline_vec === nothing
        throw(ArgumentError(
            "DS_SCC_METHOD=$scc_method needs baselines: solve_scc_ordered! was " *
            "called without `baseline_vec`. This is a caller bug, not a configuration one."
        ))
    end
    if minimize && baseline_vec === nothing
        throw(ArgumentError(
            "DS_SCC_METHOD=minimize needs baselines: solve_scc_ordered! was " *
            "called without `baseline_vec`, so the gamma||x - x0||^2 term " *
            "cannot be formed. This is a caller bug, not a configuration one."
        ))
    end

    # Staging buffers for the Jacobi scheme, allocated ONCE and grown to the
    # largest component rather than per component — the Gauss-Seidel path must
    # not pay an allocation per SCC for a scheme it does not use. There is one
    # reaction per target node (`convert_to_reaction_network` groups every edge
    # by child), so staged writes never collide.
    stage_idx = Int[]
    stage_val = Float64[]

    # Type-aware loop handling. Positive/recycling SCCs want full convergence
    # (resolves spuriously-sustained signal). Negative-feedback / switch SCCs
    # are hurt by full relaxation — the fully-relaxed steady state buffers or
    # inverts the transient response a perturbation experiment measures. For
    # those we cap the iteration at DS_SCC_NEG_ITERS sweeps ("transient" read)
    # instead of relaxing to the feedback-corrected fixed point.
    #   DS_SCC_NEG_MODE: "converge" (default, plain v2) | "transient"
    #   DS_SCC_NEG_FRAC: fraction of internal edges that must be negative for an
    #                    SCC to count as negative-feedback-dominant (default 0.5)
    #   DS_SCC_NEG_ITERS: sweeps for negative SCCs in transient mode (default 1)
    neg_mode = get(ENV, "DS_SCC_NEG_MODE", "converge")
    neg_frac_thresh = _float_env("DS_SCC_NEG_FRAC", 0.5)
    neg_iters = parse(Int, get(ENV, "DS_SCC_NEG_ITERS", "1"))

    # Reactions grouped by their target node's component.
    comp_rxns = [Int[] for _ in 1:n_comp]
    @inbounds for ri in eachindex(rxns_idx)
        c = comp_id[rxns_idx[ri].target_idx]
        push!(comp_rxns[c], ri)
    end
    # Node count per component → distinguishes acyclic singletons from real SCCs.
    comp_size = zeros(Int, n_comp)
    @inbounds for c in comp_id
        c >= 1 && (comp_size[c] += 1)
    end

    # Components containing a SELF-loop must be iterated, not evaluated once.
    # Tarjan correctly assigns a self-looping node its own single-node
    # component, so `comp_size == 1` alone would classify it as acyclic and
    # evaluate it exactly once — reading its own pre-update value and reporting
    # a converged solve with residual 0 (measured 8.25x off the damped
    # fixed-point answer on a 3-node network). A node that feeds itself is a
    # cycle regardless of component size.
    comp_has_self_loop = falses(n_comp)
    @inbounds for r in rxns_idx
        t = r.target_idx
        c = comp_id[t]
        c >= 1 || continue
        comp_has_self_loop[c] && continue
        if t in r.activator_indices || t in r.inhibitor_indices || t in r.depletion_indices ||
           t in r.substrate_indices
            comp_has_self_loop[c] = true
        end
    end

    # Per-component internal edge sign census (only when transient mode is on).
    # Activators count as positive; inhibitors + depletions as negative. Only
    # intra-SCC edges (source and target in the same component) are counted.
    comp_neg_frac = zeros(Float64, n_comp)
    if neg_mode != "converge"
        comp_pos = zeros(Int, n_comp); comp_neg = zeros(Int, n_comp)
        @inbounds for r in rxns_idx
            c = comp_id[r.target_idx]
            comp_size[c] > 1 || continue
            for a in r.activator_indices
                comp_id[a] == c && (comp_pos[c] += 1)
            end
            for i in r.inhibitor_indices
                comp_id[i] == c && (comp_neg[c] += 1)
            end
            for dpt in r.depletion_indices
                comp_id[dpt] == c && (comp_neg[c] += 1)
            end
        end
        @inbounds for c in 1:n_comp
            tot = comp_pos[c] + comp_neg[c]
            comp_neg_frac[c] = tot > 0 ? comp_neg[c] / tot : 0.0
        end
    end

    # Pool bookkeeping (specs/017): member node lists and an internal-negative
    # census per component, built only when pooling is on.
    comp_nodes = pooling ? [Int[] for _ in 1:n_comp] : Vector{Int}[]
    comp_has_neg = falses(n_comp)
    if pooling
        @inbounds for (i, c) in enumerate(comp_id)
            c >= 1 && push!(comp_nodes[c], i)
        end
        @inbounds for r in rxns_idx
            c = comp_id[r.target_idx]
            c >= 1 || continue
            comp_has_neg[c] && continue
            if any(i -> comp_id[i] == c, r.inhibitor_indices) ||
               any(i -> comp_id[i] == c, r.depletion_indices)
                comp_has_neg[c] = true
            end
        end
    end
    n_pooled = 0; n_iterated = 0; n_fb_neg = 0; n_fb_inc = 0; n_pooled_nodes = 0

    total_iters = 0
    last_change = 0.0

    # Topological order: upstream (high comp_id) before downstream (low).
    @inbounds for c in n_comp:-1:1
        rs = comp_rxns[c]
        isempty(rs) && continue

        if comp_size[c] == 1 && !comp_has_self_loop[c]
            # Acyclic node: exact single evaluation (upstream already final).
            for ri in rs
                t = rxns_idx[ri].target_idx
                t in obs_set && continue
                # A closure edge into an acyclic node reads the current state,
                # which IS its entry value (upstream is final, downstream is
                # still at its initial value).
                x[t] = compute_reaction_output_vec(x, rxns_idx[ri]; supply = use_supply ? x : nothing, config=config)
            end
        else
            if pooling
                status, nw = pool_component!(x, rs, rxns_idx, comp_id, c, comp_nodes[c],
                                             obs_set, baseline_vec, config, parity,
                                             pool_all ? false : comp_has_neg[c])
                if status == :pooled
                    n_pooled += 1
                    n_pooled_nodes += nw
                    continue
                elseif status == :negative
                    n_fb_neg += 1
                else
                    n_fb_inc += 1
                end
            end
            n_iterated += 1
            # Genuine loop: damped fixed point confined to this component.
            # Freeze the entry state for the optional catalyst-break layer.
            supply = use_supply ? copy(x) : nothing
            # Negative-feedback SCCs get a bounded (transient) relaxation; all
            # others relax to convergence.
            is_neg = neg_mode != "converge" && comp_neg_frac[c] >= neg_frac_thresh
            cap = is_neg ? neg_iters : max_inner
            comp_residual = 0.0
            if jacobi && length(rs) > length(stage_idx)
                resize!(stage_idx, length(rs))
                resize!(stage_val, length(rs))
            end
            for it in 1:cap
                total_iters += 1
                maxch = 0.0
                n_staged = 0
                for ri in rs
                    r = rxns_idx[ri]
                    t = r.target_idx
                    t in obs_set && continue
                    fwd = compute_reaction_output_vec(x, r; supply=supply, config=config)
                    nv = (1.0 - λ) * x[t] + λ * fwd
                    # Measure the UNDAMPED model residual |F(x) - x|, not the
                    # damped step |nv - x|. The damped step is λ·|F(x) - x|, so
                    # with the default λ = 0.5 a component could break at a true
                    # residual of up to 2·tolerance and then be reported
                    # non-converged by the |F(x) - x| < tolerance test below —
                    # the stopping rule and the convergence verdict were on
                    # different scales.
                    ch = abs(fwd - x[t])
                    ch > maxch && (maxch = ch)
                    if jacobi
                        n_staged += 1
                        stage_idx[n_staged] = t
                        stage_val[n_staged] = nv
                    else
                        x[t] = nv
                    end
                end
                if jacobi
                    for s in 1:n_staged
                        x[stage_idx[s]] = stage_val[s]
                    end
                end
                comp_residual = maxch
                maxch < tol && break
            end
            # Warm-started minimisation of the specified objective. The
            # fixed-point sweep above produced the starting point; this moves
            # from "a root the sweep happened to reach" to "the minimiser of
            # the objective", which is defined independently of visit order.
            if minimize
                comp_residual = minimize_component!(
                    x, rs, rxns_idx, obs_set, baseline_vec, params, config, supply)
            end
            # Report the worst final residual across all loop components, so a
            # non-converged upstream loop isn't masked by a later converged one.
            last_change = max(last_change, comp_residual)
        end
    end
    stats = (method = scc_method, pooled = n_pooled, iterated = n_iterated,
             fallback_negative = n_fb_neg, fallback_inconsistent = n_fb_inc,
             pooled_nodes = n_pooled_nodes)
    return total_iters, last_change, stats
end

"""
Resolve one cyclic component as a conserved pool (specs/017). Returns
`(:pooled, n_written)`, `(:negative, 0)` (has an internal negative edge and
parity is off) or `(:inconsistent, 0)` (parity on, but a member is reached
from one entry with both signs -- an odd negative cycle). On the two fallback
statuses `x` is untouched and the caller iterates the component as before.

The external fold of a member reaction is its output with every non-pinned
member node at its own baseline, divided by the target's baseline. Reactions
whose target is pinned are not entries: the pinned value itself carries the
signal into the reactions that read it. Entry factors are multiplied in sorted
order so the product is bit-identical under relabelling.

Known consequences of the rule, measured (specs/017 research.md):
- a pooled state is not a fixed point of F, so `converged` reads false for
  any solve that pooled a component with a real entry;
- a pinned member does not sever the pool: entries on both sides of a pin
  multiply into one fold (the fixed point treats a pin as a boundary);
- alternative routes through one component (RAS isoforms in one GTPase
  cycle) are multiplied as if co-required, and 0 absorbs, so one isoform's
  knockout zeroes the cycle -- the mechanism behind the -53 held-out on MET,
  SCF-KIT and DAP12;
- DS_SCC_BREAK_CATALYST does not apply inside a pooled component.
"""
function pool_component!(x::Vector{Float64}, rs::Vector{Int}, rxns_idx::Vector{IndexedReaction},
                         comp_id::Vector{Int}, c::Int, members::Vector{Int}, obs_set::Set{Int},
                         baseline_vec::Vector{Float64}, config::ReactionEvalConfig,
                         parity::Bool, has_neg::Bool)
    (has_neg && !parity) && return :negative, 0
    xpool = copy(x)
    @inbounds for m in members
        m in obs_set || (xpool[m] = baseline_vec[m])
    end
    entries = Tuple{Int, Float64}[]          # (entry node, external fold)
    @inbounds for ri in rs
        r = rxns_idx[ri]
        t = r.target_idx
        t in obs_set && continue
        # closure edges (DS_SCC_BREAK_ROLES / legacy catalyst break) read the entry state
        f = compute_reaction_output_vec(xpool, r; supply=x, config=config) / baseline_vec[t]
        # hill_sat's smooth cap returns bl*(1 + 3e-15) at pure baseline, so an
        # exact `!= 1` test made ~97% of member reactions "entries" (TP53:
        # 808 of 837). A relative tolerance well below any real perturbation
        # and well above that noise keeps the notion meaningful.
        abs(f - 1.0) > 1e-9 && push!(entries, (t, f))
    end
    # Sign of every member relative to every entry (parity only).
    signs = Dict{Int, Dict{Int, Int8}}()     # entry node => (member => +-1)
    if parity && has_neg && !isempty(entries)
        adj = Dict{Int, Vector{Tuple{Int, Int8}}}()
        @inbounds for ri in rs
            r = rxns_idx[ri]
            t = r.target_idx
            for a in r.activator_indices
                comp_id[a] == c && push!(get!(adj, a, Tuple{Int, Int8}[]), (t, Int8(1)))
            end
            for i in r.inhibitor_indices
                comp_id[i] == c && push!(get!(adj, i, Tuple{Int, Int8}[]), (t, Int8(-1)))
            end
            for i in r.depletion_indices
                comp_id[i] == c && push!(get!(adj, i, Tuple{Int, Int8}[]), (t, Int8(-1)))
            end
        end
        for (e, _) in entries
            haskey(signs, e) && continue
            s = Dict{Int, Int8}(e => Int8(1))
            queue = [e]
            while !isempty(queue)
                u = popfirst!(queue)
                for (v, sg) in get(adj, u, Tuple{Int, Int8}[])
                    sv = s[u] * sg
                    if haskey(s, v)
                        s[v] == sv || return :inconsistent, 0
                    else
                        s[v] = sv
                        push!(queue, v)
                    end
                end
            end
            signs[e] = s
        end
    end
    n_written = 0
    factors = Float64[]
    @inbounds for m in members
        m in obs_set && continue
        empty!(factors)
        for (e, f) in entries
            sg = parity && has_neg ? get(signs[e], m, Int8(1)) : Int8(1)
            push!(factors, sg == 1 ? f : (f == 0.0 ? 100.0 : 1.0 / f))
        end
        sort!(factors)
        pool = 1.0
        for v in factors
            pool *= v
        end
        pool = clamp(pool, 0.0, 100.0)
        x[m] = clamp(baseline_vec[m] * pool, 0.0, 1.0)
        n_written += 1
    end
    return :pooled, n_written
end

"""
Steady-state solver: damped fixed-point iteration with observations pinned
as hard constraints. This is the right tool for the interactive perturbation
case ("user sets node X to value Y, propagate") because observations are
treated as ground truth and the network's fixed-point equilibrium relaxes
around them.

The spec's L_SS penalty formulation (soft observations, soft consistency)
is more appropriate for the *parameter-learning* mode where you're fitting
θ to many noisy multi-condition observations. We'll bring that back when we
implement the learning pass; for single-condition inference, the fixed-point
form gives biologically clean answers without the local-minima traps that
the penalty form has.

When DS_SCC_SOLVE is set, the iteration is replaced by solve_scc_ordered!
(SCC-condensation); otherwise the flat feed-forward (optionally DS_DAMPING)
loop below runs.
"""
function solve_steady_state_penalty(
    reactions::Vector{Reaction},
    observations::Dict{String, Tuple{Float64, Float64}},
    x0::Dict{String, Float64},
    baseline_activities::Dict{String, Float64},
    params::SteadyStateParams,
    start_time::Float64,
    gene_uuids::Set{String} = Set{String}();
    self_shared::Dict{Tuple{String, String}, Vector{String}} =
        Dict{Tuple{String, String}, Vector{String}}(),
)::SolverResult

    all_nodes = collect(keys(x0))
    n = length(all_nodes)
    uuid_to_idx = Dict(uuid => i for (i, uuid) in enumerate(all_nodes))
    index_stats = Dict{String, Any}()
    rxns_idx, comp_id, n_comp = index_reactions(reactions, uuid_to_idx, baseline_activities, gene_uuids;
                                               stats = index_stats, self_shared = self_shared)

    # Resolve the per-reaction DS_* knobs ONCE here (not once per reaction per
    # iteration inside compute_reaction_output_vec). Threaded into every forward
    # evaluation below.
    eval_config = resolve_reaction_eval_config()

    # Initial state vector.
    x = Vector{Float64}(undef, n)
    @inbounds for (i, uuid) in enumerate(all_nodes)
        x[i] = x0[uuid]
    end

    # Observation pinning: indices + values to enforce as hard constraints
    # after each iteration. Confidence is honored only at the >0 threshold —
    # any confident observation pins; confidence < tol observations are
    # ignored. (Future: weighted soft-pin when learning across conditions.)
    obs_indices = Int[]
    obs_values = Float64[]
    for (uuid, (val, conf)) in observations
        haskey(uuid_to_idx, uuid) || continue
        conf > OBS_CONFIDENCE_TOL || continue
        push!(obs_indices, uuid_to_idx[uuid])
        push!(obs_values, val / 100.0)  # UI 0-100 → internal 0-1
    end

    # Re-pin observed nodes to start.
    @inbounds for k in eachindex(obs_indices)
        x[obs_indices[k]] = obs_values[k]
    end

    # Feed-forward propagation with observations pinned. No damping between
    # reactions — each reaction's output is its deterministic F_r(inputs).
    # We iterate synchronously: compute F(x) for all reactions using the
    # current x, then assign. For acyclic networks this converges in
    # network-depth iterations; for cyclic networks (negative feedback loops)
    # it may oscillate, which we accept and report via the consistency residual.
    obs_set = Set(obs_indices)
    converged = false
    iters = 0
    scc_stats = (method = "flat", pooled = 0, iterated = 0, fallback_negative = 0,
                 fallback_inconsistent = 0, pooled_nodes = 0)
    # A closure activator with no `supply` reads the TARGET's baseline (the
    # legacy catalyst-break semantics); on the flat path, in the final residual
    # and in influence scores that zeroed the influence of every other input
    # under assembly-limiting. Pass the live state instead wherever closures
    # are marked: it IS the entry value there.
    flat_supply = _bool_env("DS_SCC_BREAK_CATALYST", false) || !isempty(_break_roles_env())
    max_change = 0.0

    # Optional damping for loop convergence. The original feed-forward
    # iteration (damping = 0) converges in network-depth steps for acyclic
    # networks but oscillates in cyclic ones (negative-feedback loops,
    # bistable switches). DS_DAMPING > 0 blends x ← (1-λ)·x + λ·F(x), which
    # is a contraction for bounded operators and converges to the same
    # fixed-point as the un-damped iteration when one exists. Observations
    # are still hard-pinned each step. Default 0 = original behaviour.
    damping = _float_env("DS_DAMPING", 0.0)
    if !(0.0 <= damping < 1.0)
        throw(ArgumentError(
            "DS_DAMPING=$damping is out of range; must be in [0, 1). " *
            "1.0 would freeze the iteration entirely."
        ))
    end

    scc_solve = _bool_env("DS_SCC_SOLVE", true)  # SCC-condensation solve is the default; set DS_SCC_SOLVE=0 for legacy flat iteration
    if scc_solve && n_comp > 0
        # SCC-condensation solve (see solve_scc_ordered!): solves the acyclic
        # majority exactly in topological order and confines damped iteration to
        # the strongly-connected components (real loops), so feedback converges
        # instead of settling on an oscillating last-iterate.
        # Baselines indexed like x, for the gamma||x - x0||^2 term of the
        # specified objective (DS_SCC_METHOD=minimize). Built here because
        # this is where the uuid -> index mapping lives.
        baseline_vec = Vector{Float64}(undef, n)
        @inbounds for (i, uuid) in enumerate(all_nodes)
            baseline_vec[i] = get(baseline_activities, uuid, 0.01)
        end
        iters, max_change, scc_stats = solve_scc_ordered!(
            x, rxns_idx, comp_id, n_comp, obs_set, params, eval_config, baseline_vec)
        converged = max_change < params.tolerance
    else
        for it in 1:params.max_iters
            x_fwd = forward_model_vec(x, rxns_idx; supply = flat_supply ? x : nothing, config=eval_config)

            max_change = 0.0
            @inbounds for i in 1:n
                if i in obs_set
                    continue  # pinned observation, don't update
                end
                new_val = damping > 0.0 ? (1.0 - damping) * x[i] + damping * x_fwd[i] : x_fwd[i]
                change = abs(new_val - x[i])
                change > max_change && (max_change = change)
                x[i] = new_val
            end

            iters = it
            if max_change < params.tolerance
                converged = true
                break
            end
        end
    end

    # Final consistency check on the stable state.
    x_fwd_final = forward_model_vec(x, rxns_idx; supply = flat_supply ? x : nothing, config=eval_config)
    consistency_inf = n == 0 ? 0.0 : maximum(abs.(x .- x_fwd_final))

    # Honest residual: max |x - F(x)| over the FREE nodes only.
    #
    # Two reasons not to report what was reported before:
    #  - `max_change` is written only inside the multi-node-SCC branch of
    #    solve_scc_ordered!. A network whose components are all singletons (the
    #    common case) therefore returned max_change = 0.0 and converged = true
    #    unconditionally, regardless of the actual residual.
    #  - the all-node `consistency_inf` is systematically large by construction,
    #    because a pinned observation is a hard constraint and is deliberately
    #    NOT at its model value; reporting it would call every perturbation
    #    solve non-converged.
    # Excluding pinned nodes measures the thing that actually matters: whether
    # the free variables reached a fixed point of the forward model.
    free_residual = 0.0
    saw_nonfinite = false
    @inbounds for i in 1:n
        i in obs_set && continue
        d = abs(x[i] - x_fwd_final[i])
        if !isfinite(d)
            # Track separately rather than skipping: filtering non-finite values
            # out of the max would leave free_residual finite and report
            # `converged = true` for a state containing NaN/Inf (whose activities
            # then serialize as JSON null). The code this replaced got that right
            # only incidentally, via `NaN < tolerance == false`.
            saw_nonfinite = true
        elseif d > free_residual
            free_residual = d
        end
    end

    converged = !saw_nonfinite && free_residual < params.tolerance
    safe_residual = saw_nonfinite ? Float64(consistency_inf) : Float64(free_residual)

    result_dict = Dict{String, Float64}(all_nodes[i] => clamp(x[i], 0.0, 1.0) for i in 1:n)
    solve_time = time() - start_time

    return SolverResult(
        result_dict,
        converged,
        iters,
        safe_residual,
        solve_time,
        Dict{String, Any}(
            "method" => "feed_forward_pinned",
            # Retained for continuity with earlier runs: the last per-iteration
            # delta from the loop components (0.0 when there were none).
            "max_change_final" => isfinite(max_change) ? Float64(max_change) : -1.0,
            "free_residual_inf" => safe_residual,
            "model_residual_inf" => isfinite(consistency_inf) ? Float64(consistency_inf) : -1.0,
            # specs/017: how the cyclic components were resolved, so an arm cannot
            # silently measure the old solver.
            "scc_method" => scc_stats.method,
            "scc_pooled" => scc_stats.pooled,
            "scc_iterated" => scc_stats.iterated,
            "scc_fallback_negative" => scc_stats.fallback_negative,
            "scc_fallback_inconsistent" => scc_stats.fallback_inconsistent,
            "scc_pooled_nodes" => scc_stats.pooled_nodes,
            # specs/018: recycling closures by role and the component census
            "scc_break_roles" => get(ENV, "DS_SCC_BREAK_ROLES", ""),
            "scc_closures_catalyst" => get(index_stats, "scc_closures_catalyst", 0),
            "scc_closures_assembly" => get(index_stats, "scc_closures_assembly", 0),
            "scc_closures_depletion" => get(index_stats, "scc_closures_depletion", 0),
            "scc_cyclic_before" => get(index_stats, "scc_cyclic_before", 0),
            # specs/022: inhibitor slots damped by DS_SELF_INHIBITOR_WEIGHT
            # (0 when off), so an arm cannot silently measure the old model.
            "self_inhibitors" => get(index_stats, "self_inhibitors", 0),
            "scc_cyclic_after" => get(index_stats, "scc_cyclic_after", 0),
            "scc_largest_after" => get(index_stats, "scc_largest_after", 0),
        ),
    )
end

"""
Compute influence scores for explainability.
Returns ranking of input nodes by their influence on the solution.
"""
function compute_influence_scores(
    result::SolverResult,
    reactions::Vector{Reaction};
    network::Union{Nothing, ReactionNetwork} = nothing,
)::Dict{String, Float64}

    # Influence = Σ_r |∂F_r/∂x_input| at the solved operating point, computed
    # with the SAME forward model the solver used (compute_reaction_output_vec,
    # which honors the active DS_* config — divide inhibition, hill_log AND,
    # depletion, assembly-limiting). The previous implementation used the
    # dict-form compute_reaction_output, which reads NONE of that config and
    # defaults inhibitor β=0, so it reported ZERO influence for every inhibitor
    # regardless of the config actually solved — misleading the moment any
    # explainability view surfaces it.
    all_nodes = collect(keys(result.node_activities))
    uuid_to_idx = Dict(uuid => i for (i, uuid) in enumerate(all_nodes))

    # Baselines are the universal spec x₀ = 0.01 (see tsv_parser). Passing an
    # empty dict lets index_reactions default every target to 0.01, matching
    # exactly what the solver built.
    # Pass `network` so the self-inhibitor rule (specs/022) is applied here too;
    # without it the scores describe the old double-counting model.
    self_shared = network === nothing ? Dict{Tuple{String, String}, Vector{String}}() :
                  first(self_inhibitor_setup(network))
    indexed, _, _ = index_reactions(reactions, uuid_to_idx, Dict{String, Float64}();
                                    self_shared = self_shared)

    x = Vector{Float64}(undef, length(all_nodes))
    @inbounds for (i, uuid) in enumerate(all_nodes)
        x[i] = result.node_activities[uuid]
    end

    eval_config = resolve_reaction_eval_config()
    infl_supply = _bool_env("DS_SCC_BREAK_CATALYST", false) || !isempty(_break_roles_env())
    influence_scores = Dict{String, Float64}()
    h = 1e-6
    for rxn in indexed
        base = compute_reaction_output_vec(x, rxn; supply = infl_supply ? x : nothing, config=eval_config)
        # Every input the reaction's output depends on: activators, inhibitors,
        # depletion (catalyst→substrate), and substrates.
        input_idxs = vcat(rxn.activator_indices, rxn.inhibitor_indices,
                          rxn.depletion_indices, rxn.substrate_indices)
        for idx in input_idxs
            x_saved = x[idx]
            # Step inward from the [0,1] boundary so a saturated input still
            # yields a finite-difference slope (forward step would clamp to 0).
            dir = (x_saved + h <= 1.0) ? h : -h
            x[idx] = x_saved + dir
            plus = compute_reaction_output_vec(x, rxn; supply = infl_supply ? x : nothing, config=eval_config)
            x[idx] = x_saved
            deriv = (plus - base) / dir
            uuid = all_nodes[idx]
            influence_scores[uuid] = get(influence_scores, uuid, 0.0) + abs(deriv)
        end
    end

    return influence_scores
end
