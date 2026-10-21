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
    return steady_state_params()
end

"""
    steady_state_params(; mu, gamma, max_iters, tolerance)

Solver parameters: an argument given here wins, else its `DS_*` variable, else
the code default; every value is validated the same way whichever it came from.
The CLI passes its flags here (code review 2026-10-02: it built
`SteadyStateParams` from its own defaults, so it ignored `DS_MAX_ITERS`,
`DS_TOLERANCE`, `DS_MU` and `DS_GAMMA`, which the API honours, and skipped
their validation).
"""
function steady_state_params(; mu = nothing, gamma = nothing, max_iters = nothing, tolerance = nothing)
    # DS_MAX_ITERS overrides the per-SCC iteration budget. Added because the
    # budget was hardcoded with no way to distinguish "this component is
    # oscillating" from "this component is converging but needs more sweeps" —
    # a distinction that decides whether the fixed-point iteration can work at
    # all here. Weaker damping made non-convergence WORSE (122 -> 152 TP53
    # cases at lambda=0.1), which points at slow convergence rather than
    # oscillation, but only a larger budget settles it.
    # An error names where the value came from: the argument (the CLI's flag)
    # or the DS_* variable.
    src(given, flag, var) = given ? flag : var
    n_it = src(max_iters !== nothing, "--max-iters", "DS_MAX_ITERS")
    n_mu = src(mu !== nothing, "--mu", "DS_MU")
    n_ga = src(gamma !== nothing, "--gamma", "DS_GAMMA")
    n_to = src(tolerance !== nothing, "--tolerance", "DS_TOLERANCE")
    max_iters = max_iters === nothing ? round(Int, _float_env("DS_MAX_ITERS", 500.0)) : max_iters
    if max_iters < 1
        throw(ArgumentError("$n_it=$max_iters must be at least 1."))
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
    mu = mu === nothing ? _float_env("DS_MU", 1.0) : Float64(mu)
    gamma = gamma === nothing ? _float_env("DS_GAMMA", 1e-6) : Float64(gamma)
    if !(mu > 0.0) || !isfinite(mu)
        throw(ArgumentError("$n_mu=$mu must be finite and > 0."))
    end
    if gamma < 0.0 || !isfinite(gamma)
        throw(ArgumentError("$n_ga=$gamma must be finite and >= 0."))
    end
    tolerance = tolerance === nothing ? _float_env("DS_TOLERANCE", 1e-6) : Float64(tolerance)
    if !(tolerance > 0.0) || !isfinite(tolerance)
        throw(ArgumentError("$n_to=$tolerance must be finite and > 0."))
    end
    return SteadyStateParams(mu, gamma, max_iters, tolerance, "penalty")
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
    # specs/081: recycle edges are structural only; no rule below reads them
    network, recycle_pairs = split_recycle_edges(network)

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
    # specs/048: pathogen-derived nodes, the same hold for an uninfected cell.
    pmode = pathogen_mode()
    pathogens = pmode == "inert" ? pathogen_uuids(network) : Set{String}()
    pathogen_rule = pmode == "propagate" ? "propagate" :
                    network.pathogen_stids === nothing ? "inert: no pathogen table" : "inert"
    # specs/049: a carrier loop's free form reads its supply, not its own release.
    # First, so every rule resolved on the network below sees the same edges.
    carrier_open = carrier_loop_open()   # specs/077
    network, carrier_loop_rule, carrier_loop_forms, carrier_loop_edges, carrier_loop_twins, carrier_open_edges =
        apply_carrier_loops(network, carrier_loop_mode(); open = carrier_open)
    # specs/059: a recycled copy of an entity also reads its root copies.
    network, recycled_twin_rule, recycled_twin_edges = apply_recycled_twins(network, recycled_twin_mode())
    # specs/063: a dead-end binding complex depletes its free partner.
    network, sequestration_rule, sequestration_edges = apply_sequestration(network, sequestration_mode())
    # specs/035: loop participants released unchanged by their loop, held at
    # baseline like cofactors. Resolved on the network as given (before bridges).
    cmode = conserved_mode()
    conserved = cmode == "inert" ? conserved_uuids(network) : Set{String}()

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
                                  network.containment, network.drug_stids, network.pools,
                                  network.pathogen_stids, network.carrier_loops, network.node_leaves)
    end

    # specs/052: under a pin to 0, nodes built from a protein no unpinned source
    # can supply are held at 0. Resolved on the network as solved (after carrier
    # loops and bridges), from the pins the caller gave.
    zmode = zero_supply_mode()
    zpinned = Dict{String, Float64}(u => v for (u, (v, c)) in observations if c > OBS_CONFIDENCE_TOL)
    zheld, zero_supply_rule, zero_supply_leaves = zero_supply_uuids(network, zpinned, zmode)

    # Convert network to reactions
    reactions = convert_to_reaction_network(network)

    # specs/022: inhibitors that contain their own reaction's input. Built unless
    # DS_SELF_INHIBITOR_WEIGHT=off. specs/040 rule B widens the pair test to
    # shared leaves; `self_leaf_pairs` counts the slots only that clause flags.
    self_shared, self_rule, self_leaf_pairs = self_inhibitor_setup(network)
    # specs/083: activators built from their own reaction's input
    self_act, self_act_rule = self_activator_setup(network)
    # specs/060 amendment 1: depleters whose entity contains their target's.
    own_entity = own_entity_depletion_pairs(network)
    # specs/061: entity-own pairs gated at solve time on their source copies.
    own_source = own_source_depletion_sources(network)
    # specs/040 rule A: resolved here so a misspelt mode fails before any solve.
    sfmode = self_fed_mode()
    # specs/066: resolved here so a misspelt mode fails before any solve.
    ncmode = nonconv_fallback_mode()
    # specs/073: reported as skip whenever the fallback itself is off.
    ncpools = ncmode == "lm" ? nonconv_pools_mode() : "skip"
    # specs/076: resolved here so a bad value fails before any solve.
    ncrestarts = ncmode == "lm" ? _nonconv_restarts() : 0
    _nonconv_restart_sweeps()
    
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
        # An observation that PINS (a measured drug, or a cofactor pin) wins; see
        # above. One the confidence gate would discard does not: skipping the
        # hold for it let a confidence-0 row silently unpin the drug (review).
        # (The cofactor block above still uses haskey; changing it would move
        # the default model, so it is left to its own change.)
        dpins = Dict{String, Tuple{Float64, Float64}}(
            u => (network.nodes[u].baseline * 100.0, 1.0)
            for u in drugs
            if !(haskey(observations, u) && observations[u][2] > OBS_CONFIDENCE_TOL))
        drugs_held = length(dpins)
        observations = merge(observations, dpins)
    end
    pathogens_held = 0
    if !isempty(pathogens)
        # The drug rule's pin precedence: an observation that pins wins.
        ppins = Dict{String, Tuple{Float64, Float64}}(
            u => (network.nodes[u].baseline * 100.0, 1.0)
            for u in pathogens
            if !(haskey(observations, u) && observations[u][2] > OBS_CONFIDENCE_TOL))
        pathogens_held = length(ppins)
        observations = merge(observations, ppins)
    end
    zero_supply_held = 0
    if !isempty(zheld)
        # Held at 0, not baseline: the leaf is absent. `zheld` excludes the
        # caller's pins, so an observation wins; a cofactor, drug or pathogen
        # hold merged above also wins, and is not counted here.
        zpins = Dict{String, Tuple{Float64, Float64}}(
            u => (0.0, 1.0) for u in zheld
            if haskey(network.nodes, u) &&
               !(haskey(observations, u) && observations[u][2] > OBS_CONFIDENCE_TOL))
        zero_supply_held = length(zpins)
        observations = merge(observations, zpins)
    end
    conserved_held = 0
    if !isempty(conserved)
        # An observation that PINS wins (a perturbed GAP such as NF1 stays
        # perturbed); one the confidence gate would discard does not.
        cpins = Dict{String, Tuple{Float64, Float64}}(
            u => (network.nodes[u].baseline * 100.0, 1.0)
            for u in conserved
            if haskey(network.nodes, u) &&
               !(haskey(observations, u) && observations[u][2] > OBS_CONFIDENCE_TOL))
        conserved_held = length(cpins)
        observations = merge(observations, cpins)
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
    # specs/039: interconversion pools solved at steady state (pi = pi P).
    ymode = cycle_mode()
    ypools = ymode == "balance" ? network.pools : nothing
    # Pools are solved only by the SCC fixed point; the other methods and the
    # legacy flat iteration never build them, so say so instead of "balance"
    # (code review 2026-10-02: they reported "balance" with 0 pools solved).
    smethod = get(ENV, "DS_SCC_METHOD", "fixed_point")
    y_inert = !_bool_env("DS_SCC_SOLVE", true) ? "DS_SCC_SOLVE=0" :
              smethod != "fixed_point" ? "DS_SCC_METHOD=$smethod" : nothing
    cycle_rule = ymode == "off" ? "off" : network.pools === nothing ? "balance: no pool table" :
                 y_inert !== nothing ? "balance: inert under $y_inert" : "balance"
    node_keys = Dict{String, String}(u => something(nd.reactome_id, "") for (u, nd) in network.nodes)
    supply_reach = ymode == "balance" ? pool_supply_reach() : "all"
    signal_pairs = supply_reach == "material" ? pool_signal_edge_pairs(network) : nothing
    crelease = ymode == "balance" && cycle_carriers() ? cycle_carrier_release() : "listed"
    # specs/075: resolved here so a misspelt mode fails before any solve.
    rcarriers = ymode == "balance" && cycle_carriers() ? cycle_reeval_carriers() : "off"
    penalty(os) = solve_steady_state_penalty(reactions, observations, x0, baseline_activities, params, start_time, gene_uuids;
                                             self_shared = self_shared,
                                             self_act = self_act,
                                             own_entity = own_entity,
                                             own_source = os,
                                             cycle_pools = ypools,
                                             cycle_phi = ymode == "balance" ? cycle_phi() : 0.1,
                                             cycle_carriers = ymode == "balance" ? cycle_carriers() : true,
                                             carrier_release = crelease,
                                             node_keys = node_keys,
                                             pool_signal_pairs = signal_pairs,
                                             structural = recycle_pairs)
    own_mode = get(ENV, "DS_DEPLETION_OWN_PRODUCT", DS_DEPLETION_OWN_PRODUCT_DEFAULT)
    if own_mode == "suppress_only_source"
        # specs/061 amendments 4-5: no cap is applied during an iteration. Pass 1
        # solves with every own-product depletion free (as `full`); an
        # entity-own edge that would de-repress (its depleter fell) and one of
        # whose source copies fell is removed; pass 2 solves without those
        # edges, again uncapped. Read live, the
        # `min(factor, 1)` cap and the source switch made TP53 and TGF-beta
        # oscillate (R5).
        pfull(rx) = solve_steady_state_penalty(rx, observations, x0, baseline_activities, params, start_time, gene_uuids;
                                               self_shared = self_shared, self_act = self_act, own_override = "full",
                                               cycle_pools = ypools,
                                               cycle_phi = ymode == "balance" ? cycle_phi() : 0.1,
                                               cycle_carriers = ymode == "balance" ? cycle_carriers() : true,
                                               carrier_release = crelease,
                                               node_keys = node_keys,
                                               pool_signal_pairs = signal_pairs,
                                               structural = recycle_pairs)
        first_pass = pfull(reactions)
        down(u) = haskey(first_pass.node_activities, u) &&
                  first_pass.node_activities[u] < get(baseline_activities, u, 0.01) * (1 - 1e-6)
        # amendment 5: node-copy own products are left free (as `full`); removing
        # or capping them made TP53 oscillate in every form tried.
        removed = Set{Tuple{String, String}}()
        for (k, srcs) in own_source
            down(k[1]) && any(down, srcs) && push!(removed, k)
        end
        if isempty(removed)
            result = first_pass
        else
            result = pfull(convert_to_reaction_network(without_depletions(network, removed)))
        end
        result.diagnostics["own_source_first_pass_converged"] = first_pass.converged
        result.diagnostics["own_source_depletions"] = length(removed)
        result.diagnostics["own_product_removed"] = sort!(collect(removed))
    else
        result = penalty(own_source)
    end
    delete!(result.diagnostics, "own_node_pairs")
    # Say which model ran: "on", "off", or inert because nothing told the solver
    # what contains what (a POSTed network or an older bundle).
    result.diagnostics["self_inhibitor_rule"] = self_rule
    result.diagnostics["self_activator_rule"] = self_act_rule   # specs/083
    result.diagnostics["cycle_rule"] = cycle_rule
    result.diagnostics["cycle_carriers"] = ymode == "balance" ? (cycle_carriers() ? "on" : "off") : "off"
    result.diagnostics["pool_supply_reach"] = supply_reach
    result.diagnostics["cycle_carrier_release"] = crelease
    result.diagnostics["cycle_reeval_carrier_rule"] = rcarriers
    result.diagnostics["drug_rule"] = drug_rule
    result.diagnostics["drugs_held"] = drugs_held
    result.diagnostics["pathogen_rule"] = pathogen_rule
    result.diagnostics["pathogens_held"] = pathogens_held
    result.diagnostics["carrier_loop_rule"] = carrier_loop_rule
    result.diagnostics["carrier_loop_forms"] = carrier_loop_forms
    result.diagnostics["carrier_loop_edges_dropped"] = carrier_loop_edges
    result.diagnostics["carrier_loop_twin_edges"] = carrier_loop_twins
    result.diagnostics["carrier_loop_open_rule"] = carrier_open
    result.diagnostics["carrier_loop_open_edges"] = carrier_open_edges
    result.diagnostics["recycled_twin_rule"] = recycled_twin_rule
    result.diagnostics["recycled_twin_edges"] = recycled_twin_edges
    result.diagnostics["sequestration_rule"] = sequestration_rule
    result.diagnostics["sequestration_edges"] = sequestration_edges
    result.diagnostics["own_entity_depletions"] = length(own_entity)
    # specs/061 (week review 2026-10-15): which own-product rule ran, now that
    # it is a default, so a run can show it measured the model it claims.
    result.diagnostics["depletion_own_product"] = get(ENV, "DS_DEPLETION_OWN_PRODUCT", DS_DEPLETION_OWN_PRODUCT_DEFAULT)
    result.diagnostics["zero_supply_rule"] = zero_supply_rule
    result.diagnostics["zero_supply_leaves"] = zero_supply_leaves
    result.diagnostics["zero_supply_held"] = zero_supply_held
    result.diagnostics["conserved_rule"] = cmode
    result.diagnostics["conserved_held"] = conserved_held
    # specs/040: which self-feedback rules ran. The counts come from the solve.
    result.diagnostics["self_fed_rule"] = sfmode
    result.diagnostics["nonconv_fallback"] = ncmode
    result.diagnostics["nonconv_pools"] = ncpools
    result.diagnostics["nonconv_restarts_max"] = ncrestarts
    result.diagnostics["self_inhibitor_leaves"] = self_inhibitor_leaves() ? "on" : "off"
    result.diagnostics["self_inhibitor_leaf_pairs"] = self_leaf_pairs
    return result
end

"""
(map, status) for DS_SELF_INHIBITOR_WEIGHT on this network. status is "off",
"on", or "inert: no containment table" -- the last means the default model is
NOT the one being run, which a caller must be able to see.
"""
function self_inhibitor_setup(network::ReactionNetwork)
    empty = Dict{Tuple{String, String}, Vector{String}}()
    leaves = self_inhibitor_leaves()        # validated even when the rule is off
    _self_inhibitor_weight_env() < 0 && return empty, "off", 0
    isempty(network.containment) && return empty, "inert: no containment table", 0
    strict = self_contained_inhibitor_map(network; leaves = false)
    leaves || return strict, "on", 0
    full = self_contained_inhibitor_map(network; leaves = true)
    # specs/040 rule B: slots the leaf clause flags that stId containment did not.
    n_leaf = count(k -> !haskey(strict, k), keys(full))
    return full, "on", n_leaf
end

"""
specs/083. `(map, rule)`: (activator uuid, target uuid) -> the target's `input`
activators that a positive `regulator` or `catalyst` of that target CONTAINS
(containment table) and is computed from (reachable from the input) -- a TF:gene
complex regulating that gene's transcription. The propagator counts the shared
input once (`DS_SELF_ACTIVATOR_WEIGHT`). Rule "off", "on" or "inert: no
containment table".
"""
function self_activator_setup(network::ReactionNetwork)
    out = Dict{Tuple{String, String}, Vector{String}}()
    _self_activator_weight_env() < 0 && return out, "off"
    isempty(network.containment) && return out, "inert: no containment table"
    sid(u) = (n = get(network.nodes, u, nothing); n === nothing ? nothing : n.reactome_id)
    inputs = Dict{String, Vector{String}}()
    fwd = Dict{String, Vector{String}}()
    for e in network.edges
        push!(get!(fwd, e.parent_uuid, String[]), e.child_uuid)
        e.is_positive && e.edge_type == "input" && push!(get!(inputs, e.child_uuid, String[]), e.parent_uuid)
    end
    reach_cache = Dict{String, Set{String}}()
    reaches(a) = get!(reach_cache, a) do
        seen = Set{String}(); stack = copy(get(fwd, a, String[]))
        while !isempty(stack)
            n = pop!(stack); n in seen && continue
            push!(seen, n); append!(stack, get(fwd, n, String[]))
        end
        seen
    end
    for e in network.edges
        (e.is_positive && e.edge_type in ("regulator", "catalyst")) || continue
        A = e.parent_uuid; s = sid(A)
        s === nothing && continue
        inside = get(network.containment, s, nothing)
        inside === nothing && continue
        shared = String[a for a in unique(get(inputs, e.child_uuid, String[]))
                        if a != A && (as = sid(a)) !== nothing && as != s && as in inside && A in reaches(a)]
        isempty(shared) || (out[(A, e.child_uuid)] = sort(shared))
    end
    return out, "on"
end

"""
specs/060 amendment 1. Under `DS_DEPLETION_OWN_PRODUCT=suppress_only_entity`,
the (depleter uuid, target uuid) pairs of depletion edges whose depleter's
entity CONTAINS the target's entity (the network's containment table): the
depleter is built from the target, so its absence cannot raise the target.
Empty under any other value, or with no containment table.
"""
function own_entity_depletion_pairs(network::ReactionNetwork; mode::String = "suppress_only_entity")::Set{Tuple{String, String}}
    out = Set{Tuple{String, String}}()
    get(ENV, "DS_DEPLETION_OWN_PRODUCT", DS_DEPLETION_OWN_PRODUCT_DEFAULT) == mode || return out
    isempty(network.containment) && return out
    sid(u) = (n = get(network.nodes, u, nothing); n === nothing ? nothing : n.reactome_id)
    for e in network.edges
        (e.is_positive || e.edge_type != "depletion") && continue
        d, t = sid(e.parent_uuid), sid(e.child_uuid)
        (d === nothing || t === nothing) && continue
        t in get(network.containment, d, Set{String}()) && push!(out, (e.parent_uuid, e.child_uuid))
    end
    return out
end

"""
specs/061. Under `DS_DEPLETION_OWN_PRODUCT=suppress_only_source` (default), each
entity-own depletion pair (as `own_entity_depletion_pairs`) -> its SOURCE
copies: nodes with the target's stable id reached walking backward from the
depleter along material edges (not depletion, catalyst or regulator;
amendment 2), stopping at each such copy (the copies the depleter is built
from). The target itself is never a source
(amendment 1: that test was circular); the walk passes through it. The propagator caps the
edge at 1 while any source reads below baseline. Pairs with no source are
left out (not capped). Empty under any other value.
"""
const OWN_SOURCE_SKIP_EDGES = ("depletion", "catalyst", "regulator")

"""
specs/071. The (source, target) pairs whose every positive edge is a catalyst
or regulator edge: routes along which no material moves. Read by the pool
supply test under `DS_POOL_SUPPLY_REACH=material`. A pair that is also an
input (an entity that is both catalyst and substrate) carries material and
is not listed.
"""
function pool_signal_edge_pairs(network::ReactionNetwork)::Set{Tuple{String, String}}
    signal = Set{Tuple{String, String}}(); material = Set{Tuple{String, String}}()
    for e in network.edges
        e.is_positive || continue
        e.edge_type == "depletion" && continue
        push!(e.edge_type in ("catalyst", "regulator") ? signal : material, (e.parent_uuid, e.child_uuid))
    end
    return setdiff(signal, material)
end

"""`network` without the depletion edges (depleter, target) in `pairs` (specs/061 amendment 4)."""
function without_depletions(network::ReactionNetwork, pairs::Set{Tuple{String, String}})
    edges = [e for e in network.edges
             if !(e.edge_type == "depletion" && (e.parent_uuid, e.child_uuid) in pairs)]
    return ReactionNetwork(network.nodes, edges, network.set_mappings, network.cofactor_stids,
                           network.containment, network.drug_stids, network.pools, network.pathogen_stids,
                           network.carrier_loops, network.node_leaves)
end

function own_source_depletion_sources(network::ReactionNetwork)::Dict{Tuple{String, String}, Vector{String}}
    out = Dict{Tuple{String, String}, Vector{String}}()
    pairs = own_entity_depletion_pairs(network; mode = "suppress_only_source")
    isempty(pairs) && return out
    back = Dict{String, Vector{String}}()
    for e in network.edges
        # amendment 2: only edges that carry material; a copy that merely
        # catalyses or regulates the depleter's reaction is not built into it.
        e.edge_type in OWN_SOURCE_SKIP_EDGES && continue
        push!(get!(back, e.child_uuid, String[]), e.parent_uuid)
    end
    for (d, t) in pairs
        tid = network.nodes[t].reactome_id
        seen = Set{String}([d])
        stack = copy(get(back, d, String[]))
        src = String[]
        while !isempty(stack)
            n = pop!(stack)
            n in seen && continue
            push!(seen, n)
            nd = get(network.nodes, n, nothing)
            # amendment 1: the target is never its own source (its value
            # includes the depletion being decided); walk through it.
            if n != t && nd !== nothing && nd.reactome_id == tid
                push!(src, n)
                continue
            end
            append!(stack, get(back, n, String[]))
        end
        isempty(src) || (out[(d, t)] = sort!(src))
    end
    return out
end

"""
(inhibitor uuid, target uuid) -> the target's activator uuids that the
inhibitor CONTAINS, according to the network's containment table (specs/022).
A node's identity is its own `reactome_id`; activators are collected per target
node, so an inhibitor is never matched against a sibling variant's input.
Depletion edges are not inhibitors here and are never included.
"""
function self_contained_inhibitor_map(network::ReactionNetwork;
                                      leaves::Bool = false)::Dict{Tuple{String, String}, Vector{String}}
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
    # specs/040 rule B: containment at the leaf level. A leaf is a contained
    # stId that contains nothing else; cofactor leaves (ATP, GTP ...) are not
    # shared species. The inhibitor and the input are "built from the same
    # species" when a carrier of a shared leaf reaches both.
    cont = network.containment
    leaf_cache = Dict{String, Set{String}}()
    function leaves_of(st)
        st === nothing && return Set{String}()
        get!(leaf_cache, st) do
            cands = union(get(cont, st, Set{String}()), Set([st]))
            Set{String}(c for c in cands
                        if isempty(setdiff(get(cont, c, Set{String}()), Set([c]))) &&
                           !(c in network.cofactor_stids))
        end
    end
    leaf_shared(a_st, i_st) = !isempty(intersect(leaves_of(a_st), leaves_of(i_st)))
    carriers = Dict{String, Vector{String}}()
    if leaves
        for u in sort(collect(keys(network.nodes)))
            st = network.nodes[u].reactome_id
            st === nothing && continue
            for l in leaves_of(st); push!(get!(carriers, l, String[]), u); end
        end
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
        I = e.parent_uuid
        function common_cause(a, a_st)
            for l in intersect(leaves_of(a_st), leaves_of(s))
                for n in get(carriers, l, String[])
                    n == I && continue
                    (n == a || a in reaches(n)) && I in reaches(n) && return true
                end
            end
            return false
        end
        shared = String[a for a in unique(get(acts, e.child_uuid, String[]))
                        if (as = sid(a)) !== nothing &&
                           ((as in inside && I in reaches(a)) ||
                            (leaves && leaf_shared(as, s) && (I in reaches(a) || common_cause(a, as))))]
        isempty(shared) || (out[(I, e.child_uuid)] = sort(shared))
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
function lm_minimize!(z::Vector{Float64}, p;
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
    supply::Union{Nothing,Vector{Float64}};
    lm_cap::Union{Nothing,Int} = nothing,
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
            cap = lm_cap === nothing ? max(1, min(params.max_iters, _lm_iter_cap())) : lm_cap
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
Every node the pool update `update_pools!(x, pools, ...)` reads (specs/073):
the pools' own nodes and carriers, the inputs of their step copies, of the
re-evaluated drives, of the states' and carriers' own updates (their
modifiers), and the states' supply producers.
"""
function pool_update_inputs(pools, rxns_idx)::Set{Int}
    out = Set{Int}()
    for p in pools
        union!(out, p.forms, p.flux_nodes, p.copy_nodes, p.carriers, p.reeval)
        rxs = Int[]
        for st in p.path_steps, (ris, src) in st
            append!(rxs, ris); push!(out, src)
        end
        append!(rxs, p.reeval_rxn, p.carrier_rxn)
        append!(rxs, filter(!=(0), p.form_rxn))
        for prods in p.supply_rxns; union!(out, prods); end
        for ri in rxs
            r = rxns_idx[ri]
            union!(out, r.activator_indices, r.inhibitor_indices, r.depletion_indices, r.substrate_indices)
        end
    end
    return out
end

"""
A component holding pool-managed nodes, as a least-squares problem over all
its free nodes (specs/073). A managed node's row is `G(x) - x`, `G` the pool
update; any other node's row is `F(x) - x`, as in `ComponentProblem`. No prior.
"""
struct PoolComponentProblem
    base::ComponentProblem
    pools::Vector{IndexedPool}
    shared::Set{Int}
    managed::Set{Int}
    obs_set::Set{Int}
    carriers::Bool
    gcols::Vector{Int}          # free nodes G reads, as node indices
    gbuf::Vector{Float64}
    gamma::Float64
end

function PoolComponentProblem(base::ComponentProblem, pools, shared, managed, obs_set, rxns_idx, carriers)
    gcols = sort!([j for j in pool_update_inputs(pools, rxns_idx) if haskey(base.slot, j)])
    PoolComponentProblem(base, pools, shared, managed, obs_set, carriers, gcols, copy(base.xbuf), 0.0)
end

# G(xbuf) into gbuf: a copy of xbuf with the pool update applied.
function _pool_map!(p::PoolComponentProblem)
    copyto!(p.gbuf, p.base.xbuf)
    update_pools!(p.gbuf, p.pools, p.shared, p.base.rxns, p.base.baseline, p.obs_set, p.base.config;
                  carriers = p.carriers)
    return p.gbuf
end

function component_residual_jacobian(p::PoolComponentProblem, z::Vector{Float64}; gamma::Float64 = 0.0)
    b = p.base
    @inbounds for (k, t) in enumerate(b.free); b.xbuf[t] = z[k]; end
    rows = Int[]            # rs positions with a free target, in order
    @inbounds for (k, ri) in enumerate(b.rs)
        haskey(b.slot, b.rxns[ri].target_idx) && push!(rows, k)
    end
    m = length(rows)
    r = zeros(Float64, m)
    I = Int[]; Jc = Int[]; V = Float64[]
    g0 = copy(_pool_map!(p))
    grow = Dict{Int, Int}()     # managed target -> row
    @inbounds for (row, k) in enumerate(rows)
        rx = b.rxns[b.rs[k]]
        t = rx.target_idx
        if t in p.managed
            r[row] = g0[t] - b.xbuf[t]
            grow[t] = row
        else
            d = b.deps[k]
            r[row] = compute_reaction_output_vec(b.xbuf, rx; supply=b.supply, config=b.config) - b.xbuf[t]
            if !isempty(d)
                local_in = Float64[b.xbuf[i] for i in d]
                fr = function (v)
                    yy = OverlayVector{eltype(v), Float64}(b.xbuf, d, v)
                    compute_reaction_output_vec(yy, rx; supply=b.supply, config=b.config)
                end
                gloc = ForwardDiff.gradient(fr, local_in)
                for (mm, i) in enumerate(d)
                    push!(I, row); push!(Jc, b.slot[i]); push!(V, gloc[mm])
                end
            end
        end
        push!(I, row); push!(Jc, b.slot[t]); push!(V, -1.0)
    end
    # G's columns by forward differences (backward at the upper bound).
    if !isempty(grow)
        @inbounds for j in p.gcols
            xj = b.xbuf[j]
            h = 1e-8 + 1e-6 * abs(xj)
            xj + h > 1.0 && (h = -h)
            b.xbuf[j] = xj + h
            g1 = _pool_map!(p)
            b.xbuf[j] = xj
            for (t, row) in grow
                dg = (g1[t] - g0[t]) / h
                dg == 0.0 && continue
                push!(I, row); push!(Jc, b.slot[j]); push!(V, dg)
            end
        end
    end
    return r, sparse(I, Jc, V, m, length(b.free))
end

"""Worst residual of a pool component: `|F(x) - x|` over its free unmanaged
targets and `|G(x) - x|` over the nodes its pools write (specs/073), the two
quantities the final convergence check reads."""
function _component_residual_pools(x, rs, rxns_idx, obs_set, supply, config, managed, pools, shared,
                                   baseline_vec, carriers)::Float64
    resid = 0.0
    @inbounds for ri in rs
        r = rxns_idx[ri]
        t = r.target_idx
        (t in obs_set || t in managed) && continue
        d = abs(compute_reaction_output_vec(x, r; supply=supply, config=config) - x[t])
        d > resid && (resid = d)
    end
    d = update_pools!(copy(x), pools, shared, rxns_idx, baseline_vec, obs_set, config; carriers = carriers)
    return isfinite(d) ? max(resid, d) : Inf
end

"""
Levenberg-Marquardt over a pool component's free nodes (specs/073), from `x`.
Writes the result into `x` and returns `_component_residual_pools`.
"""
function minimize_component_pools!(x, rs, rxns_idx, obs_set, baseline_vec, config, supply,
                                   pools, shared, managed, carriers; lm_cap::Int, tol::Float64)::Float64
    base = ComponentProblem(x, rs, rxns_idx, obs_set, baseline_vec, 1.0, 0.0, config, supply)
    if !isempty(base.free)
        prob = PoolComponentProblem(base, pools, shared, managed, obs_set, rxns_idx, carriers)
        z = Float64[x[t] for t in base.free]
        try
            lm_minimize!(z, prob; max_iter=lm_cap, g_tol=tol * 1e-2, gamma=0.0)
            @inbounds for (k, t) in enumerate(base.free)
                x[t] = clamp(z[k], 0.0, 1.0)
            end
        catch err
            @warn "pool component minimisation failed; keeping the iterate" exception=(err, catch_backtrace()) n_free=length(base.free)
        end
    end
    return _component_residual_pools(x, rs, rxns_idx, obs_set, supply, config, managed, pools, shared,
                                     baseline_vec, carriers)
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

"""
`DS_NONCONV_FALLBACK` (specs/066): `lm` (default since 2026-10-16) | `off`. Under `lm`, a cyclic
component whose damped fixed point ends at or above tolerance is re-solved by
Levenberg-Marquardt on the pure residual |F(x) - x|^2 (no baseline prior),
started from the last iterate, for `DS_NONCONV_LM_ITERS` (default 400)
iterations. The result is kept only if its residual is lower. Converged
components are never touched. Components holding cycle-balance pool nodes are
skipped here: those nodes are solved by their pool (pi = piP), not by F.
`DS_NONCONV_POOLS=joint` (default, specs/073) re-solves them with the pool
update in the residual.

Why: some components are Goodwin rings (a chain of lags closed by a high-gain
negative feedback, e.g. a complex depleting each of its own free subunits
across an AND). Their fixed point exists but is unstable under any damped
relaxation, so the iteration settles into a limit cycle and the scored value
is the phase the iteration cap lands on.
"""
const DS_VALID_NONCONV_FALLBACK = ("off", "lm")

function nonconv_fallback_mode()::String
    v = get(ENV, "DS_NONCONV_FALLBACK", "lm")
    v in DS_VALID_NONCONV_FALLBACK || throw(ArgumentError(
        "DS_NONCONV_FALLBACK=\"$v\" is not a recognised mode; expected \"off\" or \"lm\"."))
    return v
end

function _nonconv_lm_iters()::Int
    n = round(Int, _float_env("DS_NONCONV_LM_ITERS", 400.0))
    n < 1 && throw(ArgumentError("DS_NONCONV_LM_ITERS=$n must be at least 1."))
    n
end

"""
`DS_NONCONV_POOLS` (specs/073): `joint` (default since 2026-10-18) | `skip`. What the specs/066
fallback does with a non-converged component that holds cycle-balance pool
nodes. Under `skip` it is left at the last iterate and counted
(`nonconv_skipped_pools`). Under `joint` it is re-solved by
Levenberg-Marquardt over all its free nodes at once: a node a pool writes
contributes the row `G(x) - x`, where `G` is the pool update (`update_pools!`),
and every other node the row `F(x) - x`. That is the residual the final
convergence check reads. `G`'s Jacobian is taken by forward differences over
the inputs the pool update reads (`pool_update_inputs`). The result is kept
only if that residual is lower.
"""
const DS_VALID_NONCONV_POOLS = ("skip", "joint")

"""
`DS_NONCONV_RESTARTS` (specs/076): how many times a component the specs/066
fallback leaves above tolerance is iterated `DS_NONCONV_RESTART_SWEEPS`
(default 250) more sweeps from its last iterate and the fallback retried.
Default 4 since 2026-10-19 (specs/076, adopted by Adam; 0 restores the
previous solver). LM's result depends on where it starts: in TP53 on
the faithful build it stalls from the cap-500 iterate of MDM4 overexpression
and converges, to one state, from the iterates at 750 to 2,000. The best state
over all attempts is kept.
"""
function _nonconv_restarts()::Int
    n = round(Int, _float_env("DS_NONCONV_RESTARTS", 4.0))
    n < 0 && throw(ArgumentError("DS_NONCONV_RESTARTS=$n must be at least 0."))
    n
end

function _nonconv_restart_sweeps()::Int
    n = round(Int, _float_env("DS_NONCONV_RESTART_SWEEPS", 250.0))
    n < 1 && throw(ArgumentError("DS_NONCONV_RESTART_SWEEPS=$n must be at least 1."))
    n
end

function nonconv_pools_mode()::String
    v = get(ENV, "DS_NONCONV_POOLS", "joint")
    v in DS_VALID_NONCONV_POOLS || throw(ArgumentError(
        "DS_NONCONV_POOLS=\"$v\" is not a recognised mode; expected \"skip\" or \"joint\"."))
    return v
end

"""
`DS_SCC_ORDER`: the order reactions are swept in inside a cyclic component
under Gauss-Seidel (specs/013 amendment 2).

- `label` (historical; the default until 2026-10-03): the order of the indexed reactions, which comes
  from `Dict` iteration over node uuids. A rebuild re-mints uuids, so it
  re-orders the sweep: canonical and an identical rebuild differ by 158 TP53
  predictions (specs/013 amendment 1).
- `structure`: sorted by a label-free colour of each reaction's target node --
  its stable id refined by its neighbours' (Weisfeiler-Lehman, see
  `structural_node_colours`). Only structurally identical positions tie; a tie
  keeps the label order and is counted (`scc_order_ties`). Label-free, but it
  converged 1.3-7% less often than the label order (amendment 2).
- `flow` (**default**, amendment 3: rebuilds score identically, convergence
  and held-out unchanged): downstream first -- by breadth-first distance inside
  the component from its entries (members fed from outside it, or pinned),
  then by the structural colour. Gauss-Seidel converges fastest sweeping with
  the signal.
"""
const SCC_ORDERS = ("label", "structure", "flow")

function scc_order_mode()::String
    m = get(ENV, "DS_SCC_ORDER", "flow")
    m in SCC_ORDERS || throw(ArgumentError(
        "DS_SCC_ORDER=$m is not an order; expected one of $(join(SCC_ORDERS, ", "))."))
    return m
end

"""
    structural_node_colours(rxns_idx, node_key; rounds = 4) -> Vector{UInt64}

A colour per node that depends only on the network's structure and its stable
ids, never on uuids: start from `hash(stable id)`, then `rounds` times rehash
each node with the sorted colours of its in- and out-neighbours, each tagged by
the edge's role (AND / OR activator, inhibitor, depletion, substrate). Roles
are integer codes: `hash(::Symbol)` is by object identity, which differs
between processes. `hash` of strings, integers, tuples and vectors of them is
deterministic for a given Julia version.
"""
function structural_node_colours(rxns_idx::Vector{IndexedReaction}, node_key::Vector{String};
                                 rounds::Int = 4)::Vector{UInt64}
    n = length(node_key)
    col = UInt64[hash(k) for k in node_key]
    ins = [Tuple{Int, Int}[] for _ in 1:n]
    outs = [Tuple{Int, Int}[] for _ in 1:n]
    link(role, s, t) = (push!(ins[t], (role, s)); push!(outs[s], (role, t)))
    for r in rxns_idx
        t = r.target_idx
        for (j, s) in enumerate(r.activator_indices)
            link(r.activator_is_and[j] ? 1 : 2, s, t)
        end
        for s in r.inhibitor_indices; link(3, s, t); end
        for s in r.depletion_indices; link(4, s, t); end
        for s in r.substrate_indices; link(5, s, t); end
    end
    for _ in 1:rounds
        new = similar(col)
        @inbounds for i in 1:n
            a = sort!(UInt64[hash((role, col[s])) for (role, s) in ins[i]])
            b = sort!(UInt64[hash((role, col[t])) for (role, t) in outs[i]])
            new[i] = hash((col[i], a, b))
        end
        col = new
    end
    return col
end

function solve_scc_ordered!(
    x::Vector{Float64},
    rxns_idx::Vector{IndexedReaction},
    comp_id::Vector{Int},
    n_comp::Int,
    obs_set::Set{Int},
    params::SteadyStateParams,
    config::ReactionEvalConfig = resolve_reaction_eval_config(),
    baseline_vec::Union{Nothing,Vector{Float64}} = nothing;
    pools::Vector{IndexedPool} = IndexedPool[],
    pool_stats::Dict{String, Any} = Dict{String, Any}(),
    pool_carriers::Bool = true,
    pool_signal_pairs::Union{Nothing, Set{Tuple{Int, Int}}} = nothing,
    entry_supply::Union{Nothing, Vector{Float64}} = nothing,
    order_key::Union{Nothing, Vector{UInt64}} = nothing,
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
    use_supply = _bool_env("DS_SCC_BREAK_CATALYST", false) || (!isempty(_break_roles_env()) || self_fed_mode() == "entry")
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
    nonconv_lm = nonconv_fallback_mode() == "lm"        # specs/066
    nonconv_joint = nonconv_lm && nonconv_pools_mode() == "joint"   # specs/073
    nonconv_restarts = nonconv_lm ? _nonconv_restarts() : 0          # specs/076
    nonconv_restart_sweeps = _nonconv_restart_sweeps()
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
    # specs/013 amendment 2: a label-free sweep order. MergeSort is stable, so a
    # tie keeps the label order; ties are counted so a run says how much order
    # is still left to labels.
    if order_key !== nothing
        flow = get(pool_stats, "scc_order", "structure") == "flow"
        n_tied = 0
        for c in 1:n_comp
            rs_c = comp_rxns[c]
            length(rs_c) > 1 || continue
            dist = Dict{Int, Int}()
            if flow
                # entries: a target fed from outside the component, or pinned
                succ = Dict{Int, Vector{Int}}()
                queue = Int[]
                for ri in rs_c
                    r = rxns_idx[ri]; t = r.target_idx
                    srcs = vcat(r.activator_indices, r.inhibitor_indices,
                                r.depletion_indices, r.substrate_indices)
                    if t in obs_set || any(s -> comp_id[s] != c, srcs)
                        haskey(dist, t) || (dist[t] = 0; push!(queue, t))
                    end
                    for s_ in srcs
                        comp_id[s_] == c && push!(get!(succ, s_, Int[]), t)
                    end
                end
                head = 1
                while head <= length(queue)
                    u = queue[head]; head += 1
                    for v in get(succ, u, Int[])
                        haskey(dist, v) || (dist[v] = dist[u] + 1; push!(queue, v))
                    end
                end
            end
            key(ri) = (get(dist, rxns_idx[ri].target_idx, typemax(Int)), order_key[rxns_idx[ri].target_idx])
            sort!(rs_c; alg = MergeSort, by = key)
            for k in 2:length(rs_c)
                key(rs_c[k - 1]) == key(rs_c[k]) && (n_tied += 1)
            end
        end
        pool_stats["scc_order_ties"] = n_tied
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
    self_fed_on = self_fed_mode() == "entry"
    # Node lists per component: for pooling, for rule A, and for the entry
    # snapshot the final residual reads (specs/040).
    need_nodes = pooling || self_fed_on || entry_supply !== nothing
    comp_nodes = need_nodes ? [Int[] for _ in 1:n_comp] : Vector{Int}[]
    comp_has_neg = falses(n_comp)
    if need_nodes
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

    # specs/039: a pool is solved in closed form when all its nodes (forms and
    # transition reactions) sit in ONE cyclic component; its nodes then skip the
    # ordinary update. Any other pool is left to the iteration and counted, as is
    # every pool under a component method that does not run the sweep below
    # (pool*, which returns before it; minimize, which overwrites it).
    comp_pools = Dict{Int, Vector{IndexedPool}}()
    managed = Set{Int}()
    n_pool_unmanaged = 0
    # Mass-flow out-edges (activator source -> target) of every update, for the
    # supply test below. Depletion and inhibitor edges do not carry the pool's
    # protein: following them marked RAS's own maturation supply as fed by the
    # pool (RAS:GTP:BRAP -| mature RAS), so KRAS 80x read NORMAL (review of
    # amendment 2, finding 1).
    pool_out = isempty(pools) ? Dict{Int, Vector{Int}}() : begin
        d = Dict{Int, Vector{Int}}()
        for r in rxns_idx, s in r.activator_indices
            push!(get!(d, s, Int[]), r.target_idx)
        end
        d
    end
    # Every node some pool writes (states, intermediates, step copies): a carrier
    # must not be one of them, or two updates write it (finding 2).
    pool_written = Set{Int}()
    for p in pools
        union!(pool_written, p.forms, p.flux_nodes, p.copy_nodes)
    end
    n_carrier_conflicts = 0; n_pinned_fallback = 0; n_supply_depl = 0; n_held_inputs = 0
    n_supply_kept = 0
    rxn_of = isempty(pools) ? Dict{Int, Int}() : Dict{Int, Int}(r.target_idx => ri for (ri, r) in enumerate(rxns_idx))
    n_carriers_held = 0; n_multistep = 0
    carriers_seen = Set{Int}()   # a carrier shared by k pools is counted once
    for p in pools
        nodes = vcat(p.forms, p.flux_nodes, p.copy_nodes)
        c0 = comp_id[nodes[1]]
        # A pinned intermediate or step copy has no place in the chain (its value
        # is not the flux through it): the pool falls back to the iteration and
        # is counted, as pre-registered (finding 5). A pinned STATE is handled:
        # it sets the supply.
        if any(v -> v in obs_set, p.flux_nodes) || any(v -> v in obs_set, p.copy_nodes)
            n_pool_unmanaged += 1; n_pinned_fallback += 1
            continue
        end
        if baseline_vec !== nothing && !pooling && !minimize &&
           all(v -> comp_id[v] == c0, nodes) && comp_size[c0] > 1
            # Supply excludes a producer the pool's own forms reach, since it
            # can be fed by the pool itself (review of specs/039, finding 6).
            # Reachability does not pass through pinned nodes (cofactors,
            # drugs, observations): a pinned node cannot carry the pool's
            # feedback, and a component welded only through GDP must not hide
            # the protein's supply (Fable review, finding 2).
            reach = Set{Int}(p.forms); stack = copy(p.forms)
            while !isempty(stack)
                v = pop!(stack)
                for w in get(pool_out, v, Int[])
                    (w in reach || w in obs_set || comp_id[w] != c0) && continue
                    push!(reach, w); push!(stack, w)
                end
            end
            # specs/071: under DS_POOL_SUPPLY_REACH=material the supply test
            # follows only material edges (pool_signal_pairs are the
            # catalyst/regulator-only ones), so a producer the pool merely
            # regulates stays a supply.
            sreach = reach
            if pool_signal_pairs !== nothing
                sreach = Set{Int}(p.forms); stack = copy(p.forms)
                while !isempty(stack)
                    v = pop!(stack)
                    for w in get(pool_out, v, Int[])
                        (w in sreach || w in obs_set || comp_id[w] != c0 || (v, w) in pool_signal_pairs) && continue
                        push!(sreach, w); push!(stack, w)
                    end
                end
            end
            for prods in p.supply_rxns
                filter!(a -> !(a in sreach), prods)
                n_supply_kept += count(a -> a in reach && !(a in sreach), prods)
            end
            # The supply is exogenous to the pool. A DERIVED depletion edge from a
            # node the pool feeds onto the supply chain (the activator ancestry of
            # the supply producers inside the component) is the pool's protein
            # consuming its own precursor, which the pool's conservation already
            # holds; it is not applied, as the catalyst -| source-state edge is
            # not (amendment 3, review finding: RAS:GTP:BRAP -| mature RAS made
            # the RAS pool suppress its own supply, SOS1 80x 4.38 not 8.99).
            # Curated inhibitor edges are untouched.
            anc = Set{Int}(); stk = Int[]
            for prods in p.supply_rxns, a in prods
                a in anc || (push!(anc, a); push!(stk, a))
            end
            while !isempty(stk)
                v = pop!(stk)
                ri = get(rxn_of, v, 0)
                ri == 0 && continue
                for a in rxns_idx[ri].activator_indices
                    (a in anc || a in reach || a in obs_set || comp_id[a] != c0) && continue
                    pool_signal_pairs !== nothing && (a, v) in pool_signal_pairs && continue
                    push!(anc, a); push!(stk, a)
                end
            end
            for v in sort!(collect(anc))
                ri = get(rxn_of, v, 0)
                ri == 0 && continue
                r = rxns_idx[ri]
                drop = findall(i -> i in reach, r.depletion_indices)
                isempty(drop) && continue
                n_supply_depl += length(drop)
                length(r.depletion_break) == length(r.depletion_indices) && deleteat!(r.depletion_break, drop)
                length(r.depletion_own_product) == length(r.depletion_indices) && deleteat!(r.depletion_own_product, drop)
                length(r.depletion_own_src) == length(r.depletion_indices) && deleteat!(r.depletion_own_src, drop)
                deleteat!(r.depletion_indices, drop)
            end
            # amendment 5 (post hoc): pool-fed inputs of the pool's own steps,
            # and the pool-fed nodes they depend on, are re-evaluated in the
            # drives with the pool at baseline. Order: BFS from the pool over
            # sorted ids (fixed, label-independent up to the id order).
            poolset = Set(vcat(p.forms, p.flux_nodes, p.copy_nodes))
            allcar = Set(p.carriers)
            fed(v) = v in reach && !(v in poolset) && !(v in allcar)
            need = Set{Int}(); stk = Int[]
            for st in p.path_steps, (ris, _) in st, ri in ris
                r = rxns_idx[ri]
                for i in Iterators.flatten((r.activator_indices, r.inhibitor_indices, r.depletion_indices))
                    fed(i) && !(i in need) && (push!(need, i); push!(stk, i); n_held_inputs += 1)
                end
            end
            while !isempty(stk)
                v = pop!(stk)
                ri = get(rxn_of, v, 0)
                ri == 0 && continue
                r = rxns_idx[ri]
                for i in Iterators.flatten((r.activator_indices, r.inhibitor_indices, r.depletion_indices))
                    fed(i) && !(i in need) && (push!(need, i); push!(stk, i))
                end
            end
            # BFS order from the pool's states through activator edges
            order = Int[]; seen = Set{Int}(p.forms); q = sort(copy(p.forms))
            while !isempty(q)
                v = popfirst!(q)
                v in need && push!(order, v)
                for w in sort(get(pool_out, v, Int[]))
                    (w in seen || !(w in reach)) && continue
                    push!(seen, w); push!(q, w)
                end
            end
            for v in order
                ri = get(rxn_of, v, 0)
                ri == 0 && continue
                push!(p.reeval, v); push!(p.reeval_rxn, ri)
            end
            # carriers are managed only inside the pool's component, and never
            # a node the pool already writes
            keepc = [comp_id[c] == c0 && !(c in pool_written) for c in p.carriers]
            n_carrier_conflicts += count(c -> c in pool_written, p.carriers)
            deleteat!(p.carrier_rxn, findall(!, keepc)); deleteat!(p.carrier_release, findall(!, keepc))
            deleteat!(p.carriers, findall(!, keepc))
            push!(get!(comp_pools, c0, IndexedPool[]), p)
            union!(managed, nodes)
            if pool_carriers
                union!(managed, p.carriers)
                for c in p.carriers
                    (c in obs_set || c in carriers_seen) && continue
                    push!(carriers_seen, c); n_carriers_held += 1
                end
            end
            any(st -> length(st) > 1, p.path_steps) && (n_multistep += 1)
        else
            n_pool_unmanaged += 1
        end
    end
    # specs/049: nodes two pools of one component both write are written once,
    # from both
    comp_shared = Dict(c => shared_pool_nodes(ps) for (c, ps) in comp_pools)
    # specs/075: a pool-fed input the drives re-evaluate that is the carrier of a
    # pool in the same component reads by the carrier rule under `carrier`.
    # Counted under both modes, so `full` reports what `carrier` would change.
    reeval_carrier_mode = isempty(comp_pools) ? "full" : cycle_reeval_carriers()
    n_reeval_carriers = 0
    for ps in values(comp_pools)
        crel_c = Dict{Int, Set{Int}}()
        if pool_carriers
            for q in ps, (c, rel) in zip(q.carriers, q.carrier_release)
                union!(get!(crel_c, c, Set{Int}()), rel)
            end
        end
        for p in ps
            empty!(p.reeval_excl)
            for v in p.reeval
                ex = get(crel_c, v, nothing)
                ex === nothing || (n_reeval_carriers += 1)
                push!(p.reeval_excl, reeval_carrier_mode == "carrier" ? ex : nothing)
            end
        end
    end
    pool_stats["cycle_reeval_carriers"] = n_reeval_carriers
    pool_stats["cycle_shared_nodes"] = sum(length(v) for v in values(comp_shared); init = 0)
    pool_stats["cycle_pools_solved"] = sum(length(v) for v in values(comp_pools); init = 0)
    pool_stats["cycle_pools_unmanaged"] = n_pool_unmanaged
    pool_stats["cycle_pools_multistep"] = n_multistep
    pool_stats["cycle_carriers_held"] = n_carriers_held
    pool_stats["cycle_carriers_flag"] = pool_carriers
    pool_stats["cycle_carrier_conflicts"] = n_carrier_conflicts
    pool_stats["cycle_pools_pinned_fallback"] = n_pinned_fallback
    pool_stats["cycle_supply_depletions_held"] = n_supply_depl
    pool_stats["cycle_supply_regulated_kept"] = n_supply_kept
    pool_stats["cycle_self_fed_inputs"] = n_held_inputs
    # For the final consistency check: managed nodes are held by their pool, not
    # by the forward model (like pinned nodes), and are checked against it.
    pool_stats["cycle_managed"] = managed
    pool_stats["cycle_solved"] = reduce(vcat, collect(values(comp_pools)); init = IndexedPool[])

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
            # specs/040 rule A: self-fed inputs of this component read the
            # entry value. Computed ONCE here, from the entry state; the flags
            # live on this solve's IndexedReactions only. See core/self_fed.jl.
            if self_fed_on && baseline_vec !== nothing
                nodes_c = comp_nodes[c]
                in_c = Set(nodes_c)
                act_out = Dict{Int, Vector{Int}}()
                is_pool = Set{Int}()
                entry = Set{Int}()
                for ri in rs
                    r = rxns_idx[ri]; t = r.target_idx
                    is_set_pool_reaction(r) && push!(is_pool, t)
                    for s in r.activator_indices
                        s in in_c && push!(get!(act_out, s, Int[]), t)
                    end
                    if t in obs_set
                        abs(x[t] - baseline_vec[t]) > SELF_FED_TOL && push!(entry, t)
                    elseif any(s -> !(s in in_c) && abs(x[s] - baseline_vec[s]) > SELF_FED_TOL,
                               r.activator_indices)
                        push!(entry, t)
                    end
                end
                # specs/039 owns pool states, intermediates and step copies:
                # never candidates; an entry iff already off baseline.
                skip = copy(is_pool)
                for p in get(comp_pools, c, IndexedPool[])
                    union!(skip, p.forms); union!(skip, p.flux_nodes); union!(skip, p.copy_nodes)
                end
                for u in nodes_c
                    if (u in obs_set || u in skip) && abs(x[u] - baseline_vec[u]) > SELF_FED_TOL
                        push!(entry, u)
                    end
                end
                blocked = Set{Int}(u for u in nodes_c if u in obs_set && !(u in entry))
                sf = self_fed_nodes(nodes_c, act_out, is_pool, skip, entry, blocked)
                n_held = 0
                for ri in rs
                    r = rxns_idx[ri]
                    for k in eachindex(r.activator_indices)
                        s = r.activator_indices[k]
                        if s in sf && s != r.target_idx && k <= length(r.activator_break)
                            r.activator_break[k] || (n_held += 1)
                            r.activator_break[k] = true
                        end
                    end
                end
                pool_stats["self_fed_nodes"] = get(pool_stats, "self_fed_nodes", 0) + length(sf)
                pool_stats["self_fed_edges_held"] = get(pool_stats, "self_fed_edges_held", 0) + n_held
            end
            # Genuine loop: damped fixed point confined to this component.
            # Freeze the entry state for the optional catalyst-break layer.
            supply = use_supply ? copy(x) : nothing
            # The final residual must read held edges at the SAME value the
            # iteration did, or it reads the loop live and reports a residual
            # of up to 1 on a converged solve (specs/040 trace: 3 of 8 RAF
            # cases). Record the entry state of this component's nodes.
            if entry_supply !== nothing
                for u in comp_nodes[c]; entry_supply[u] = x[u]; end
            end
            # Negative-feedback SCCs get a bounded (transient) relaxation; all
            # others relax to convergence.
            is_neg = neg_mode != "converge" && comp_neg_frac[c] >= neg_frac_thresh
            cap = is_neg ? neg_iters : max_inner
            comp_residual = 0.0
            if jacobi && length(rs) > length(stage_idx)
                resize!(stage_idx, length(rs))
                resize!(stage_val, length(rs))
            end
            # specs/076: under DS_NONCONV_RESTARTS=k, a component the fallback
            # leaves above tolerance is iterated `DS_NONCONV_RESTART_SWEEPS` more
            # sweeps from the last ITERATE (not the LM state) and the fallback
            # retried, up to k times. LM's answer depends on where it starts; the
            # best state over all attempts is kept.
            budget = cap
            restarts_left = nonconv_restarts
            best_res = Inf; best_x = Float64[]; iterate_x = Float64[]
            attempted = false
            while true
                for it in 1:budget
                    total_iters += 1
                    maxch = 0.0
                    n_staged = 0
                    for ri in rs
                        r = rxns_idx[ri]
                        t = r.target_idx
                        t in obs_set && continue
                        t in managed && continue       # specs/039: solved by its pool below
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
                    if haskey(comp_pools, c)
                        ch = update_pools!(x, comp_pools[c], comp_shared[c], rxns_idx, baseline_vec, obs_set, config;
                                           carriers = pool_carriers)
                        ch > maxch && (maxch = ch)
                    end
                    comp_residual = maxch
                    maxch < tol && break
                end
                # specs/066: a component still above tolerance after the cap is
                # usually a limit cycle (a Goodwin ring), whose fixed point exists
                # but repels damped relaxation. Find it by least squares on the
                # residual instead, and keep it only if it is closer to F(x) = x.
                # A component bounded on purpose (DS_SCC_NEG_MODE) is left alone.
                if nonconv_lm && !minimize && !is_neg && comp_residual >= tol
                    # specs/076: the specs/066 and 073 counters describe the first
                    # attempt; restarts are counted separately
                    first = restarts_left == nonconv_restarts
                    first && (pool_stats["nonconv_components"] = get(pool_stats, "nonconv_components", 0) + 1)
                    iterate_x = x[comp_nodes[c]]
                    attempted = true
                    if haskey(comp_pools, c) && nonconv_joint
                        # specs/073: pool nodes join the least-squares problem as
                        # rows G(x) - x
                        first && (pool_stats["nonconv_pool_components"] = get(pool_stats, "nonconv_pool_components", 0) + 1)
                        saved = x[comp_nodes[c]]
                        r0 = _component_residual_pools(x, rs, rxns_idx, obs_set, supply, config, managed,
                                                       comp_pools[c], comp_shared[c], baseline_vec, pool_carriers)
                        r1 = minimize_component_pools!(x, rs, rxns_idx, obs_set, baseline_vec, config, supply,
                                                       comp_pools[c], comp_shared[c], managed, pool_carriers;
                                                       lm_cap = _nonconv_lm_iters(), tol = tol)
                        if r1 < r0
                            comp_residual = r1
                            first && (pool_stats["nonconv_pool_improved"] = get(pool_stats, "nonconv_pool_improved", 0) + 1)
                            first && r1 < tol && (pool_stats["nonconv_pool_rescued"] = get(pool_stats, "nonconv_pool_rescued", 0) + 1)
                        else
                            x[comp_nodes[c]] = saved
                        end
                    elseif haskey(comp_pools, c) || any(u -> u in managed, comp_nodes[c])
                        pool_stats["nonconv_skipped_pools"] = get(pool_stats, "nonconv_skipped_pools", 0) + 1
                        attempted = false
                    else
                        saved = x[comp_nodes[c]]
                        r0 = _component_residual(x, rs, rxns_idx, obs_set, supply, config)
                        pz = SteadyStateParams(1.0, 0.0, params.max_iters, tol, params.method)
                        r1 = minimize_component!(x, rs, rxns_idx, obs_set, baseline_vec, pz, config, supply;
                                                 lm_cap = _nonconv_lm_iters())
                        if r1 < r0
                            comp_residual = r1
                            first && (pool_stats["nonconv_improved"] = get(pool_stats, "nonconv_improved", 0) + 1)
                            first && r1 < tol && (pool_stats["nonconv_rescued"] = get(pool_stats, "nonconv_rescued", 0) + 1)
                        else
                            x[comp_nodes[c]] = saved
                        end
                    end
                end
                # specs/076: keep the best state; restart from the iterate.
                if attempted && comp_residual >= tol && nonconv_restarts > 0
                    if comp_residual < best_res
                        best_res = comp_residual; best_x = x[comp_nodes[c]]
                    end
                    if restarts_left > 0
                        restarts_left -= 1
                        pool_stats["nonconv_restarts"] = get(pool_stats, "nonconv_restarts", 0) + 1
                        x[comp_nodes[c]] = iterate_x
                        budget = nonconv_restart_sweeps
                        continue
                    end
                    x[comp_nodes[c]] = best_x; comp_residual = best_res
                elseif nonconv_restarts > 0 && restarts_left < nonconv_restarts && comp_residual < tol
                    pool_stats["nonconv_restart_rescued"] = get(pool_stats, "nonconv_restart_rescued", 0) + 1
                end
                break
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
    self_act::Dict{Tuple{String, String}, Vector{String}} = Dict{Tuple{String, String}, Vector{String}}(),
    self_shared::Dict{Tuple{String, String}, Vector{String}} =
        Dict{Tuple{String, String}, Vector{String}}(),
    own_entity::Set{Tuple{String, String}} = Set{Tuple{String, String}}(),
    own_source::Dict{Tuple{String, String}, Vector{String}} =
        Dict{Tuple{String, String}, Vector{String}}(),
    own_override::Union{Nothing, String} = nothing,
    cycle_pools = nothing,
    cycle_phi::Float64 = 0.1,
    cycle_carriers::Bool = true,
    carrier_release::String = "listed",
    node_keys::Union{Nothing, Dict{String, String}} = nothing,
    pool_signal_pairs::Union{Nothing, Set{Tuple{String, String}}} = nothing,
    structural::Vector{Tuple{String, String}} = Tuple{String, String}[],
)::SolverResult

    all_nodes = collect(keys(x0))
    n = length(all_nodes)
    uuid_to_idx = Dict(uuid => i for (i, uuid) in enumerate(all_nodes))
    index_stats = Dict{String, Any}()
    rxns_idx, comp_id, n_comp = index_reactions(reactions, uuid_to_idx, baseline_activities, gene_uuids;
                                               stats = index_stats, self_shared = self_shared, self_act = self_act,
                                               own_entity = own_entity, own_source = own_source,
                                               structural = structural)
    index_stats["recycle_edges"] = length(structural)
    order_mode = scc_order_mode()
    index_stats["scc_order"] = order_mode
    order_key = nothing
    if order_mode in ("structure", "flow")
        node_keys === nothing && throw(ArgumentError(
            "DS_SCC_ORDER=$order_mode needs the network's stable ids; this solve was given none."))
        order_key = structural_node_colours(rxns_idx, String[get(node_keys, u, "") for u in all_nodes])
    end

    # Resolve the per-reaction DS_* knobs ONCE here (not once per reaction per
    # iteration inside compute_reaction_output_vec). Threaded into every forward
    # evaluation below.
    eval_config = resolve_reaction_eval_config()
    # specs/061 amendment 4: the two-pass own-product rule solves uncapped
    own_override === nothing || (eval_config = with_own_product(eval_config, own_override))

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
    flat_supply = _bool_env("DS_SCC_BREAK_CATALYST", false) || (!isempty(_break_roles_env()) || self_fed_mode() == "entry")
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
        ipools, n_skipped = cycle_pools === nothing ? (IndexedPool[], 0) :
                            index_pools(cycle_pools, uuid_to_idx, rxns_idx, cycle_phi;
                                        release_keys = carrier_release == "reaction" && node_keys !== nothing ?
                                            Dict{Int, String}(uuid_to_idx[u] => k for (u, k) in node_keys
                                                              if !isempty(k) && haskey(uuid_to_idx, u)) : nothing,
                                        stats = index_stats)
        index_stats["cycle_pools_skipped"] = n_skipped
        # Component-entry values, for the final residual below: a held edge
        # (specs/018 closure, specs/040 self-fed) reads these, not the live state.
        # Only iterated components record one (NaN elsewhere); every other node
        # is checked against its live value, as before (review of specs/040).
        entry_supply = fill(NaN, length(x))
        iters, max_change, scc_stats = solve_scc_ordered!(
            x, rxns_idx, comp_id, n_comp, obs_set, params, eval_config, baseline_vec;
            pools = ipools, pool_stats = index_stats, pool_carriers = cycle_carriers,
            pool_signal_pairs = pool_signal_pairs === nothing ? nothing :
                Set{Tuple{Int, Int}}((uuid_to_idx[a], uuid_to_idx[b]) for (a, b) in pool_signal_pairs
                                     if haskey(uuid_to_idx, a) && haskey(uuid_to_idx, b)),
            entry_supply = entry_supply, order_key = order_key)
        @inbounds for i in eachindex(entry_supply)
            isnan(entry_supply[i]) && (entry_supply[i] = x[i])
        end
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
    x_fwd_final = forward_model_vec(x, rxns_idx;
                                    supply = flat_supply ? (scc_solve && n_comp > 0 ? entry_supply : x) : nothing,
                                    config=eval_config)
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
    # specs/039: a pool-managed node is at its pool's steady state, not at its
    # forward-model value; its residual is how far one more pool update moves it.
    managed_final = get(index_stats, "cycle_managed", Set{Int}())
    solved_final = get(index_stats, "cycle_solved", IndexedPool[])
    if !isempty(solved_final)
        xp = copy(x)
        bvec = Vector{Float64}(undef, n)
        @inbounds for (i, uuid) in enumerate(all_nodes)
            bvec[i] = get(baseline_activities, uuid, 0.01)
        end
        d = update_pools!(xp, solved_final, shared_pool_nodes(solved_final), rxns_idx, bvec, obs_set,
                          eval_config; carriers = cycle_carriers)
        isfinite(d) || (saw_nonfinite = true)
        isfinite(d) && d > free_residual && (free_residual = d)
    end
    @inbounds for i in 1:n
        i in obs_set && continue
        i in managed_final && continue
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
            "self_activators" => get(index_stats, "self_activators", 0),   # specs/083
            # specs/061: entity-own depletion edges gated on their source copies
            "own_source_depletions" => get(index_stats, "own_source_depletions", 0),
            # specs/061 amendment 4: node-copy own-product pairs (uuid), read by
            # the two-pass rule and removed from the result before it returns
            "own_node_pairs" => get(index_stats, "own_node_pairs", Tuple{String, String}[]),
            "scc_cyclic_after" => get(index_stats, "scc_cyclic_after", 0),
            "scc_largest_after" => get(index_stats, "scc_largest_after", 0),
            # specs/039: interconversion pools solved in closed form, left to the
            # iteration (not within one cyclic component), or skipped at indexing.
            "cycle_pools_solved" => get(index_stats, "cycle_pools_solved", 0),
            "recycle_edges" => get(index_stats, "recycle_edges", 0),   # specs/081
            "cycle_pools_unmanaged" => get(index_stats, "cycle_pools_unmanaged", 0),
            "cycle_pools_skipped" => get(index_stats, "cycle_pools_skipped", 0),
            "cycle_pools_multistep" => get(index_stats, "cycle_pools_multistep", 0),
            # specs/049: intermediates and step copies written from several pools
            "cycle_shared_nodes" => get(index_stats, "cycle_shared_nodes", 0),
            "cycle_carriers_held" => get(index_stats, "cycle_carriers_held", 0),
            "cycle_carrier_conflicts" => get(index_stats, "cycle_carrier_conflicts", 0),
            "cycle_pools_pinned_fallback" => get(index_stats, "cycle_pools_pinned_fallback", 0),
            "cycle_supply_depletions_held" => get(index_stats, "cycle_supply_depletions_held", 0),
            "cycle_supply_regulated_kept" => get(index_stats, "cycle_supply_regulated_kept", 0),
            "cycle_carrier_release_added" => get(index_stats, "cycle_carrier_release_added", 0),
            # specs/075: re-evaluated pool-fed inputs that are carriers (pairs).
            "cycle_reeval_carriers" => get(index_stats, "cycle_reeval_carriers", 0),
            "cycle_self_fed_inputs" => get(index_stats, "cycle_self_fed_inputs", 0),
            # specs/040 rule A: nodes found self-fed and edges held at entry (0 when off).
            "self_fed_nodes" => get(index_stats, "self_fed_nodes", 0),
            # specs/013 amendment 2: the sweep order and how many adjacent
            # reactions in a cyclic component it could not tell apart.
            "scc_order" => get(index_stats, "scc_order", "label"),
            "scc_order_ties" => get(index_stats, "scc_order_ties", 0),
            "self_fed_edges_held" => get(index_stats, "self_fed_edges_held", 0),
            # specs/066: non-converged components, and what the fallback did.
            "nonconv_components" => get(index_stats, "nonconv_components", 0),
            "nonconv_improved" => get(index_stats, "nonconv_improved", 0),
            "nonconv_rescued" => get(index_stats, "nonconv_rescued", 0),
            "nonconv_skipped_pools" => get(index_stats, "nonconv_skipped_pools", 0),
            # specs/073: pool components re-solved jointly (0 under skip).
            "nonconv_pool_components" => get(index_stats, "nonconv_pool_components", 0),
            # specs/076: restarts used, and components brought below tolerance on one.
            "nonconv_restarts" => get(index_stats, "nonconv_restarts", 0),
            "nonconv_restart_rescued" => get(index_stats, "nonconv_restart_rescued", 0),
            "nonconv_pool_improved" => get(index_stats, "nonconv_pool_improved", 0),
            "nonconv_pool_rescued" => get(index_stats, "nonconv_pool_rescued", 0),
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

    # Index as solve_steady_state does (silo bridges aside: an off-by-default
    # arm whose bridged edges this pass does not see): the network's own baselines
    # (an empty dict defaulted every one to 0.01, so a network with baselines
    # 0.5 read influence 0 where the derivative is 1 -- code review
    # 2026-10-02), its gene uuids, and the self-inhibitor rule (specs/022).
    # Without `network` the scores describe the default-baseline model.
    #
    # Known limit: a node set by the specs/039 pool solve, or a supply edge the
    # cyclic-component solve held, is differentiated through its reaction's
    # forward model here, which is not how the solve set it. Carrier loops
    # (specs/049) are not a limit: callers pass `solved_carrier_network`.
    self_shared = network === nothing ? Dict{Tuple{String, String}, Vector{String}}() :
                  first(self_inhibitor_setup(network))
    self_act = network === nothing ? Dict{Tuple{String, String}, Vector{String}}() :
               first(self_activator_setup(network))   # specs/083
    baselines = Dict{String, Float64}()
    gene_uuids = Set{String}()
    if network !== nothing
        gene_stids = gene_stid_set()
        for (u, node) in network.nodes
            baselines[u] = node.baseline
            if !isempty(gene_stids) && node.reactome_id !== nothing && node.reactome_id in gene_stids
                push!(gene_uuids, u)
            end
        end
    end
    # specs/061 amendment 4: a two-pass solve ran uncapped on the network without
    # the edges it removed; differentiate that model, not the live cap.
    removed = get(result.diagnostics, "own_product_removed", nothing)
    two_pass = removed !== nothing && network !== nothing
    if two_pass
        reactions = convert_to_reaction_network(without_depletions(network, Set{Tuple{String, String}}(removed)))
    end
    indexed, _, _ = index_reactions(reactions, uuid_to_idx, baselines, gene_uuids;
                                    self_shared = self_shared, self_act = self_act,
                                    own_entity = network === nothing || two_pass ? Set{Tuple{String, String}}() :
                                                 own_entity_depletion_pairs(network),
                                    own_source = network === nothing || two_pass ? Dict{Tuple{String, String}, Vector{String}}() :
                                                 own_source_depletion_sources(network))

    x = Vector{Float64}(undef, length(all_nodes))
    @inbounds for (i, uuid) in enumerate(all_nodes)
        x[i] = result.node_activities[uuid]
    end

    eval_config = resolve_reaction_eval_config()
    two_pass && (eval_config = with_own_product(eval_config, "full"))
    infl_supply = _bool_env("DS_SCC_BREAK_CATALYST", false) || (!isempty(_break_roles_env()) || self_fed_mode() == "entry")
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
