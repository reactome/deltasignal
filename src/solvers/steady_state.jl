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
    return SteadyStateParams(
        1.0,    # mu
        0.1,    # gamma
        500,    # max_iters
        1e-6,   # tolerance
        "penalty"  # method
    )
end

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
    
    # Convert network to reactions
    reactions = convert_to_reaction_network(network)
    
    # Initialize node activities
    all_nodes = collect(keys(network.nodes))
    n_nodes = length(all_nodes)
    node_to_idx = Dict(uuid => i for (i, uuid) in enumerate(all_nodes))
    
    # Initial guess: use observations where available, baseline elsewhere
    x0 = Dict{String, Float64}()
    baseline_activities = Dict{String, Float64}()
    
    # Gene-entity UUIDs (nodes whose Reactome stable_id is a gene) — used to
    # detect transcriptional autoregulation loops. Empty unless DS_GENE_STIDS_FILE
    # is set, so this is a no-op by default.
    gene_stids = gene_stid_set()
    gene_uuids = Set{String}()
    for (uuid, node) in network.nodes
        if haskey(observations, uuid)
            x0[uuid] = observations[uuid][1] / 100.0  # Convert from UI scale 0-100 to internal 0-1
        else
            x0[uuid] = node.baseline
        end
        baseline_activities[uuid] = node.baseline
        if !isempty(gene_stids) && node.reactome_id !== nothing && node.reactome_id in gene_stids
            push!(gene_uuids, uuid)
        end
    end

    if params.method == "penalty"
        return solve_steady_state_penalty(reactions, observations, x0, baseline_activities, params, start_time, gene_uuids)
    else
        return solve_steady_state_fixed_point(reactions, observations, x0, baseline_activities, params, start_time)
    end
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
function solve_scc_ordered!(
    x::Vector{Float64},
    rxns_idx::Vector{IndexedReaction},
    comp_id::Vector{Int},
    n_comp::Int,
    obs_set::Set{Int},
    params::SteadyStateParams,
)
    λ = parse(Float64, get(ENV, "DS_SCC_DAMPING", "0.5"))
    use_supply = get(ENV, "DS_SCC_BREAK_CATALYST", "0") != "0"
    max_inner = params.max_iters
    tol = params.tolerance

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
    neg_frac_thresh = parse(Float64, get(ENV, "DS_SCC_NEG_FRAC", "0.5"))
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

    total_iters = 0
    last_change = 0.0

    # Topological order: upstream (high comp_id) before downstream (low).
    @inbounds for c in n_comp:-1:1
        rs = comp_rxns[c]
        isempty(rs) && continue

        if comp_size[c] == 1
            # Acyclic node: exact single evaluation (upstream already final).
            for ri in rs
                t = rxns_idx[ri].target_idx
                t in obs_set && continue
                x[t] = compute_reaction_output_vec(x, rxns_idx[ri])
            end
        else
            # Genuine loop: damped fixed point confined to this component.
            # Freeze the entry state for the optional catalyst-break layer.
            supply = use_supply ? copy(x) : nothing
            # Negative-feedback SCCs get a bounded (transient) relaxation; all
            # others relax to convergence.
            is_neg = neg_mode != "converge" && comp_neg_frac[c] >= neg_frac_thresh
            cap = is_neg ? neg_iters : max_inner
            for it in 1:cap
                total_iters += 1
                maxch = 0.0
                for ri in rs
                    r = rxns_idx[ri]
                    t = r.target_idx
                    t in obs_set && continue
                    fwd = compute_reaction_output_vec(x, r; supply=supply)
                    nv = (1.0 - λ) * x[t] + λ * fwd
                    ch = abs(nv - x[t])
                    ch > maxch && (maxch = ch)
                    x[t] = nv
                end
                last_change = maxch
                maxch < tol && break
            end
        end
    end
    return total_iters, last_change
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
    gene_uuids::Set{String} = Set{String}(),
)::SolverResult

    all_nodes = collect(keys(x0))
    n = length(all_nodes)
    uuid_to_idx = Dict(uuid => i for (i, uuid) in enumerate(all_nodes))
    rxns_idx, comp_id, n_comp = index_reactions(reactions, uuid_to_idx, baseline_activities, gene_uuids)

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
        conf > 1e-6 || continue
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
    max_change = 0.0

    # Optional damping for loop convergence. The original feed-forward
    # iteration (damping = 0) converges in network-depth steps for acyclic
    # networks but oscillates in cyclic ones (negative-feedback loops,
    # bistable switches). DS_DAMPING > 0 blends x ← (1-λ)·x + λ·F(x), which
    # is a contraction for bounded operators and converges to the same
    # fixed-point as the un-damped iteration when one exists. Observations
    # are still hard-pinned each step. Default 0 = original behaviour.
    damping = parse(Float64, get(ENV, "DS_DAMPING", "0.0"))

    scc_solve = get(ENV, "DS_SCC_SOLVE", "1") != "0"  # SCC-condensation solve is the default; set DS_SCC_SOLVE=0 for legacy flat iteration
    if scc_solve && n_comp > 0
        # SCC-condensation solve (see solve_scc_ordered!): solves the acyclic
        # majority exactly in topological order and confines damped iteration to
        # the strongly-connected components (real loops), so feedback converges
        # instead of settling on an oscillating last-iterate.
        iters, max_change = solve_scc_ordered!(
            x, rxns_idx, comp_id, n_comp, obs_set, params)
        converged = max_change < params.tolerance
    else
        for it in 1:params.max_iters
            x_fwd = forward_model_vec(x, rxns_idx)

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
    x_fwd_final = forward_model_vec(x, rxns_idx)
    consistency_inf = maximum(abs.(x .- x_fwd_final))

    # Sanitize residual to a finite Float64 — JSON can't serialize Inf/NaN.
    safe_residual = if isfinite(max_change)
        Float64(max_change)
    else
        Float64(consistency_inf)
    end

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
            "max_change_final" => safe_residual,
            "model_residual_inf" => isfinite(consistency_inf) ? Float64(consistency_inf) : -1.0,
        ),
    )
end

"""
Fixed-point iteration method.
Repeatedly applies x^{t+1} = λF(x^t) + (1-λ)x^t until convergence.
"""
function solve_steady_state_fixed_point(
    reactions::Vector{Reaction},
    observations::Dict{String, Tuple{Float64, Float64}},
    x0::Dict{String, Float64},
    baseline_activities::Dict{String, Float64},
    params::SteadyStateParams,
    start_time::Float64
)::SolverResult
    
    x_current = copy(x0)
    damping = 0.3  # Lower damping for better feedback loop dynamics
    
    for iter in 1:params.max_iters
        # Apply forward model
        x_forward = forward_model(x_current, reactions)
        
        # Update with damping
        x_new = Dict{String, Float64}()
        max_change = 0.0
        
        for uuid in keys(x_current)
            if haskey(x_forward, uuid)
                # Update: x^{t+1} = (1-λ)x^t + λF(x^t)
                x_new[uuid] = (1.0 - damping) * x_current[uuid] + damping * x_forward[uuid]
            else
                x_new[uuid] = x_current[uuid]
            end
            
            # Apply observations as constraints (soft)
            if haskey(observations, uuid)
                obs_value = observations[uuid][1] / 100.0  # Convert to internal scale
                confidence = observations[uuid][2]
                
                # Gentler blending to preserve feedback dynamics
                blend_factor = min(confidence * 0.3, 0.4)  # Much gentler blending
                x_new[uuid] = (1.0 - blend_factor) * x_new[uuid] + blend_factor * obs_value
            end
            
            # Clamp to valid range
            x_new[uuid] = clamp(x_new[uuid], 0.0, 1.0)
            
            # Track convergence
            change = abs(x_new[uuid] - x_current[uuid])
            max_change = max(max_change, change)
        end
        
        x_current = x_new
        
        # Check convergence
        if max_change < params.tolerance
            solve_time = time() - start_time
            
            return SolverResult(
                x_current,
                true,
                iter,
                max_change,
                solve_time,
                Dict("method" => "fixed_point", "damping" => damping)
            )
        end
        
        if iter % 50 == 0
            println("Fixed-point iteration $iter, max_change = $max_change")
        end
    end
    
    # Max iterations reached
    solve_time = time() - start_time
    
    return SolverResult(
        x_current,
        false,
        params.max_iters,
        NaN,
        solve_time,
        Dict("method" => "fixed_point", "warning" => "max_iterations_reached")
    )
end

"""
Compute influence scores for explainability.
Returns ranking of input nodes by their influence on the solution.
"""
function compute_influence_scores(
    result::SolverResult,
    reactions::Vector{Reaction}
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
    indexed, _, _ = index_reactions(reactions, uuid_to_idx, Dict{String, Float64}())

    x = Vector{Float64}(undef, length(all_nodes))
    @inbounds for (i, uuid) in enumerate(all_nodes)
        x[i] = result.node_activities[uuid]
    end

    influence_scores = Dict{String, Float64}()
    h = 1e-6
    for rxn in indexed
        base = compute_reaction_output_vec(x, rxn)
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
            plus = compute_reaction_output_vec(x, rxn)
            x[idx] = x_saved
            deriv = (plus - base) / dir
            uuid = all_nodes[idx]
            influence_scores[uuid] = get(influence_scores, uuid, 0.0) + abs(deriv)
        end
    end

    return influence_scores
end

"""
Generate upstream explanation: what minimal changes to inputs would achieve target outputs?
"""
function explain_upstream_drivers(
    network::ReactionNetwork,
    result::SolverResult,
    target_changes::Dict{String, Float64}  # Desired changes in target nodes
)::Dict{String, Float64}
    
    # This is a simplified version - a full implementation would use constrained optimization
    upstream_suggestions = Dict{String, Float64}()
    
    reactions = convert_to_reaction_network(network)
    influence_scores = compute_influence_scores(result, reactions)
    
    # Find nodes with high influence that are not targets
    target_nodes = Set(keys(target_changes))
    
    potential_drivers = [(uuid, score) for (uuid, score) in influence_scores 
                        if !(uuid in target_nodes) && score > 1e-6]
    
    # Sort by influence (descending)
    sort!(potential_drivers, by=x->x[2], rev=true)
    
    # Simple heuristic: suggest changes proportional to influence
    total_influence = sum(score for (_, score) in potential_drivers[1:min(5, length(potential_drivers))])
    
    for (i, (uuid, score)) in enumerate(potential_drivers[1:min(5, length(potential_drivers))])
        # Suggest change proportional to influence and desired target change
        avg_target_change = sum(values(target_changes)) / length(target_changes)
        suggested_change = (score / total_influence) * avg_target_change * 0.5  # Scale down for safety
        
        upstream_suggestions[uuid] = suggested_change
    end
    
    return upstream_suggestions
end