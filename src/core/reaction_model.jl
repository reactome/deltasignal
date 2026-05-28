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
"""
struct Reaction
    target_uuid::String
    activator_uuids::Vector{String}
    inhibitor_uuids::Vector{String}
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
    activators = String[]
    activator_is_and = Bool[]
    inhibitors = String[]
    inhibitor_is_and = Bool[]
    substrates = String[]
    products = String[]

    for edge in edges
        if edge.is_positive
            push!(activators, edge.parent_uuid)
            push!(activator_is_and, edge.is_and)
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
        substrates,
        products,
        params,
        is_and_gate,
        activator_is_and,
        inhibitor_is_and,
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
    inhibitor_indices::Vector{Int}
    substrate_indices::Vector{Int}
    params::ReactionParams
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
)::Vector{IndexedReaction}
    indexed = IndexedReaction[]
    for r in reactions
        haskey(uuid_to_idx, r.target_uuid) || continue

        act_indices = Int[]
        act_is_and = Bool[]
        for (k, uuid) in enumerate(r.activator_uuids)
            haskey(uuid_to_idx, uuid) || continue
            push!(act_indices, uuid_to_idx[uuid])
            # Defensive: fall back to is_and_gate if per-edge vector is
            # the wrong length for any reason.
            push!(act_is_and, k <= length(r.activator_is_and) ? r.activator_is_and[k] : r.is_and_gate)
        end

        push!(indexed, IndexedReaction(
            uuid_to_idx[r.target_uuid],
            get(baselines, r.target_uuid, 0.01),
            act_indices,
            act_is_and,
            [uuid_to_idx[u] for u in r.inhibitor_uuids if haskey(uuid_to_idx, u)],
            [uuid_to_idx[u] for u in r.substrate_uuids if haskey(uuid_to_idx, u)],
            r.params,
        ))
    end
    return indexed
end

"""
Compute reaction output from a flat activity vector. Same math as the dict
version (sensitivity transform → geomean activators × Hill-suppression
inhibitors × geomean substrates → Hill output), but generic-typed over the
element type so ForwardDiff Duals propagate.
"""
function compute_reaction_output_vec(x::AbstractVector{T}, rxn::IndexedReaction) where {T<:Real}
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
        transformed = apply_sensitivity_transform(
            x[i];
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
            else
                geometric_mean_aggregator(and_vals, and_wts)
            end
        end
        # OR cluster = biological alternatives ("any one of {A,B,C}").
        #   "max"  — reaction proceeds at the best-available alternative. Clean
        #            logical-OR reading, but masks single-member knockouts:
        #            max(0, x₀, x₀) = x₀ ⇒ no change registered.
        #   "mean" — each alternative carries weight, so knocking one out lowers
        #            the cluster. Matches the curator's assumption that a single
        #            ligand knockout has an effect even when paralogs exist.
        or_result = if isempty(or_vals)
            nothing
        else
            or_mode = get(ENV, "DS_OR_MODE", "max")
            or_mode == "mean" ? sum(or_vals) / length(or_vals) : maximum(or_vals)
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
        result = one(T)
        @inbounds for k in 1:length(rxn.inhibitor_indices)
            x_inh = clamp(x[rxn.inhibitor_indices[k]], zero(T), one(T))
            result *= (bl + eps_T) / (x_inh + eps_T)
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

    # Same rationale as the dict-version compute_reaction_output above:
    # drop the spec's output Hill step in favor of clamped linear propagation
    # to keep baseline as a stable fixed point. Hill primitives remain
    # available for explicit per-reaction use after training.
    return clamp(A * H * L, zero(T), one(T))
end

"""
Forward model F(x; θ) in vector form. Returns a new vector where each
target node carries the output of its driving reaction; nodes with no
incoming reaction keep their current value.
"""
function forward_model_vec(x::AbstractVector{T}, reactions::Vector{IndexedReaction}) where {T<:Real}
    y = copy(x)  # preserves element type, including Dual
    for rxn in reactions
        y[rxn.target_idx] = compute_reaction_output_vec(x, rxn)
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