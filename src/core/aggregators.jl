# Multi-Input Aggregation Functions

"""
Geometric mean aggregator for activators with weights.
A = exp(Σ wᵢ * log(x̃ᵢ + ε))

This is equivalent to: A = ∏(x̃ᵢ + ε)^wᵢ
where Σ wᵢ = 1 (weights form a simplex)

Parameters:
- inputs: vector of transformed activator signals x̃ᵢ ∈ [0,1]
- weights: vector of weights wᵢ ≥ 0, Σ wᵢ = 1
- ε: small constant for numerical stability
"""
function geometric_mean_aggregator(
    inputs::AbstractVector{T},
    weights::AbstractVector{<:Real};
    ε::Real=1e-12,
) where {T<:Real}

    if length(inputs) != length(weights)
        error("Number of inputs must match number of weights")
    end

    if length(inputs) == 0
        return one(T)  # Identity for no inputs
    end

    w_sum = sum(weights)
    if w_sum <= 0
        error("Sum of weights must be positive")
    end

    # exp(Σ wᵢ · log(x̃ᵢ + ε)). Use a generic-typed accumulator so autodiff
    # Dual types propagate through.
    log_sum = zero(T)
    for i in eachindex(inputs)
        x_safe = clamp(inputs[i], zero(T), one(T)) + ε
        log_sum += (weights[i] / w_sum) * log(x_safe)
    end

    return clamp(exp(log_sum), zero(T), one(T))
end

"""
Product aggregator for AND logic.
A = ∏ x̃ᵢ
"""
function product_aggregator(inputs::Vector{Float64})::Float64
    if length(inputs) == 0
        return 1.0
    end
    
    result = 1.0
    for x in inputs
        result *= clamp(x, 0.0, 1.0)
    end
    
    return result
end

"""
Maximum aggregator for OR logic.
A = max(x̃ᵢ)
"""
function maximum_aggregator(inputs::Vector{Float64})::Float64
    if length(inputs) == 0
        return 0.0
    end
    
    return maximum(clamp.(inputs, 0.0, 1.0))
end

"""
Soft maximum aggregator (smooth approximation to max).
A = (Σ xᵢᵏ)^(1/k) where k is large
"""
function soft_maximum_aggregator(inputs::Vector{Float64}; k::Float64=10.0)::Float64
    if length(inputs) == 0
        return 0.0
    end
    
    # Compute soft max: (Σ xᵢᵏ)^(1/k)
    sum_powers = 0.0
    for x in inputs
        x_safe = clamp(x, 0.0, 1.0)
        sum_powers += x_safe^k
    end
    
    result = sum_powers^(1.0/k)
    return clamp(result, 0.0, 1.0)
end

"""
Hill-based inhibition aggregator.
H = ∏ⱼ [1 / (1 + βⱼ * xⱼ^mⱼ)]

Parameters:
- inhibitors: vector of inhibitor signals xⱼ ∈ [0,1]
- betas: vector of inhibition strengths βⱼ ≥ 0
- ms: vector of Hill coefficients mⱼ > 0
"""
function inhibition_aggregator(
    inhibitors::AbstractVector{T},
    betas::AbstractVector{<:Real},
    ms::AbstractVector{<:Real},
) where {T<:Real}

    if !(length(inhibitors) == length(betas) == length(ms))
        error("All parameter vectors must have same length")
    end

    if length(inhibitors) == 0
        return one(T)  # No inhibition
    end

    result = one(T)
    for i in eachindex(inhibitors)
        if betas[i] <= 0
            continue  # neutral, no inhibition from this input
        end
        x_safe = clamp(inhibitors[i], zero(T), one(T))
        result *= one(T) / (one(T) + betas[i] * x_safe^ms[i])
    end

    return clamp(result, zero(T), one(T))
end

"""
Substrate availability aggregator (soft AND).
L = exp(Σ uₖ * log(xₖ + ε))

Similar to geometric mean but for substrate availability.
"""
function substrate_availability_aggregator(
    substrates::AbstractVector{T},
    weights::AbstractVector{<:Real};
    ε::Real=1e-12,
) where {T<:Real}

    # Same implementation as activator geometric mean (spec 2.2 substrate gate).
    return geometric_mean_aggregator(substrates, weights; ε=ε)
end

"""
Choose appropriate aggregator based on logic type and gate parameters.
"""
function apply_logic_aggregation(
    inputs::Vector{Float64},
    is_and_gate::Bool,
    weights::Union{Vector{Float64}, Nothing} = nothing
)::Float64
    
    if length(inputs) == 0
        return is_and_gate ? 1.0 : 0.0
    end
    
    if is_and_gate
        if weights !== nothing
            return geometric_mean_aggregator(inputs, weights)
        else
            return product_aggregator(inputs)
        end
    else  # OR gate
        if weights !== nothing
            # Weighted soft maximum
            return soft_maximum_aggregator(inputs)
        else
            return maximum_aggregator(inputs)
        end
    end
end

