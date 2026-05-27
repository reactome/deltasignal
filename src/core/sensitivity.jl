# Sensitivity Transform Functions
#
# Per spec section 2.1: each incoming activator signal x can be reshaped by a
# Hill-based sensitivity exponent α(x) before aggregation. Default is neutral
# (s=0 ⇒ α=1 ⇒ x̃ = x); learned per-edge values produce switch-like or graded
# behavior. Types are generic over Real so autodiff (ForwardDiff Duals) propagates.

"""
Sensitivity exponent α(x) = 1 + s · xⁿ/(xⁿ + K_α^n).
"""
function sensitivity_exponent(x::Real; s::Real=0.0, n::Real=2.0, K_α::Real=0.1)
    if s <= 0
        return one(typeof(x))
    end
    x_safe = clamp(x, 1e-12, one(typeof(x)))
    x_n = x_safe^n
    K_n = K_α^n
    return one(typeof(x)) + s * x_n / (x_n + K_n)
end

"""
Transformed input: x̃ = x^α(x). Identity when sensitivity is neutral (s=0).
"""
function apply_sensitivity_transform(x::Real; s::Real=0.0, n::Real=2.0, K_α::Real=0.1)
    if s <= 0
        return x
    end
    T = typeof(x)
    x_safe = clamp(x, T(1e-12), one(T))
    α = sensitivity_exponent(x_safe; s=s, n=n, K_α=K_α)
    return clamp(x_safe^α, zero(T), one(T))
end
