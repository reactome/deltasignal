# Hill Function Operations
#
# Spec section 2.3: reaction output y_r = s^h / (s^h + K^h). The current
# solver does NOT apply this Hill on the reaction's output (see
# reaction_model.jl — using linear identity for baseline stability), but
# the primitive is kept for explicit per-reaction use after training and
# for the inhibition aggregator (a Hill-suppression variant lives in
# aggregators.jl directly).

"""
Standard Hill function. y = s^h / (s^h + K^h). Generic-typed for autodiff.
"""
function hill_activation(s::Real, h::Real, K::Real)
    s_safe = clamp(s, 0, 1)
    K_safe = clamp(K, 1e-6, 1)
    h_safe = max(h, 1)

    s_h = s_safe^h_safe
    K_h = K_safe^h_safe

    return clamp(s_h / (s_h + K_h), 0, 1)
end
