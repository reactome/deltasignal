# specs/066: limit cycles (DS_NONCONV_FALLBACK). A Goodwin ring, the shape
# traced in Fanconi Anemia: a pinned subunit P and nine free subunits S1..S9
# AND-assemble (rA) into C1, which passes along a chain of single-input
# reactions to D; D depletes every free subunit. Nine depletions multiplied by
# one AND give a loop gain of 9 across nine lags, so the fixed point repels
# damped relaxation and the iteration cycles. With hill_sat's exact fold
# product and divide-form depletion (s = 1/d), the fixed point is
# d = P * s^9 = P * d^-9, i.e. d = P^(1/10) and s = P^(-1/10).
using Test
using DeltaSignal
const DS = DeltaSignal

function with_env(f, pairs...)
    old = Dict(k => get(ENV, k, nothing) for (k, _) in pairs)
    for (k, v) in pairs
        v === nothing ? delete!(ENV, k) : (ENV[k] = v)
    end
    try
        f()
    finally
        for (k, v) in old
            v === nothing ? delete!(ENV, k) : (ENV[k] = v)
        end
    end
end

const BL = 0.01
const SUBS = ["S$i" for i in 1:9]
const CHAIN = ["C1", "C2", "C3", "D"]
e(s, t, et = "input") = DS.LogicNetworkEdge(s, t, et != "output", et != "depletion", 1.0, et)

# ring = false drops the depletion edges: the same nodes, no cycle.
# deplete = k keeps only the first k depletion edges (loop gain k).
function ring(; ring = true, deplete = 9)
    names = vcat(["P"], SUBS, CHAIN, ["rA", "rB", "rC", "rD"])
    nodes = Dict(u => DS.NetworkNode(u, "R-$u", "protein", nothing, u, BL) for u in names)
    edges = [e("P", "rA"); [e(s, "rA") for s in SUBS]; e("rA", "C1", "output");
             e("C1", "rB"); e("rB", "C2", "output"); e("C2", "rC"); e("rC", "C3", "output");
             e("C3", "rD"); e("rD", "D", "output")]
    ring && append!(edges, [e("D", s, "depletion") for s in SUBS[1:deplete]])
    DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                       Dict{String, Set{String}}(), nothing, nothing, nothing, nothing, nothing)
end
# p is on the UI scale (80 = 80x baseline).
solve(net, mode; p = 80.0, iters = 500) =
    with_env("DS_NONCONV_FALLBACK" => mode, "DS_DEPLETION_OWN_PRODUCT" => "full") do
        DS.solve_steady_state(net, Dict("P" => (p, 1.0)), DS.SteadyStateParams(1.0, 0.1, iters, 1e-6, "penalty"))
    end
fold(r, u) = r.node_activities[u] / BL

@testset "limit cycles (specs/066)" begin

@testset "default is off; bad values are errors" begin
    with_env("DS_NONCONV_FALLBACK" => nothing) do
        @test DS.nonconv_fallback_mode() == "off"
    end
    with_env("DS_NONCONV_FALLBACK" => "LM") do
        @test_throws ArgumentError DS.nonconv_fallback_mode()
    end
    with_env("DS_NONCONV_FALLBACK" => "lm", "DS_NONCONV_LM_ITERS" => "0") do
        @test_throws ArgumentError DS._nonconv_lm_iters()
    end
end

@testset "the ring is a limit cycle under damped relaxation" begin
    r500 = solve(ring(), "off"; iters = 500)
    r501 = solve(ring(), "off"; iters = 501)
    @test !r500.converged
    @test r500.diagnostics["nonconv_fallback"] == "off"
    @test r500.diagnostics["nonconv_components"] == 0      # off counts nothing
    # the scored value depends on where the cap stops
    @test abs(log(fold(r500, "D")) - log(fold(r501, "D"))) > 1e-3
end

@testset "lm finds the fixed point, independent of the cap" begin
    r = solve(ring(), "lm"; iters = 500)
    r2 = solve(ring(), "lm"; iters = 2000)
    @test r.converged
    @test r.diagnostics["nonconv_fallback"] == "lm"
    @test r.diagnostics["nonconv_components"] == 1
    @test r.diagnostics["nonconv_improved"] == 1
    @test r.diagnostics["nonconv_rescued"] == 1
    @test r.diagnostics["nonconv_skipped_pools"] == 0
    @test isapprox(fold(r, "D"), 80.0^(1 / 10); rtol = 1e-4)
    for s in SUBS
        @test isapprox(fold(r, s), 80.0^(-1 / 10); rtol = 1e-4)
    end
    @test isapprox(fold(r2, "D"), fold(r, "D"); rtol = 1e-6)
    # a knockout: the AND reads 0 whatever the subunits do (each is capped at
    # the 10x de-repression ceiling), so D falls to hill_sat's floor.
    k = solve(ring(), "lm"; p = 0.0)
    @test fold(k, "D") < 1e-3
end

@testset "converged solves are untouched" begin
    # no cycle: identical under off and lm, bit for bit
    a = solve(ring(ring = false), "off")
    b = solve(ring(ring = false), "lm")
    @test a.converged && b.converged
    @test a.node_activities == b.node_activities
    @test b.diagnostics["nonconv_components"] == 0
    # a ring of gain 1 (one depletion edge) relaxes: the fallback never runs
    c = solve(ring(deplete = 1), "lm")
    d = solve(ring(deplete = 1), "off")
    @test d.converged
    @test c.node_activities == d.node_activities
    @test c.diagnostics["nonconv_components"] == 0
end

end
