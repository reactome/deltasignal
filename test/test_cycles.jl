# specs/039: interconversion cycles as conserved pools, solved as pi = pi P
# (docs/MODEL.md §4). Fixture = a two-form pool A ⇄ B (RAS-GDP ⇄ RAS-GTP):
# G -> P -> A supplies the protein; F: A -> B catalysed by EF (a GEF); R: B -> A
# catalysed by EB (a GAP); B -> W -> D is a downstream readout. Expected folds
# at phi0 = 0.1 are the specs/039 fixtures (derived independently, twice).
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
node(u) = u => DS.NetworkNode(u, "R-" * u, "protein", nothing, u, BL)
edge(s, t, et; and = true) = DS.LogicNetworkEdge(s, t, and, true, 1.0, et)

function fixture(; pools = :default, rl = identity, extra = [])
    ids = ["G", "P", "A", "B", "EF", "EB", "F", "R", "W", "D", "Z", "Y", "F2", "K", "Q", "M"]
    nodes = Dict(rl(k) => DS.NetworkNode(rl(k), "R-" * k, "protein", nothing, rl(k), BL) for k in ids)
    E = [("G", "P", "input", true), ("P", "A", "output", false),
         ("A", "F", "input", true), ("EF", "F", "catalyst", true), ("F", "B", "output", false),
         ("B", "R", "input", true), ("EB", "R", "catalyst", true), ("R", "A", "output", false),
         ("B", "W", "input", true), ("W", "D", "output", false)]
    edges = [edge(rl(s), rl(t), et; and = a) for (s, t, et, a) in E]
    for (s, t, et, pos) in extra
        push!(edges, DS.LogicNetworkEdge(rl(s), rl(t), true, pos, 1.0, et))
    end
    p = pools === :default ?
        [DS.CyclePool("pool1", sort([rl("A"), rl("B")]), rl("A"),
                      [(rl("A"), rl("B"), rl("F")), (rl("B"), rl("A"), rl("R"))])] : pools
    DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                       Dict{String, Set{String}}(), nothing, p)
end

function solve(net; mode = "balance", phi = "0.1", obs = Dict{String, Float64}(), rl = identity,
               method = nothing)
    with_env("DS_CYCLE_MODE" => mode, "DS_CYCLE_PHI" => phi, "DS_SELF_INHIBITOR_WEIGHT" => "off",
             "DS_SCC_METHOD" => method) do
        o = Dict(rl(k) => (1.0, 1.0) for k in ("G", "EF", "EB", "Z"))
        for (k, v) in obs; o[rl(k)] = (v, 1.0); end
        DS.solve_steady_state(net, o, DS.SteadyStateParams(1.0, 0.1, 500, 1e-9, "penalty"))
    end
end
fold(r, u) = r.node_activities[u] / BL

@testset "interconversion cycles (specs/039)" begin

@testset "the specs/039 fixtures at phi0 = 0.1" begin
    net = fixture()
    for (obs, a, b, f) in ((Dict{String,Float64}(), 1.0, 1.0, 1.0),
                           (Dict("EF" => 80.0), 0.11236, 8.98876, 8.98876),
                           (Dict("EF" => 0.0), 1.11111, 0.0, 0.0),
                           (Dict("EB" => 0.0), 0.0, 10.0, 0.0),
                           (Dict("EB" => 80.0), 1.10957, 0.01387, 1.10957),
                           (Dict("G" => 80.0), 80.0, 80.0, 80.0))
        r = solve(net; obs = obs)
        @test r.converged
        @test fold(r, "A") ≈ a atol = 1e-4
        @test fold(r, "B") ≈ b atol = 1e-4
        @test fold(r, "F") ≈ f atol = 1e-4
        @test fold(r, "D") ≈ b atol = 1e-4          # the readout follows the modified form
        @test r.diagnostics["cycle_rule"] == "balance"
        @test r.diagnostics["cycle_pools_solved"] == 1
    end
end

@testset "phi0 = 0.5 compresses, as specs/039 warns" begin
    r = solve(fixture(); phi = "0.5", obs = Dict("EF" => 80.0))
    @test fold(r, "B") ≈ 2 * 80 / 81 atol = 1e-4      # s·r/(φ r + 1 − φ)
end

@testset "default off, and a missing table is reported" begin
    with_env("DS_CYCLE_MODE" => nothing) do
        @test DS.cycle_mode() == "off"
    end
    r = solve(fixture(); mode = "off")
    @test r.diagnostics["cycle_rule"] == "off"
    @test r.diagnostics["cycle_pools_solved"] == 0
    r = solve(fixture(; pools = nothing))
    @test r.diagnostics["cycle_rule"] == "balance: no pool table"
    @test r.diagnostics["cycle_pools_solved"] == 0
end

@testset "a pinned form sets the pool" begin
    r = solve(fixture(); obs = Dict("B" => 5.0))
    @test fold(r, "B") ≈ 5.0 atol = 1e-9
    @test fold(r, "A") ≈ 5.0 atol = 1e-4              # unperturbed split: supply 5x
end

@testset "review fixes (specs/039 amendment 1)" begin
    # finding 5: every drive at zero leaves the split at pi0, not wherever a
    # constant added to Q puts it
    r = solve(fixture(); obs = Dict("EF" => 0.0, "EB" => 0.0))
    @test r.converged
    @test fold(r, "A") ≈ 1.0 atol = 1e-6
    @test fold(r, "B") ≈ 1.0 atol = 1e-6
    # finding 6: a producer of a form inside the pool's own component is not
    # supply. Y: B -> A is such a producer; the supply stays P's.
    net = fixture(; extra = [("B", "Y", "input", true), ("Y", "A", "output", true)])
    r = solve(net; obs = Dict("G" => 80.0))
    @test r.converged
    @test fold(r, "A") ≈ 80.0 atol = 1e-4
    @test fold(r, "B") ≈ 80.0 atol = 1e-4
    # finding 4: a depletion edge into a form from outside the pool still acts
    net = fixture(; extra = [("Z", "B", "depletion", false)])
    r0 = solve(net)
    @test fold(r0, "B") ≈ 1.0 atol = 1e-6               # Z at baseline: no effect
    r = solve(net; obs = Dict("Z" => 10.0))
    @test r.converged
    @test fold(r, "B") < 0.99
    @test fold(r, "D") ≈ fold(r, "B") atol = 1e-6         # and the readout sees it
    # ... but the pool's own catalyst ⊣ source-form depletion is not applied
    # again (pre-registered: the balance already contains it)
    net = fixture(; extra = [("EF", "A", "depletion", false)])
    r = solve(net; obs = Dict("EF" => 80.0))
    @test fold(r, "A") ≈ 0.11236 atol = 1e-4
    @test fold(r, "B") ≈ 8.98876 atol = 1e-4
    # Fable finding 5: an uncatalysed step beside a catalysed one for the same
    # pair carries a negligible baseline share, so GEF knockout empties the
    # modified form instead of halving it
    xtra = [("A", "F2", "input", true), ("F2", "B", "output", true)]
    pools2(uncat) = [DS.CyclePool("pool1", ["A", "B"], "A",
                                  sort([("A", "B", "F"), ("A", "B", "F2"), ("B", "A", "R")]), uncat)]
    r = solve(fixture(; extra = xtra, pools = pools2(Set(["F2"]))); obs = Dict("EF" => 0.0))
    @test r.converged
    @test fold(r, "B") < 0.01
    r = solve(fixture(; extra = xtra, pools = pools2(Set{String}())); obs = Dict("EF" => 0.0))
    @test 0.4 < fold(r, "B") < 0.7                       # without the flag: the old equal split
    r = solve(fixture(; extra = xtra, pools = pools2(Set(["F2"]))))
    @test fold(r, "B") ≈ 1.0 atol = 1e-6                 # baseline is still exact
    # Fable finding 2: a producer that shares the pool's component only through
    # a pinned cofactor (K, released by R and consumed by M) is still supply
    net = fixture(; extra = [("R", "K", "output", true), ("K", "M", "input", true),
                             ("Q", "M", "input", true), ("M", "A", "output", true)])
    r = solve(net; obs = Dict("K" => 1.0, "Q" => 80.0))
    @test r.converged
    @test fold(r, "A") ≈ 40.5 atol = 1e-3                # mean of P (1x) and M (80x)
    @test fold(r, "B") ≈ 40.5 atol = 1e-3
    # finding 9: a component method that bypasses the sweep manages no pool,
    # and says so
    for m in ("pool", "minimize")
        r = solve(fixture(); method = m)
        @test r.diagnostics["cycle_pools_solved"] == 0
        @test r.diagnostics["cycle_pools_unmanaged"] == 1
    end
end

@testset "label-independent" begin
    rl(u) = "zz_" * u
    r1 = solve(fixture(); obs = Dict("EF" => 3.0))
    r2 = solve(fixture(; rl = rl); obs = Dict("EF" => 3.0), rl = rl)
    for k in ("A", "B", "F", "R", "D")
        @test r1.node_activities[k] ≈ r2.node_activities[rl(k)] atol = 1e-12
    end
end

@testset "configuration errors" begin
    for bad in ("Balance", "on", "")
        with_env("DS_CYCLE_MODE" => bad) do
            @test_throws ArgumentError DS.cycle_mode()
        end
    end
    for bad in ("0", "1", "1.5", "x")
        with_env("DS_CYCLE_PHI" => bad) do
            @test_throws ArgumentError DS.cycle_phi()
        end
    end
end

@testset "pools.csv parsing" begin
    mktempdir() do dir
        ln = joinpath(dir, "logic_network.csv"); write(ln, "")
        @test DS.parse_pools(ln) === nothing
        write(joinpath(dir, "pools.csv"), "pool_id,node_uuid,stable_id,is_base\npool1,u-b,R-B,False\npool1,u-a,R-A,True\n")
        write(joinpath(dir, "pool_transitions.csv"),
              "pool_id,from_uuid,to_uuid,reaction_uuid,reaction_stid\npool1,u-a,u-b,u-f,R-F\npool1,u-b,u-a,u-r,R-R\n")
        p = DS.parse_pools(ln)
        @test length(p) == 1 && p[1].base == "u-a" && p[1].forms == ["u-a", "u-b"]
        @test p[1].transitions == [("u-a", "u-b", "u-f"), ("u-b", "u-a", "u-r")]
        write(joinpath(dir, "pools.csv"), "pool_id,node_uuid,stable_id,is_base\npool1,u-a,R-A,False\n")
        @test_throws ArgumentError DS.parse_pools(ln)   # a pool must have a base form
    end
end

end
