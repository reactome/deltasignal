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

function fixture(; pools = :default, rl = identity)
    ids = ["G", "P", "A", "B", "EF", "EB", "F", "R", "W", "D"]
    nodes = Dict(rl(k) => DS.NetworkNode(rl(k), "R-" * k, "protein", nothing, rl(k), BL) for k in ids)
    E = [("G", "P", "input", true), ("P", "A", "output", false),
         ("A", "F", "input", true), ("EF", "F", "catalyst", true), ("F", "B", "output", false),
         ("B", "R", "input", true), ("EB", "R", "catalyst", true), ("R", "A", "output", false),
         ("B", "W", "input", true), ("W", "D", "output", false)]
    edges = [edge(rl(s), rl(t), et; and = a) for (s, t, et, a) in E]
    p = pools === :default ?
        [DS.CyclePool("pool1", sort([rl("A"), rl("B")]), rl("A"),
                      [(rl("A"), rl("B"), rl("F")), (rl("B"), rl("A"), rl("R"))])] : pools
    DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                       Dict{String, Set{String}}(), nothing, p)
end

function solve(net; mode = "balance", phi = "0.1", obs = Dict{String, Float64}(), rl = identity)
    with_env("DS_CYCLE_MODE" => mode, "DS_CYCLE_PHI" => phi, "DS_SELF_INHIBITOR_WEIGHT" => "off") do
        o = Dict(rl(k) => (1.0, 1.0) for k in ("G", "EF", "EB"))
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
