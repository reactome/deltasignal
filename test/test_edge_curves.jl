#!/usr/bin/env julia
#
# specs/128: response curves on edges (specs/125 stage 1). Off by default; a table at its identity setting must leave
# every solve unchanged; each curve is 1 at rest and 0 at a knockout, and bounds what one edge can pass.

using Test

include(joinpath(@__DIR__, "..", "src", "DeltaSignal.jl"))
using .DeltaSignal

const EX = joinpath(@__DIR__, "..", "examples")
const D = DeltaSignal

function with_env(f, pairs::Pair...)
    prev = Dict(k => get(ENV, k, nothing) for (k, _) in pairs)
    for (k, v) in pairs
        v === nothing ? delete!(ENV, k) : (ENV[k] = v)
    end
    try
        return f()
    finally
        for (k, v) in prev
            v === nothing ? delete!(ENV, k) : (ENV[k] = v)
        end
    end
end

function table(rows...)
    path = tempname() * ".tsv"
    write(path, join(["class\tK\th"; collect(rows)], "\n") * "\n")
    return path
end

"""Solve `edges` (source, target, positive, type) with `pins` (node => fold); returns node => fold and diagnostics."""
function solve(edges, pins; extra=String[])
    bl = D.resting_baseline()
    names = unique(vcat([e[1] for e in edges], [e[2] for e in edges], extra))
    nodes = Dict(u => D.NetworkNode(u, u, "protein", nothing, u, bl) for u in names)
    network = D.ReactionNetwork(nodes, [D.LogicNetworkEdge(s, t, true, pos, 1.0, ty) for (s, t, pos, ty) in edges],
                                Dict{String, D.SetExpansionMapping}())
    obs = Dict(u => (f * bl * 100.0, 1.0) for (u, f) in pins)
    r = D.solve_steady_state(network, obs, D.SteadyStateParams(1.0, 0.1, 500, 1e-9, "penalty"))
    d = copy(r.diagnostics); d["converged"] = r.converged
    return Dict(u => v / bl for (u, v) in r.node_activities), d
end

function sample_activities()
    network = D.parse_complete_network(joinpath(EX, "sample_logic_network.tsv"),
                                       joinpath(EX, "sample_uuid_mapping.tsv"),
                                       joinpath(EX, "sample_set_mappings.tsv"))
    ids = sort(collect(keys(network.nodes)))
    obs = Dict(ids[1] => (0.0, 1.0), ids[end] => (5.0, 1.0))   # a knockout and a 5x pin
    r = D.solve_steady_state(network, obs, D.SteadyStateParams(1.0, 0.1, 500, 1e-6, "penalty"))
    return r.node_activities, r.diagnostics
end

@testset "specs/128 edge curves" begin
    @testset "activating curve" begin
        for (K, h) in ((9.0, 1.0), (3.0, 2.0), (0.5, 1.5))
            @test D.curve_fold(1.0, K, h) ≈ 1.0
            @test D.curve_fold(0.0, K, h) == 0.0
            # elasticity at rest h K^h / (1 + K^h)
            δ = 1e-6
            e = (log(D.curve_fold(1 + δ, K, h)) - log(D.curve_fold(1 - δ, K, h))) / (log(1 + δ) - log(1 - δ))
            @test e ≈ h * K^h / (1 + K^h) rtol = 1e-5
            @test D.curve_fold(1e9, K, h) ≈ 1 + K^h rtol = 1e-6
            fs = 0.0:0.25:50.0
            @test issorted([D.curve_fold(f, K, h) for f in fs])
        end
        for f in (0.0, 0.3, 1.0, 7.0, 80.0)
            @test D.curve_fold(f, Inf, 1.0) === f          # the identity, exactly
        end
        @test D.curve_fold(4.0, Inf, 0.5) ≈ 2.0
    end

    @testset "inhibiting curve" begin
        @test D.curve_inhibit(1.0, 9.0, 1.0) ≈ 1.0
        @test D.curve_inhibit(0.0, 9.0, 1.0) ≈ 10.0
        @test D.curve_inhibit(2.0, 9.0, 1.0) ≈ 10 / 19
        @test D.curve_inhibit(0.0, 2.0, 2.0) ≈ 5.0
        @test D.curve_inhibit(1e9, 9.0, 1.0) < 1e-6
    end

    @testset "table parsing" begin
        @test D.edge_curves_env() === D.EDGE_CURVES_OFF
        c = D.parse_edge_curves(table("input\t9\t1", "regulator_neg\t9\t1", "catalyst\tinf\t1"))
        @test c.on && c.K == (9.0, Inf, Inf) && c.h == (1.0, 1.0, 1.0) && c.inh == (9.0, 1.0) && c.dep === nothing
        @test D.curve_active(c, D.SLOT_INPUT) && !D.curve_active(c, D.SLOT_CATALYST)
        @test !D.curve_active(c, D.SLOT_IDENTITY)
        c2 = D.parse_edge_curves(table("depletion\tdivide\t1"))
        @test c2.on && c2.dep === nothing
        for bad in (("input\t-1\t1",), ("input\t0\t1",), ("input\tx\t1",), ("output\t9\t1",),
                    ("input\t9\t1", "input\t3\t1"), ("regulator_neg\tdivide\t2",), ("input\t9",),
                    ("regulator_neg\tinf\t1",))
            @test_throws ArgumentError D.parse_edge_curves(table(bad...))
        end
        bad_header = tempname(); write(bad_header, "class\tk\tm\ninput\t9\t1\n")
        @test_throws ArgumentError D.parse_edge_curves(bad_header)
        @test_throws ArgumentError D.parse_edge_curves(tempname())
        with_env("DS_EDGE_CURVES" => table("regulator_neg\t9\t1"), "DS_INHIBITION_MODE" => "hill_sat") do
            @test_throws ArgumentError D.resolve_reaction_eval_config()
        end
        with_env("DS_EDGE_CURVES" => table("depletion\t9\t1"), "DS_DEPLETION_COMBINE" => "share") do
            @test_throws ArgumentError D.resolve_reaction_eval_config()
        end
    end

    @testset "slot classes" begin
        for (ty, cls) in (("input", D.SLOT_INPUT), ("assembly", D.SLOT_INPUT), ("composition", D.SLOT_INPUT),
                          ("catalyst", D.SLOT_CATALYST), ("regulator", D.SLOT_REGULATOR_POS),
                          ("output", D.SLOT_IDENTITY), ("variant_pool", D.SLOT_IDENTITY),
                          ("specialization", D.SLOT_IDENTITY), ("recycled_supply", D.SLOT_IDENTITY))
            @test D.slot_class(ty) == cls
        end
        bl = 0.01
        nodes = Dict(u => D.NetworkNode(u, u, "protein", nothing, u, bl) for u in ("A", "B", "C", "R"))
        net = D.ReactionNetwork(nodes, [D.LogicNetworkEdge("A", "R", true, true, 1.0, "input"),
                                        D.LogicNetworkEdge("B", "R", true, true, 1.0, "catalyst"),
                                        D.LogicNetworkEdge("C", "R", true, true, 1.0, "regulator")],
                                Dict{String, D.SetExpansionMapping}())
        r = only(D.convert_to_reaction_network(net))
        @test Dict(zip(r.activator_uuids, r.activator_class)) ==
              Dict("A" => D.SLOT_INPUT, "B" => D.SLOT_CATALYST, "C" => D.SLOT_REGULATOR_POS)
        # a merged input + catalyst slot takes the substrate curve
        with_env("DS_DEDUP_ACTIVATORS" => "1") do
            net2 = D.ReactionNetwork(nodes, [D.LogicNetworkEdge("A", "R", true, true, 1.0, "catalyst"),
                                             D.LogicNetworkEdge("A", "R", true, true, 1.0, "input")],
                                     Dict{String, D.SetExpansionMapping}())
            @test only(D.convert_to_reaction_network(net2)).activator_class == [D.SLOT_INPUT]
        end
    end

    @testset "off and the identity table leave a solve unchanged" begin
        base, d0 = sample_activities()
        @test d0["edge_curve_rule"] == "off" && isempty(d0["edge_curve_slots"])
        ident = table("input\tinf\t1", "catalyst\tinf\t1", "regulator_pos\tinf\t1",
                      "regulator_neg\tdivide\t1", "depletion\tdivide\t1")
        with_env("DS_EDGE_CURVES" => ident) do
            same, d1 = sample_activities()
            @test same == base
            @test d1["edge_curve_rule"] != "off" && isempty(d1["edge_curve_slots"])
        end
    end

    @testset "curves act on their class only" begin
        chain = [("A", "R", true, "input"), ("R", "T", true, "output"), ("C", "R", true, "catalyst")]
        off, _ = solve(chain, Dict("A" => 2.0, "C" => 1.0))
        @test off["T"] ≈ 2.0
        with_env("DS_EDGE_CURVES" => table("input\t9\t1")) do
            on, d = solve(chain, Dict("A" => 2.0, "C" => 1.0))
            @test on["R"] ≈ 2 * 10 / 11       # substrate below saturation passes part of the change
            @test on["T"] ≈ on["R"]            # the output edge is identity
            @test d["edge_curve_slots"] == Dict("input" => 1)
            ko, _ = solve(chain, Dict("A" => 0.0, "C" => 1.0))
            @test ko["T"] == 0.0               # a knockout stays a knockout
            big, _ = solve(chain, Dict("A" => 80.0, "C" => 1.0))
            @test big["T"] < 10.0              # at most 1 + K
            cat, _ = solve(chain, Dict("A" => 1.0, "C" => 2.0))
            @test cat["T"] ≈ 2.0               # the catalyst stays linear
        end
    end

    @testset "inhibitor and depletion curves" begin
        inh = [("U", "R", true, "input"), ("I", "R", false, "regulator")]
        with_env("DS_EDGE_CURVES" => table("regulator_neg\t9\t1")) do
            ko, _ = solve(inh, Dict("U" => 1.0, "I" => 0.0))
            @test ko["R"] ≈ 10.0
            up, _ = solve(inh, Dict("U" => 1.0, "I" => 2.0))
            @test up["R"] ≈ 10 / 19
        end
        # two inhibitors knocked out: each passes its own bound, so 100x, where divide's reaction cap read 10x
        inh2 = [("U", "R", true, "input"), ("I", "R", false, "regulator"), ("J", "R", false, "regulator")]
        div, _ = solve(inh2, Dict("U" => 1.0, "I" => 0.0, "J" => 0.0))
        @test div["R"] ≈ 10.0
        with_env("DS_EDGE_CURVES" => table("regulator_neg\t9\t1")) do
            two, _ = solve(inh2, Dict("U" => 1.0, "I" => 0.0, "J" => 0.0))
            @test two["R"] ≈ 100.0 rtol = 1e-6
        end
        dep = [("U", "S", true, "input"), ("E", "S", false, "depletion")]
        with_env("DS_EDGE_CURVES" => table("depletion\t9\t1")) do
            ko, _ = solve(dep, Dict("U" => 1.0, "E" => 0.0))
            @test ko["S"] ≈ 10.0
            hi, _ = solve(dep, Dict("U" => 1.0, "E" => 80.0))
            @test hi["S"] ≈ 0.1                # the depletion floor is kept
        end
    end

    @testset "a positive loop contracts" begin
        # U -> A -> B -> A: the loop's gain at rest is 1 today, below 1 with a substrate curve
        loop = [("U", "A", true, "input"), ("B", "A", true, "input"), ("A", "B", true, "input")]
        with_env("DS_EDGE_CURVES" => table("input\t3\t1")) do
            a, d = solve(loop, Dict("U" => 2.0))
            @test d["converged"]
            @test 1.0 < a["A"] < 4.0 && 1.0 < a["B"] < a["A"]
        end
    end
end
