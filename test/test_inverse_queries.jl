#!/usr/bin/env julia
#
# specs/158: inverse queries — which perturbations move a chosen readout.
# The central assertion is FORWARD-PATH FIDELITY: an inverse row's numbers
# must equal running the public /api/solve + /api/readouts for that candidate
# (the fixture runs both paths and compares). The fixture computes: A → B → C
# with D inhibiting C, so "what lowers C" is A's knockout and D's
# over-expression, and the test derives those from the free solve it also
# runs — never from a pinned expectation.

using Test
using JSON3

include(joinpath(@__DIR__, "..", "src", "DeltaSignal.jl"))
using .DeltaSignal
const SRV = DeltaSignal

const BL = 0.01
post(path, body) = SRV.HTTP.Request("POST", path, ["Content-Type" => "application/json"], body)
jread(resp) = JSON3.read(String(resp.body))

node(u) = SRV.NetworkNode(u, "R-$u", "protein", nothing, u, BL)
pedge(s, t) = SRV.LogicNetworkEdge(s, t, true, true, 1.0, "input")
nedge(s, t) = SRV.LogicNetworkEdge(s, t, true, false, 1.0, "input")

# A → B → C, D ⊣ C. Roots: A and D (B and C are produced, so the pin
# protocol cannot perturb them — exactly the candidate rule under test).
const NET_ID = "inverse-fixture"
network = SRV.ReactionNetwork(
    Dict(u => node(u) for u in ("A", "B", "C", "D")),
    [pedge("A", "B"), pedge("B", "C"), nedge("D", "C")],
    Dict{String, SRV.SetExpansionMapping}())
SRV.cache_network!(NET_ID, network)
SRV.cache_pin_index!(NET_ID, SRV.PinIndex(
    Dict("R-A" => ["A"], "R-B" => ["B"], "R-C" => ["C"], "R-D" => ["D"]),
    Dict{String, String}()))
SRV.cache_readout_index!(NET_ID, SRV.ReadoutIndex(
    Dict{String, Vector{String}}(), Set{String}(), false, false, String[], nothing))

inverse(body) = jread(SRV.inverse_handler(post("/api/inverse", JSON3.write(body))))
target_c = Dict("stable_ids" => ["R-C"])

@testset "inverse queries (specs/158)" begin

    @testset "the fixture's physics, from the free solve" begin
        free = SRV.solve_handler(post("/api/solve",
            JSON3.write(Dict("network_id" => NET_ID, "observations" => Dict()))))
        @test jread(free)[:status] == "success"
    end

    @testset "what lowers C: A knockout and D over-expression, forward-equal" begin
        r = inverse(Dict("network_id" => NET_ID, "target" => target_c,
                         "direction" => "down"))
        @test r[:status] == "success"
        got = Dict((String(row[:candidates][1]), String(row[:dose])) => row for row in r[:results])
        @test haskey(got, ("R-A", "knockout"))
        @test haskey(got, ("R-D", "overexpression"))
        @test !any(k -> k[1] == "R-B" || k[1] == "R-C", keys(got))
        for ((cid, dose), row) in got
            @test String(row[:call]) == "down"
            # FR-001: bit-equality with the public forward path
            pin = String(row[:pins][1])
            act = Float64(row[:activity])
            sol = jread(SRV.solve_handler(post("/api/solve", JSON3.write(Dict(
                "network_id" => NET_ID,
                "observations" => Dict(pin => [act, 1.0]))))))
            ro = jread(SRV.readouts_handler(post("/api/readouts", JSON3.write(Dict(
                "network_id" => NET_ID, "activities" => sol[:node_activities],
                "pinned" => [pin], "converged" => sol[:converged],
                "queries" => [Dict("id" => "t", "stable_ids" => ["R-C"])])))))
            @test Float64(row[:fold]) == Float64(ro[:results][1][:fold])
            @test String(row[:call]) == String(ro[:results][1][:call])
        end
        # counts are honest: A and D are candidates; B and C have no root pin
        @test Int(r[:counts][:candidates_total]) == 2
        @test Int(r[:counts][:reachable]) == 2
        @test Int(r[:counts][:solved]) == 2
        @test Int(r[:counts][:capped]) == 0
        @test String(r[:ordering]) == "entity_id"
    end

    @testset "what raises C: the opposite doses" begin
        r = inverse(Dict("network_id" => NET_ID, "target" => target_c,
                         "direction" => "up"))
        got = Set((String(row[:candidates][1]), String(row[:dose])) for row in r[:results])
        @test ("R-A", "overexpression") in got
        @test ("R-D", "knockout") in got
        @test all(String(row[:call]) == "up" for row in r[:results])
    end

    @testset "direction any returns every change, ranked by effect" begin
        r = inverse(Dict("network_id" => NET_ID, "target" => target_c))
        @test length(r[:results]) >= 4
        effects = [abs(log2(Float64(row[:fold]))) for row in r[:results]]
        @test issorted(effects; rev = true)
    end

    @testset "self-exclusion and the honest zero" begin
        # target D: D itself is a candidate but is the target (excluded and
        # counted); A has no path to D, so nothing is solvable — an honest
        # empty success, never an error.
        r = inverse(Dict("network_id" => NET_ID,
                         "target" => Dict("stable_ids" => ["R-D"])))
        @test r[:status] == "success"
        @test Int(r[:counts][:self_excluded]) == 1
        @test Int(r[:counts][:reachable]) == 0
        @test Int(r[:counts][:solved]) == 0
        @test isempty(r[:results])
    end

    @testset "the cap is honest" begin
        r = inverse(Dict("network_id" => NET_ID, "target" => target_c, "cap" => 1))
        @test Int(r[:counts][:solved]) == 1
        @test Int(r[:counts][:capped]) == 1
    end

    @testset "refusals are specific" begin
        bad(body) = SRV.inverse_handler(post("/api/inverse", JSON3.write(body)))
        r = bad(Dict("network_id" => NET_ID, "target" => target_c, "cap" => 0))
        @test r.status == 400
        @test jread(r)[:message] == SRV.INVERSE_CAP_ERROR
        r = bad(Dict("network_id" => NET_ID, "target" => target_c, "cap" => 501))
        @test jread(r)[:message] == SRV.INVERSE_CAP_ERROR
        r = bad(Dict("network_id" => NET_ID))
        @test jread(r)[:message] == SRV.INVERSE_QUERY_ERROR
        r = bad(Dict("network_id" => NET_ID,
                     "target" => Dict("stable_ids" => ["R-NOWHERE"])))
        @test jread(r)[:message] == SRV.INVERSE_TARGET_ERROR
        r = bad(Dict("network_id" => "no-such-network", "target" => target_c))
        @test jread(r)[:message] == SRV.UNKNOWN_NETWORK_ERROR
    end
end

println("test_inverse_queries.jl passed")
