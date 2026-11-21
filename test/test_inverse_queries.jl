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

# A row's own fields; every OTHER field must be the forward readout entry's,
# and every forward field except the target's resolution must be in the row.
const ROW_OWN = Set([:candidates, :names, :pins, :dose, :activity, :converged])
const TARGET_FIELDS = Set([:id, :mode, :basis, :groups, :pinned])
function same_fields(row, fwd)
    for (k, v) in pairs(fwd)
        k in TARGET_FIELDS && continue
        @test haskey(row, k)
        haskey(row, k) && @test row[k] == v
    end
    for k in keys(row)
        k in ROW_OWN && continue
        @test haskey(fwd, k)
    end
end

"""The public forward path for one pin at one dose: /api/solve, then
/api/readouts with whatever the solve returned (switch_components passed
only when present), exactly as a client runs it."""
function forward(nid, pin_uuids, act, target)
    sol = jread(SRV.solve_handler(post("/api/solve", JSON3.write(Dict(
        "network_id" => nid,
        "observations" => Dict(u => [act, 1.0] for u in pin_uuids))))))
    body = Dict{String, Any}("network_id" => nid, "activities" => sol[:node_activities],
                             "pinned" => pin_uuids, "converged" => sol[:converged],
                             "queries" => [merge(Dict("id" => "t"), target)])
    haskey(sol, :switch_components) && (body["switch_components"] = sol[:switch_components])
    ro = jread(SRV.readouts_handler(post("/api/readouts", JSON3.write(body))))
    return sol, ro[:results][1]
end

# Review 2026-11-20 (C1): the specs/059 recycled-twin shape. Entity X is split
# into a ROOT copy R and a RECYCLED copy Y (recruited by rA into P, released
# back to Y by rD); P signals through rT to the readout T. Over the PARSED
# edges R reaches only its side reaction (rS -> C2), so the old pre-filter
# skipped R; the solver adds R -> Y (`recycled_supply`) and T follows R.
const TWIN_ID = "inverse-twin-fixture"
tnode(u, rid) = u => SRV.NetworkNode(u, rid, "protein", nothing, u, BL)
tedge(s, t, et = "input") = SRV.LogicNetworkEdge(s, t, et != "output", true, 1.0, et)
twin = SRV.ReactionNetwork(
    Dict(tnode("R", "R-X"), tnode("Y", "R-X"), tnode("P", "R-PX"), tnode("T", "R-T"),
         tnode("C2", "R-C2"), tnode("rS", "R-rS"), tnode("rA", "R-rA"), tnode("rD", "R-rD"),
         tnode("rT", "R-rT")),
    [tedge("R", "rS"), tedge("rS", "C2", "output"),
     tedge("Y", "rA"), tedge("rA", "P", "output"), tedge("P", "rD"), tedge("rD", "Y", "output"),
     tedge("P", "rT"), tedge("rT", "T", "output")],
    Dict{String, SRV.SetExpansionMapping}())
SRV.cache_network!(TWIN_ID, twin)
SRV.cache_pin_index!(TWIN_ID, SRV.PinIndex(
    Dict("R-X" => ["R", "Y"], "R-PX" => ["P"], "R-T" => ["T"], "R-C2" => ["C2"]),
    Dict{String, String}()))
SRV.cache_readout_index!(TWIN_ID, SRV.ReadoutIndex(
    Dict{String, Vector{String}}(), Set{String}(), false, false, String[], nothing))

# Review 2026-11-20 (C4): the first fixture plus a ROOT reaction node RX
# (node kind `reaction`) feeding C, whose stable id the pin index resolves
# (as nodes.csv indexes a reaction by its diagram_entity_id), and an entity
# id E that (artificially) resolves to RX: neither may become a candidate.
const RXN_ID = "inverse-reaction-fixture"
rxnet = SRV.ReactionNetwork(
    Dict(u => node(u) for u in ("A", "B", "C", "D", "RX", "E")),
    [pedge("A", "B"), pedge("B", "C"), nedge("D", "C"),
     SRV.LogicNetworkEdge("RX", "C", false, true, 1.0, "output")],
    Dict{String, SRV.SetExpansionMapping}(), Set{String}(), Dict{String, Set{String}}(),
    nothing, nothing, nothing, nothing,
    Dict{String, SRV.NodeLeaves}("A" => ("simple_entity", ["R-A"]), "B" => ("simple_entity", ["R-B"]),
                                 "C" => ("simple_entity", ["R-C"]), "D" => ("simple_entity", ["R-D"]),
                                 "E" => ("simple_entity", ["R-E"]), "RX" => ("reaction", String[])))
SRV.cache_network!(RXN_ID, rxnet)
SRV.cache_pin_index!(RXN_ID, SRV.PinIndex(
    Dict("R-A" => ["A"], "R-B" => ["B"], "R-C" => ["C"], "R-D" => ["D"],
         "R-RX" => ["RX"], "R-E" => ["RX"]),
    Dict{String, String}()))
SRV.cache_readout_index!(RXN_ID, SRV.ReadoutIndex(
    Dict{String, Vector{String}}(), Set{String}(), false, false, String[], nothing))

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
            _, fwd = forward(NET_ID, [pin], act, target_c)
            @test Float64(row[:fold]) == Float64(fwd[:fold])
            @test String(row[:call]) == String(fwd[:call])
            # review 2026-11-20 (C2): EVERY served field, honesty fields
            # included, not only fold and call
            same_fields(row, fwd)
        end
        # counts are honest: A and D are candidates; B and C have no root pin
        @test Int(r[:counts][:candidates_total]) == 2
        @test Int(r[:counts][:reachable]) == 2
        @test Int(r[:counts][:solved]) == 2
        @test Int(r[:counts][:capped]) == 0
        @test String(r[:ordering]) == "entity_id"
        @test String(r[:reachability]) == "solve_edges"
        # review 2026-11-20: `moved` counts CANDIDATES — A and D each move C at
        # both doses, which is 2 candidates (it used to report 4 solves)
        @test Int(r[:counts][:moved]) == 2
        @test Int(r[:counts][:reaction_ids_excluded]) == 0
    end

    @testset "C1: reachability reads the edges the solve adds" begin
        # the parsed edges alone give R no path to T (the old pre-filter)
        @test !any(e -> e.parent_uuid == "R" && e.child_uuid == "Y", twin.edges)
        # the solver's own network has the recycled_supply edge
        @test any(e -> e.parent_uuid == "R" && e.child_uuid == "Y" && e.edge_type == "recycled_supply",
                  SRV.solved_carrier_network(twin).edges)
        tgt = Dict("stable_ids" => ["R-T"])
        # what the forward path says R's knockout does to T, computed, not pinned
        _, fwd = forward(TWIN_ID, ["R"], 0.0, tgt)
        @test String(fwd[:call]) == "down"
        r = inverse(Dict("network_id" => TWIN_ID, "target" => tgt, "direction" => "down"))
        @test r[:status] == "success"
        ko = [row for row in r[:results] if "R-X" in row[:candidates] && row[:dose] == "knockout"]
        @test length(ko) == 1
        if length(ko) == 1
            @test collect(String.(ko[1][:pins])) == ["R"]
            same_fields(ko[1], fwd)
        end
        # candidates: R (root) and P (root_cycle); T and C2 are produced dead
        # ends, so they pin nothing
        @test Int(r[:counts][:candidates_total]) == 2
        @test Int(r[:counts][:reachable]) == 2
    end

    @testset "C2: an absent switch_components stays absent in every row" begin
        # with detection off the forward path serves no switch field at all;
        # the inverse rows must not invent switch:false
        withenv("DS_SWITCH_DETECT" => "off") do
            r = inverse(Dict("network_id" => NET_ID, "target" => target_c))
            @test !isempty(r[:results])
            for row in r[:results]
                @test !haskey(row, :switch)
                sol, fwd = forward(NET_ID, String.(row[:pins]), Float64(row[:activity]), target_c)
                @test !haskey(sol, :switch_components)
                same_fields(row, fwd)
            end
        end
        # with detection on (the default) and a converged solve, both paths
        # carry the switch verdict, equal
        withenv("DS_SWITCH_DETECT" => "on") do
            r = inverse(Dict("network_id" => NET_ID, "target" => target_c))
            for row in r[:results]
                sol, fwd = forward(NET_ID, String.(row[:pins]), Float64(row[:activity]), target_c)
                @test haskey(row, :switch) == haskey(sol, :switch_components)
                same_fields(row, fwd)
            end
        end
    end

    @testset "C4: a reaction is never a candidate, and is counted" begin
        r = inverse(Dict("network_id" => RXN_ID, "target" => target_c))
        @test r[:status] == "success"
        cands = Set(String(c) for row in r[:results] for c in row[:candidates])
        @test !("R-RX" in cands) && !("R-E" in cands)
        @test !any(row -> "RX" in String.(row[:pins]), r[:results])
        # R-RX names a reaction node; R-E resolves to one
        @test Int(r[:counts][:reaction_ids_excluded]) == 2
        @test Int(r[:counts][:candidates_total]) == 2   # A and D, as without RX
        @test SRV._is_reaction_node(rxnet, "RX")
        @test !SRV._is_reaction_node(rxnet, "A")
        # without a node table the entity type decides
        @test SRV._is_reaction_node(SRV.ReactionNetwork(
            Dict("q" => SRV.NetworkNode("q", "R-q", "Reaction", nothing, "q", BL)),
            SRV.LogicNetworkEdge[], Dict{String, SRV.SetExpansionMapping}()), "q")
        @test !SRV._is_reaction_node(network, "A")
    end

    @testset "C3: an inner failure keeps its status class" begin
        fail(status, msg) = req -> SRV.HTTP.Response(status, SRV.JSON_HEADERS,
                                                     JSON3.write(Dict("status" => "error", "message" => msg)))
        # inner 5xx -> a server fault, reported as 500, never "Invalid input."
        err = try
            SRV._internal_post(fail(500, "Internal server error."), "/api/solve", Dict()); nothing
        catch e
            e
        end
        @test err isa SRV.InnerServerError
        @test SRV.user_facing_error(err) == (500, "Internal server error.")
        # inner 4xx with an allow-listed message -> that message, 400
        err = try
            SRV._internal_post(fail(400, SRV.READOUTS_CONTEXT_ERROR), "/api/readouts", Dict()); nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test SRV.user_facing_error(err) == (400, SRV.READOUTS_CONTEXT_ERROR)
        # end to end: a context bundle whose own-nodes table is unreadable
        # refuses readouts; /api/inverse forwards that specific refusal
        nid = "inverse-context-broken"
        SRV.cache_network!(nid, network)
        SRV.cache_pin_index!(nid, SRV.PinIndex(Dict("R-C" => ["C"], "R-A" => ["A"]), Dict{String, String}()))
        SRV.cache_readout_index!(nid, SRV.ReadoutIndex(
            Dict{String, Vector{String}}(), Set{String}(), false, false,
            String["$(SRV.READOUT_CONTEXT_TABLE): unreadable"], nothing))
        resp = SRV.inverse_handler(post("/api/inverse", JSON3.write(Dict("network_id" => nid,
                                                                        "target" => target_c))))
        @test resp.status == 400
        @test String(jread(resp)[:message]) == SRV.READOUTS_CONTEXT_ERROR
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
