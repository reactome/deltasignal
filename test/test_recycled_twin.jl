# specs/059: recycled twins (DS_RECYCLED_TWIN). One entity X is split into a
# ROOT copy R (no producer) and a RECYCLED copy Y: Y is recruited (rA) into P,
# and P is dephosphorylated (rD) back into Y, so Y's only producer is its own
# downstream. P signals on to the readout T. R feeds only a side reaction (rS
# -> C). Pinning R (the gene's root) never reaches T: the TGF-beta SMAD2 shape.
# Under `on`, Y gains an OR in-edge from R and T follows R.
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
nd(u, rid) = u => DS.NetworkNode(u, rid, "protein", nothing, u, BL)
e(s, t, et = "input") = DS.LogicNetworkEdge(s, t, et != "output", true, 1.0, et)   # outputs are OR, as the generator emits them
dep(s, t) = DS.LogicNetworkEdge(s, t, true, false, 1.0, "depletion")

# extra: more edges; table: carrier_loops; rdep: give R a depletion in-edge from C
function fixture(; extra = DS.LogicNetworkEdge[], table = nothing, rdep = false)
    nodes = Dict(nd("R", "R-X"), nd("Y", "R-X"), nd("P", "R-PX"), nd("T", "R-T"), nd("C", "R-C"),
                 nd("rS", "R-rS"), nd("rA", "R-rA"), nd("rD", "R-rD"), nd("rT", "R-rT"))
    edges = [e("R", "rS"), e("rS", "C", "output"),
             e("Y", "rA"), e("rA", "P", "output"), e("P", "rD"), e("rD", "Y", "output"),
             e("P", "rT"), e("rT", "T", "output")]
    rdep && push!(edges, dep("C", "R"))
    append!(edges, extra)
    DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                       Dict{String, Set{String}}(), nothing, nothing, nothing, table)
end

solve(net, mode, v) = with_env("DS_RECYCLED_TWIN" => mode) do
    DS.solve_steady_state(net, Dict("R" => (v, 1.0)), DS.SteadyStateParams(1.0, 0.1, 2000, 1e-9, "penalty"))
end
fold(r, u) = r.node_activities[u] / BL

@testset "recycled twins (specs/059)" begin

@testset "default is off; bad values are errors" begin
    with_env("DS_RECYCLED_TWIN" => nothing) do
        @test DS.recycled_twin_mode() == "off"
    end
    for bad in ("ON", "1", "true", "closed")
        with_env("DS_RECYCLED_TWIN" => bad) do
            @test_throws ArgumentError DS.recycled_twin_mode()
        end
    end
end

@testset "off leaves the network as given" begin
    net = fixture()
    out, rule, n = DS.apply_recycled_twins(net, "off")
    @test out === net && rule == "off" && n == 0
end

@testset "on: one recycled_supply edge, root copy to recycled copy" begin
    out, rule, n = DS.apply_recycled_twins(fixture(), "on")
    @test rule == "on" && n == 1
    added = [x for x in out.edges if x.edge_type == "recycled_supply"]
    @test length(added) == 1
    @test (added[1].parent_uuid, added[1].child_uuid) == ("R", "Y")
    @test !added[1].is_and && added[1].is_positive
end

@testset "a root copy may carry a depletion in-edge (the SMAD2 root)" begin
    @test DS.apply_recycled_twins(fixture(rdep = true), "on")[3] == 1
end

@testset "no edge when the root copy already reaches the recycled copy positively" begin
    @test DS.apply_recycled_twins(fixture(extra = [e("C", "Y", "output")]), "on")[3] == 0
end

@testset "...but a reach only through depletion does not count" begin
    @test DS.apply_recycled_twins(fixture(extra = [dep("C", "Y")]), "on")[3] == 1
end

@testset "no edge when the copy has a producer outside its own downstream" begin
    net = fixture(extra = [e("S", "Y", "output")])
    net.nodes["S"] = DS.NetworkNode("S", "R-S", "protein", nothing, "S", BL)
    @test DS.apply_recycled_twins(net, "on")[3] == 0
end

@testset "a carrier free form is left to the carrier rule" begin
    @test DS.apply_recycled_twins(fixture(table = DS.CarrierLoopEdge[("Y", "rD", true)]), "on")[3] == 0
end

@testset "solve: off, the readout ignores the root; on, it follows it" begin
    for v in (0.0, 2.0)
        off = solve(fixture(), "off", v)
        @test fold(off, "T") ≈ 1.0 atol = 1e-6
        @test off.diagnostics["recycled_twin_rule"] == "off"
        @test off.diagnostics["recycled_twin_edges"] == 0
        on = solve(fixture(), "on", v)
        @test on.diagnostics["recycled_twin_rule"] == "on"
        @test on.diagnostics["recycled_twin_edges"] == 1
        @test on.converged
        # Y = mean(regeneration, R): the cycle has gain 1, so Y settles at R
        @test fold(on, "Y") ≈ v atol = 1e-3
        @test fold(on, "T") ≈ v atol = 1e-3
    end
end

@testset "an observation on the recycled copy still wins" begin
    r = with_env("DS_RECYCLED_TWIN" => "on") do
        DS.solve_steady_state(fixture(), Dict("R" => (2.0, 1.0), "Y" => (0.5, 1.0)),
                              DS.SteadyStateParams(1.0, 0.1, 2000, 1e-9, "penalty"))
    end
    @test fold(r, "Y") ≈ 0.5 atol = 1e-9
end

@testset "influence scores see the same network" begin
    with_env("DS_RECYCLED_TWIN" => "on") do
        net = DS.solved_carrier_network(fixture())
        @test any(x -> x.edge_type == "recycled_supply", net.edges)
    end
    with_env("DS_RECYCLED_TWIN" => "off") do
        @test !any(x -> x.edge_type == "recycled_supply", DS.solved_carrier_network(fixture()).edges)
    end
end

end
