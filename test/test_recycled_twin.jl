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
function fixture(; extra = DS.LogicNetworkEdge[], table = nothing, rdep = false,
                 containment = Dict{String, Set{String}}())
    nodes = Dict(nd("R", "R-X"), nd("Y", "R-X"), nd("P", "R-PX"), nd("T", "R-T"), nd("C", "R-C"),
                 nd("rS", "R-rS"), nd("rA", "R-rA"), nd("rD", "R-rD"), nd("rT", "R-rT"))
    edges = [e("R", "rS"), e("rS", "C", "output"),
             e("Y", "rA"), e("rA", "P", "output"), e("P", "rD"), e("rD", "Y", "output"),
             e("P", "rT"), e("rT", "T", "output")]
    rdep && push!(edges, dep("C", "R"))
    append!(edges, extra)
    DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                       containment, nothing, nothing, nothing, table)
end

solve(net, mode, v) = with_env("DS_RECYCLED_TWIN" => mode) do
    DS.solve_steady_state(net, Dict("R" => (v, 1.0)), DS.SteadyStateParams(1.0, 0.1, 2000, 1e-9, "penalty"))
end
fold(r, u) = r.node_activities[u] / BL

@testset "recycled twins (specs/059)" begin

@testset "default is on (specs/061 adoption); bad values are errors" begin
    with_env("DS_RECYCLED_TWIN" => nothing) do
        @test DS.recycled_twin_mode() == "on"
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

@testset "a carrier free form is never a root copy (specs/061 amendment 2)" begin
    # F: another copy of X with an out-edge and no in-edge (as a closed free form
    # is left after the carrier rule drops its loop in-edges). Listed, it feeds
    # nothing; unlisted, it is a root like R.
    net(table) = begin
        n = fixture(extra = [e("F", "rS")], table = table)
        n.nodes["F"] = DS.NetworkNode("F", "R-X", "protein", nothing, "F", BL)
        n
    end
    added(n) = sort([(x.parent_uuid, x.child_uuid) for x in DS.apply_recycled_twins(n, "on")[1].edges
                     if x.edge_type == "recycled_supply"])
    @test added(net(DS.CarrierLoopEdge[("F", "rS", true)])) == [("R", "Y")]
    @test added(net(nothing)) == [("F", "Y"), ("R", "Y")]
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

# --- specs/060 amendment 1: own product by entity ---------------------------
# C (built from R, so its entity contains X) depletes the recycled copy Y: the
# SMAD2:SMURF2 -| SMAD2 shape. A knockout of R zeroes C.
smad_net() = fixture(extra = [dep("C", "Y")], containment = Dict("R-C" => Set(["R-X"])))

solve2(net, own, v) = with_env("DS_RECYCLED_TWIN" => "on", "DS_DEPLETION_OWN_PRODUCT" => own) do
    DS.solve_steady_state(net, Dict("R" => (v, 1.0)), DS.SteadyStateParams(1.0, 0.1, 2000, 1e-9, "penalty"))
end

@testset "own_entity pairs: only under suppress_only_entity, only by containment" begin
    for (own, n) in (("full", 0), ("suppress_only", 0), ("suppress_only_entity", 1))
        with_env("DS_DEPLETION_OWN_PRODUCT" => own) do
            @test length(DS.own_entity_depletion_pairs(smad_net())) == n
        end
    end
    with_env("DS_DEPLETION_OWN_PRODUCT" => "suppress_only_entity") do
        @test DS.own_entity_depletion_pairs(smad_net()) == Set([("C", "Y")])
        # no containment, no pair
        @test isempty(DS.own_entity_depletion_pairs(fixture(extra = [dep("C", "Y")])))
    end
end

@testset "knockout: full de-represses the recycled copy; suppress_only_entity does not" begin
    full = solve2(smad_net(), "full", 0.0)
    ent = solve2(smad_net(), "suppress_only_entity", 0.0)
    @test fold(full, "T") > 1.15          # the SMAD2/3 knockout inversion
    @test fold(ent, "T") < 0.85           # follows the knockout down
    @test ent.diagnostics["own_entity_depletions"] == 1
    @test full.diagnostics["own_entity_depletions"] == 0
    # node-copy suppress_only misses it (C is not built from Y)
    @test fold(solve2(smad_net(), "suppress_only", 0.0), "T") > 1.15
end

@testset "suppress_only_entity still lets the depleter suppress" begin
    ent = solve2(smad_net(), "suppress_only_entity", 2.0)
    @test fold(ent, "C") > 1.15
    @test fold(ent, "Y") < fold(solve2(fixture(), "suppress_only_entity", 2.0), "Y")
end

# --- specs/061: suppress_only_source ----------------------------------------
# C is built from R and a partner M (AND). A knockout of R (SMAD2 KO) zeroes
# C's source copy: C may not de-repress Y. A knockout of M (MDM2 KO) zeroes C
# while R holds: C's fall de-represses Y as under `full`.
mdm_net() = begin
    net = fixture(extra = [dep("C", "Y"), e("M", "rS")], containment = Dict("R-C" => Set(["R-X"])))
    net.nodes["M"] = DS.NetworkNode("M", "R-M", "protein", nothing, "M", BL)
    net
end
solve3(net, own, pins) = with_env("DS_RECYCLED_TWIN" => "on", "DS_DEPLETION_OWN_PRODUCT" => own) do
    DS.solve_steady_state(net, Dict(k => (v, 1.0) for (k, v) in pins),
                          DS.SteadyStateParams(1.0, 0.1, 2000, 1e-9, "penalty"))
end

@testset "own_source map: only under suppress_only_source; sources walk back from the depleter" begin
    for own in ("full", "suppress_only", "suppress_only_entity")
        with_env("DS_DEPLETION_OWN_PRODUCT" => own) do
            @test isempty(DS.own_source_depletion_sources(mdm_net()))
        end
    end
    with_env("DS_DEPLETION_OWN_PRODUCT" => "suppress_only_source") do
        @test DS.own_source_depletion_sources(mdm_net()) == Dict(("C", "Y") => ["R"])
        @test isempty(DS.own_entity_depletion_pairs(mdm_net()))   # the static cap is off
        # no containment, no pair
        @test isempty(DS.own_source_depletion_sources(fixture(extra = [dep("C", "Y")])))
        # amendment 1: a depleter built from the target itself has no source
        # (R is a root copy; C -| R would be self-sourced)
        net = mdm_net(); push!(net.edges, dep("C", "R"))
        @test DS.own_source_depletion_sources(net) == Dict(("C", "Y") => ["R"])
        # amendment 2: a copy that only catalyses or regulates the depleter's
        # reaction is not a source; one that is an input is
        for (et, want) in (("catalyst", ["R"]), ("regulator", ["R"]), ("input", ["I", "R"]))
            net = mdm_net(); push!(net.edges, e("I", "rS", et))
            net.nodes["I"] = DS.NetworkNode("I", "R-X", "protein", nothing, "I", BL)
            @test DS.own_source_depletion_sources(net) == Dict(("C", "Y") => want)
        end
    end
end

@testset "the target's own source falls (SMAD2 KO): no de-repression" begin
    src = solve3(mdm_net(), "suppress_only_source", ("R" => 0.0,))
    @test fold(src, "C") < 0.01
    @test fold(src, "T") < 0.85
    @test src.diagnostics["own_source_depletions"] == 1
    @test fold(solve3(mdm_net(), "full", ("R" => 0.0,)), "T") > 1.15
end

@testset "the partner falls (MDM2 KO): de-repression kept, as full" begin
    src = solve3(mdm_net(), "suppress_only_source", ("M" => 0.0,))
    full = solve3(mdm_net(), "full", ("M" => 0.0,))
    @test fold(src, "C") < 0.01
    @test fold(src, "Y") > 1.15
    @test fold(src, "T") ≈ fold(full, "T") atol = 1e-9
    # the static entity test forbids it
    @test fold(solve3(mdm_net(), "suppress_only_entity", ("M" => 0.0,)), "Y") ≈ 1.0 atol = 1e-3
end

@testset "suppress_only_source still lets the depleter suppress" begin
    src = solve3(mdm_net(), "suppress_only_source", ("R" => 2.0,))
    full = solve3(mdm_net(), "full", ("R" => 2.0,))
    @test fold(src, "C") > 1.15
    @test fold(src, "T") ≈ fold(full, "T") atol = 1e-9
end

end
