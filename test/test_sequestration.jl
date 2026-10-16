# specs/063: sequestration (DS_SEQUESTRATION). A sequesterer A binds partner B
# (reaction rB) into complex C, which only dissociates (sink K): a dead end.
# B also signals: B -> rS -> T (the readout). Pinning A never reaches T unless
# C depletes B: knocking A out frees B (T UP); overexpressing A binds more B
# (T DOWN). The rule was measured under DS_DEPLETION_OWN_PRODUCT=full (the
# default then); C is built from B, so under the specs/061 rule the edge is an
# own product, and the fixtures pin `full` (CI broke when the default moved).
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
const KINDS = Dict("A" => "simple_entity", "B" => "simple_entity", "C" => "simple_complex",
                   "K" => "dissociation_sink", "T" => "simple_entity", "Z" => "simple_entity",
                   "ATP" => "simple_entity", "W" => "simple_entity",
                   "rB" => "reaction", "rS" => "reaction", "rZ" => "reaction", "rD" => "reaction")
e(s, t, et = "input") = DS.LogicNetworkEdge(s, t, et != "output", true, 1.0, et)

# extra: more edges; leaves: give the network its node table or not
function fixture(; extra = DS.LogicNetworkEdge[], leaves = true, cofactors = Set{String}())
    nodes = Dict(u => DS.NetworkNode(u, "R-$u", "protein", nothing, u, BL) for u in keys(KINDS))
    edges = [e("A", "rB"), e("B", "rB"), e("rB", "C", "output"), e("C", "K", "dissociation"),
             e("B", "rS"), e("rS", "T", "output")]
    append!(edges, extra)
    nl = leaves ? Dict{String, DS.NodeLeaves}(u => (k, String[]) for (u, k) in KINDS) : nothing
    DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}(), cofactors,
                       Dict{String, Set{String}}(), nothing, nothing, nothing, nothing, nl)
end
seq_edges(net, mode) = sort([(x.parent_uuid, x.child_uuid) for x in DS.apply_sequestration(net, mode)[1].edges
                             if x.edge_type == "depletion"])
solve(net, mode, v) = with_env("DS_DEPLETION_OWN_PRODUCT" => "full", "DS_SEQUESTRATION" => mode) do
    DS.solve_steady_state(net, Dict("A" => (v, 1.0)), DS.SteadyStateParams(1.0, 0.1, 2000, 1e-9, "penalty"))
end
fold(r, u) = r.node_activities[u] / BL

@testset "sequestration (specs/063)" begin

@testset "default is off; bad values are errors" begin
    with_env("DS_DEPLETION_OWN_PRODUCT" => "full", "DS_SEQUESTRATION" => nothing) do
        @test DS.sequestration_mode() == "off"
    end
    with_env("DS_DEPLETION_OWN_PRODUCT" => "full", "DS_SEQUESTRATION" => "dead_end") do
        @test_throws ArgumentError DS.sequestration_mode()
    end
end

@testset "off leaves the network as given; no node table is reported" begin
    net = fixture()
    @test DS.apply_sequestration(net, "off")[1] === net
    out = DS.apply_sequestration(fixture(leaves = false), "deadend")
    @test out[2] == "deadend: no node table" && out[3] == 0
end

@testset "one depletion edge, complex to the partner with another use" begin
    @test seq_edges(fixture(), "deadend_sole") == [("C", "B")]
    @test seq_edges(fixture(), "deadend") == [("C", "B")]
end

@testset "sole_direct: the sequesterer depletes the partner, not the complex" begin
    @test seq_edges(fixture(), "sole_direct") == [("A", "B")]
    net = fixture(extra = [e("A", "rZ"), e("rZ", "Z", "output")])
    @test isempty(seq_edges(net, "sole_direct"))
end

@testset "sole_direct solve: knockout UP, overexpression DOWN, partner's own fall not fed back" begin
    @test fold(solve(fixture(), "sole_direct", 0.0), "T") > 1.15
    @test fold(solve(fixture(), "sole_direct", 2.0), "T") < 0.85
    # B is made from U (rU); lowering U lowers B. Under sole_direct T follows B
    # as under off; under deadend_sole the complex falls with B and props it up.
    up = [e("U", "rU"), e("rU", "B", "output")]
    pu(mode) = with_env("DS_DEPLETION_OWN_PRODUCT" => "full", "DS_SEQUESTRATION" => mode) do
        net = fixture(extra = up)
        net.nodes["U"] = DS.NetworkNode("U", "R-U", "protein", nothing, "U", BL)
        net.nodes["rU"] = DS.NetworkNode("rU", "R-rU", "protein", nothing, "rU", BL)
        net.node_leaves["U"] = ("simple_entity", String[]); net.node_leaves["rU"] = ("reaction", String[])
        DS.solve_steady_state(net, Dict("U" => (0.5, 1.0)), DS.SteadyStateParams(1.0, 0.1, 2000, 1e-9, "penalty"))
    end
    @test fold(pu("off"), "T") < 0.85
    @test fold(pu("sole_direct"), "T") ≈ fold(pu("off"), "T") atol = 1e-9
    @test abs(fold(pu("deadend_sole"), "T") - fold(pu("off"), "T")) > 0.05
end

@testset "no edge when the complex is used" begin
    net = fixture(extra = [e("C", "rD"), e("rD", "W", "output")])
    @test isempty(seq_edges(net, "deadend"))
end

@testset "deadend_sole needs a sequesterer used nowhere else" begin
    # A also feeds rZ: both inputs are partners
    net = fixture(extra = [e("A", "rZ"), e("rZ", "Z", "output")])
    @test isempty(seq_edges(net, "deadend_sole"))
    @test seq_edges(net, "deadend") == [("C", "A"), ("C", "B")]
end

@testset "a cofactor input is neither counted nor a partner" begin
    # ATP joins the binding and has another use; it must not become a partner
    net = fixture(extra = [e("ATP", "rB"), e("ATP", "rZ"), e("rZ", "Z", "output")],
                  cofactors = Set(["R-ATP"]))
    @test seq_edges(net, "deadend_sole") == [("C", "B")]
    # a binding of one protein with a cofactor is not a binding
    net = fixture(extra = [e("ATP", "rS")], cofactors = Set(["R-ATP"]))
    @test !(("T", "B") in seq_edges(net, "deadend"))
end

@testset "solve: off, A cannot reach T; on, knockout UP and overexpression DOWN" begin
    @test fold(solve(fixture(), "off", 0.0), "T") ≈ 1.0 atol = 1e-6
    @test fold(solve(fixture(), "off", 2.0), "T") ≈ 1.0 atol = 1e-6
    ko = solve(fixture(), "deadend_sole", 0.0)
    oe = solve(fixture(), "deadend_sole", 2.0)
    @test fold(ko, "T") > 1.15
    @test fold(oe, "T") < 0.85
    @test ko.diagnostics["sequestration_rule"] == "deadend_sole"
    @test ko.diagnostics["sequestration_edges"] == 1
    # an observation on the partner still wins
    obs = with_env("DS_DEPLETION_OWN_PRODUCT" => "full", "DS_SEQUESTRATION" => "deadend_sole") do
        DS.solve_steady_state(fixture(), Dict("A" => (0.0, 1.0), "B" => (1.0, 1.0)),
                              DS.SteadyStateParams(1.0, 0.1, 2000, 1e-9, "penalty"))
    end
    @test fold(obs, "B") ≈ 1.0 atol = 1e-9
end

@testset "influence scores see the same network" begin
    with_env("DS_DEPLETION_OWN_PRODUCT" => "full", "DS_SEQUESTRATION" => "deadend_sole") do
        net = DS.solved_carrier_network(fixture())
        @test any(x -> x.edge_type == "depletion" && x.parent_uuid == "C", net.edges)
    end
end

end
