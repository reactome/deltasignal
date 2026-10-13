# specs/052: zero supply (DS_ZERO_SUPPLY_MODE). A receptor E binds X into C1;
# C1 cycles through C2 (with co-input Y) and back to C1, and C2 signals to T.
# With E knocked out and Y raised (the de-repressed partners of the EGFR CBL
# cycle), the iterated cycle keeps making complexes of E. Under `hold`, every
# node built from E that no unpinned source can supply is held at 0.
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
inp(s, t) = DS.LogicNetworkEdge(s, t, true, true, 1.0, "input")
outp(s, t) = DS.LogicNetworkEdge(s, t, false, true, 1.0, "output")
cat(s, t) = DS.LogicNetworkEdge(s, t, true, true, 1.0, "catalyst")

const KINDS = Dict("r1" => "reaction", "r2" => "reaction", "r3" => "reaction", "r4" => "reaction",
                   "rs" => "reaction", "rk" => "reaction")
const LEAVES = Dict("E" => ["LE"], "X" => ["LX"], "Y" => ["LY"], "C1" => ["LE", "LX"],
                    "C2" => ["LE", "LX", "LY"], "T" => ["LT"], "M" => ["LM"], "E2" => ["LE"],
                    "K" => ["LE"], "Z" => ["LZ"], "P" => ["LE", "LQ"])

function fixture(; table = true, synthesis = false, catalysed = false, pool = false)
    ids = ["E", "X", "Y", "C1", "C2", "T", "r1", "r2", "r3", "r4"]
    edges = [inp("E", "r1"), inp("X", "r1"), outp("r1", "C1"),
             inp("C1", "r2"), inp("Y", "r2"), outp("r2", "C2"),
             inp("C2", "r3"), outp("r3", "C1"),
             inp("C2", "r4"), outp("r4", "T")]
    if synthesis      # E is also made by a reaction from M: a source the knockout does not pin
        append!(ids, ["M", "rs", "E2"])
        append!(edges, [inp("M", "rs"), outp("rs", "E2"), inp("E2", "r1")])
    end
    if catalysed      # K is built from E but only catalyses Z's production
        append!(ids, ["K", "rk", "Z"])
        append!(edges, [cat("K", "rk"), inp("X", "rk"), outp("rk", "Z")])
    end
    if pool           # an OR pool listing E among its members' leaves
        push!(ids, "P")
        push!(edges, DS.LogicNetworkEdge("C1", "P", false, true, 1.0, "variant_pool"))
    end
    kinds = Dict(u => get(KINDS, u, u == "P" ? "variant_pool" : u in ("C1", "C2") ? "simple_complex" : "simple_entity")
                 for u in ids)
    nl = table ? Dict{String, DS.NodeLeaves}(u => (kinds[u], get(LEAVES, u, String[])) for u in ids) : nothing
    DS.ReactionNetwork(Dict(node.(ids)), edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                       Dict{String, Set{String}}(), nothing, nothing, nothing, nothing, nl)
end

function solve(net; mode, obs = Dict("E" => (0.0, 1.0), "Y" => (10.0, 1.0)))
    with_env("DS_ZERO_SUPPLY_MODE" => mode) do
        DS.solve_steady_state(net, obs, DS.SteadyStateParams(1.0, 0.1, 500, 1e-6, "penalty"))
    end
end
fold(r, u) = r.node_activities[u] / BL

@testset "zero supply (specs/052)" begin

@testset "default is off; off ignores the node table; a typo is an error" begin
    with_env("DS_ZERO_SUPPLY_MODE" => nothing) do
        @test DS.zero_supply_mode() == "off"
    end
    with_env("DS_ZERO_SUPPLY_MODE" => "Hold") do
        @test_throws ArgumentError DS.zero_supply_mode()
    end
    a = solve(fixture(); mode = "off")
    b = solve(fixture(table = false); mode = "off")
    @test a.node_activities == b.node_activities
    @test a.diagnostics["zero_supply_rule"] == "off"
    @test a.diagnostics["zero_supply_held"] == 0
    # the 9-argument constructor leaves the table absent
    @test DS.ReactionNetwork(Dict(node.(["A"])), DS.LogicNetworkEdge[], Dict{String, DS.SetExpansionMapping}(),
                             Set{String}(), Dict{String, Set{String}}(), nothing, nothing, nothing, nothing).node_leaves === nothing
end

@testset "hold without a node table says so and solves as off" begin
    a = solve(fixture(table = false); mode = "hold")
    b = solve(fixture(table = false); mode = "off")
    @test a.node_activities == b.node_activities
    @test a.diagnostics["zero_supply_rule"] == "hold: no node table"
    @test a.diagnostics["zero_supply_held"] == 0
end

@testset "the defect: off keeps complexes of a knocked-out protein alive" begin
    off = solve(fixture(); mode = "off")
    @test fold(off, "C1") > 1.0          # made of E, in a cell with no E
    @test fold(off, "T") > 1.0
end

@testset "hold: complexes of E are held at 0, and the readout follows" begin
    r = solve(fixture(); mode = "hold")
    @test r.node_activities["C1"] == 0.0
    @test r.node_activities["C2"] == 0.0
    @test fold(r, "T") < 0.5
    @test r.diagnostics["zero_supply_rule"] == "hold"
    @test r.diagnostics["zero_supply_leaves"] == 1        # LE only; X and Y are supplied
    @test r.diagnostics["zero_supply_held"] == 2          # C1, C2: pinned E and reactions are not held
    @test fold(r, "X") == 1.0 && fold(r, "Y") == 10.0
end

@testset "a protein made elsewhere is still supplied: nothing held" begin
    a = solve(fixture(synthesis = true); mode = "hold")
    b = solve(fixture(synthesis = true); mode = "off")
    @test a.node_activities == b.node_activities
    @test a.diagnostics["zero_supply_held"] == 0
end

@testset "only a pin to 0 removes supply: knockdown, overexpression, no pin" begin
    for obs in (Dict("E" => (0.5, 1.0), "Y" => (10.0, 1.0)),
                Dict("E" => (2.0, 1.0)),
                Dict("Y" => (10.0, 1.0)),
                Dict("E" => (0.0, 0.0), "Y" => (10.0, 1.0)))    # below the confidence gate: not a pin
        a = solve(fixture(); mode = "hold", obs = obs)
        b = solve(fixture(); mode = "off", obs = obs)
        @test a.node_activities == b.node_activities
        @test a.diagnostics["zero_supply_held"] == 0
    end
end

@testset "an observation wins, and a node pinned above 0 supplies what it carries" begin
    obs = Dict("E" => (0.0, 1.0), "Y" => (10.0, 1.0), "C1" => (1.0, 1.0))
    r = solve(fixture(); mode = "hold", obs = obs)
    @test fold(r, "C1") == 1.0
    @test r.node_activities["C2"] > 0.0
    @test r.diagnostics["zero_supply_held"] == 0
end

@testset "catalysts do not carry; pools are not held" begin
    obs = Dict("E" => (0.0, 1.0), "K" => (0.0, 1.0), "Y" => (10.0, 1.0))
    r = solve(fixture(catalysed = true, pool = true); mode = "hold", obs = obs)
    off = solve(fixture(catalysed = true, pool = true); mode = "off", obs = obs)
    @test r.node_activities["Z"] == off.node_activities["Z"]   # K is built from E, Z is not: not held
    held = r.diagnostics["zero_supply_held"]
    @test held == 2                           # C1, C2; never the pool P, never Z
    zs, rule, n = DS.zero_supply_uuids(fixture(catalysed = true, pool = true),
                                      Dict("E" => 0.0, "K" => 0.0), "hold")
    @test zs == Set(["C1", "C2"]) && rule == "hold" && n == 1
end

@testset "nodes.csv: parse, conflicts, JSON round trip" begin
    mktempdir() do d
        write(joinpath(d, "logic_network.csv"), "source_id,target_id\n")
        write(joinpath(d, "nodes.csv"),
              "uuid,node_kind,diagram_entity_id,compartment,member_leaves,source_sets,chosen_members\n" *
              "a,simple_complex,R-1,,R-3|R-2|R-2,,\nb,reaction,R-9,,,,\na,simple_complex,R-1,,R-2|R-3,,\n")
        nl = DS.parse_node_leaves(joinpath(d, "logic_network.csv"))
        @test nl == Dict("a" => ("simple_complex", ["R-2", "R-3"]), "b" => ("reaction", String[]))
        write(joinpath(d, "nodes.csv"), "uuid,node_kind,member_leaves\na,simple_complex,R-1\na,simple_complex,R-2\n")
        @test_throws ArgumentError DS.parse_node_leaves(joinpath(d, "logic_network.csv"))
        rm(joinpath(d, "nodes.csv"))
        @test DS.parse_node_leaves(joinpath(d, "logic_network.csv")) === nothing
    end
    net = fixture()
    js = DS.JSON3.read(DS.JSON3.write(DS.node_leaves_json(net)))
    @test DS.node_leaves_from_json(js) == net.node_leaves
    @test DS.node_leaves_from_json(nothing) === nothing
    @test_throws ArgumentError DS.node_leaves_from_json(Dict("a" => 1))
end

end
