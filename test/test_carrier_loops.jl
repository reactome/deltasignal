# specs/049: carrier loops (DS_CARRIER_LOOP_MODE). A carrier E binds a partner P
# into EP and is released unchanged from it; EP signals on to T. Iterated, the
# loop E -> EP -> E has gain P: it rails or collapses with the sweep count (the
# 14-3-3 dimer in RAF). Under `closed`, E's in-edge from its own release is
# dropped, E reads baseline (its supply), and EP and T follow P.
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
edge(s, t) = DS.LogicNetworkEdge(s, t, true, true, 1.0, "input")

# fed = true adds S, a producer of E outside the loop (synthesis); the loop is then not closed
function fixture(table; fed = false)
    ids = fed ? ["E", "P", "EP", "T", "S"] : ["E", "P", "EP", "T"]
    edges = [edge("E", "EP"), edge("P", "EP"), edge("EP", "E"), edge("EP", "T")]
    fed && push!(edges, edge("S", "E"))
    DS.ReactionNetwork(Dict(node.(ids)), edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                       Dict{String, Set{String}}(), nothing, nothing, nothing, table)
end
closed_table() = DS.CarrierLoopEdge[("E", "EP", true)]
fed_table() = DS.CarrierLoopEdge[("E", "EP", false)]

function solve(net; mode, p = 2.0, iters = 500, extra = Dict{String, Tuple{Float64, Float64}}())
    with_env("DS_CARRIER_LOOP_MODE" => mode) do
        obs = merge(Dict("P" => (p, 1.0)), extra)
        DS.solve_steady_state(net, obs, DS.SteadyStateParams(1.0, 0.1, iters, 1e-6, "penalty"))
    end
end
fold(r, u) = r.node_activities[u] / BL

@testset "carrier loops (specs/049)" begin

@testset "default is off, byte-identical to a network with no carrier table" begin
    with_env("DS_CARRIER_LOOP_MODE" => nothing) do
        @test DS.carrier_loop_mode() == "off"
    end
    for p in (2.0, 0.5)
        a = solve(fixture(closed_table()); mode = nothing, p = p)
        b = solve(fixture(nothing); mode = nothing, p = p)
        @test a.node_activities == b.node_activities
        @test a.diagnostics["carrier_loop_rule"] == "off"
        @test a.diagnostics["carrier_loop_forms"] == 0
        @test a.diagnostics["carrier_loop_edges_dropped"] == 0
    end
    # the pre-existing 8-argument constructor leaves the table absent
    @test DS.ReactionNetwork(Dict(node.(["A"])), DS.LogicNetworkEdge[], Dict{String, DS.SetExpansionMapping}(),
                             Set{String}(), Dict{String, Set{String}}(), nothing, nothing, nothing).carrier_loops === nothing
end

@testset "off: the iterated loop rails or collapses instead of following P" begin
    # the defect the rule exists for: loop gain P, so 2x rails to 100x and 0.5x
    # drains to ~0, where the closed rule below reads P. (In RAF the 14-3-3 loop
    # drains slowly enough that the sweep cap decides its value, research.md P3.)
    up = solve(fixture(nothing); mode = "off", p = 2.0)
    @test fold(up, "EP") > 50
    down = solve(fixture(nothing); mode = "off", p = 0.5)
    @test fold(down, "EP") < 0.01
    # and before it settles, the value is the sweep count's
    @test fold(solve(fixture(nothing); mode = "off", p = 0.5, iters = 10), "E") > 100 * fold(down, "E")
end

@testset "closed: E reads its supply (baseline); EP and T follow P at any sweep count" begin
    for p in (2.0, 0.5, 0.0, 50.0), iters in (50, 2000)
        r = solve(fixture(closed_table()); mode = "closed", p = p, iters = iters)
        @test fold(r, "E") ≈ 1.0 rtol = 1e-9
        @test fold(r, "EP") ≈ p rtol = 1e-6
        @test fold(r, "T") ≈ p rtol = 1e-6
        @test r.converged
        @test r.diagnostics["carrier_loop_rule"] == "closed"
        @test r.diagnostics["carrier_loop_forms"] == 1
        @test r.diagnostics["carrier_loop_edges_dropped"] == 1
    end
end

@testset "a fed loop: untouched under closed, reads its outside supply under all" begin
    obs = Dict("S" => (3.0, 1.0))
    a = solve(fixture(fed_table(); fed = true); mode = "closed", p = 2.0, extra = obs)
    b = solve(fixture(fed_table(); fed = true); mode = "off", p = 2.0, extra = obs)
    @test a.node_activities == b.node_activities
    @test a.diagnostics["carrier_loop_rule"] == "closed"
    @test a.diagnostics["carrier_loop_forms"] == 0
    r = solve(fixture(fed_table(); fed = true); mode = "all", p = 2.0, extra = obs)
    @test fold(r, "E") ≈ 3.0 rtol = 1e-6            # S only; the release no longer feeds it
    @test fold(r, "EP") ≈ 6.0 rtol = 1e-6
    @test r.diagnostics["carrier_loop_rule"] == "all"
    @test r.diagnostics["carrier_loop_forms"] == 1
    # and `all` covers closed loops too
    r = solve(fixture(closed_table()); mode = "all", p = 2.0)
    @test fold(r, "E") ≈ 1.0 rtol = 1e-9
end

@testset "a missing carrier table is reported, not pretended" begin
    for mode in ("closed", "all")
        r = solve(fixture(nothing); mode = mode)
        @test r.diagnostics["carrier_loop_rule"] == "$mode: no carrier table"
        @test r.diagnostics["carrier_loop_forms"] == 0
    end
    r = solve(fixture(DS.CarrierLoopEdge[]); mode = "closed")
    @test r.diagnostics["carrier_loop_rule"] == "closed"
end

@testset "an observation on the free form wins" begin
    r = solve(fixture(closed_table()); mode = "closed", p = 2.0, extra = Dict("E" => (3.0, 1.0)))
    @test fold(r, "E") ≈ 3.0 rtol = 1e-9
    @test fold(r, "EP") ≈ 6.0 rtol = 1e-6
end

@testset "a typo in the mode is an error" begin
    for bad in ("Closed", "closed ", "on", "inert", "")
        with_env("DS_CARRIER_LOOP_MODE" => bad) do
            @test_throws ArgumentError DS.carrier_loop_mode()
        end
    end
end

@testset "carrier_loops.csv: parse, both parse return paths, silo-bridge rebuild, JSON round trip" begin
    mktempdir() do dir
        write(joinpath(dir, "logic_network.csv"),
              "source_id,target_id,pos_neg,and_or,edge_type,stoichiometry\n" *
              "u-e,u-ep,pos,and,input,1\nu-p,u-ep,pos,and,input,1\nu-ep,u-e,pos,and,input,1\n")
        write(joinpath(dir, "stid_to_uuid_mapping.csv"), "uuid,stable_id\nu-e,R-HSA-1\nu-p,R-HSA-2\nu-ep,R-HSA-3\n")
        ln, um = joinpath(dir, "logic_network.csv"), joinpath(dir, "stid_to_uuid_mapping.csv")
        @test DS.parse_carrier_loops(ln) === nothing
        @test DS.parse_complete_network(ln, um).carrier_loops === nothing
        write(joinpath(dir, "carrier_loops.csv"), "loop_id,r,closed,free_uuid,free_stid,loop_source_uuid\n")
        @test DS.parse_carrier_loops(ln) == DS.CarrierLoopEdge[]
        write(joinpath(dir, "carrier_loops.csv"),
              "loop_id,r,closed,free_uuid,free_stid,loop_source_uuid\n" *
              "loop2,1,False,u-x,R-HSA-9,u-y\nloop1,48889,True,u-e,R-HSA-1,u-ep\nloop1,48889,True,u-e,R-HSA-1,u-ep\n")
        tab = DS.parse_carrier_loops(ln)
        @test tab == [("u-e", "u-ep", true), ("u-x", "u-y", false)]     # sorted, duplicates dropped
        @test DS.parse_complete_network(ln, um).carrier_loops == tab
        # the other return path: a bundle with a cofactor list
        write(joinpath(dir, "cofactors.csv"),
              "stable_id,molecule,chebi_id,name,in_network,reactome_release\nR-HSA-2,ATP,1,x,1,97\n")
        net = DS.parse_complete_network(ln, um)
        @test net.carrier_loops == tab
        @test DS.carrier_loops_from_json(DS.JSON3.read(DS.JSON3.write(DS.carrier_loops_json(net)))) == tab
        write(joinpath(dir, "carrier_loops.csv"), "loop_id,closed,free_uuid\nloop1,True,u-e\n")
        @test_throws ArgumentError DS.parse_carrier_loops(ln)
        write(joinpath(dir, "carrier_loops.csv"),
              "loop_id,r,closed,free_uuid,free_stid,loop_source_uuid\nloop1,1,maybe,u-e,R-HSA-1,u-ep\n")
        @test_throws ArgumentError DS.parse_carrier_loops(ln)
    end
    @test DS.carrier_loops_json(fixture(nothing)) === nothing
    @test DS.carrier_loops_from_json(nothing) === nothing
    @test_throws ArgumentError DS.carrier_loops_from_json("E")
    @test_throws ArgumentError DS.carrier_loops_from_json(DS.JSON3.read("[{\"free_uuid\": \"E\"}]"))
    @test_throws ArgumentError DS.carrier_loops_from_json(
        DS.JSON3.read("[{\"free_uuid\": \"E\", \"loop_source_uuid\": \"EP\", \"closed\": \"yes\"}]"))
    # the silo-bridge rebuild keeps the table
    f = fixture(closed_table())
    nodes = merge(f.nodes, Dict("S1" => DS.NetworkNode("S1", "R-S", "protein", nothing, "S1", BL),
                                "S2" => DS.NetworkNode("S2", "R-S", "protein", nothing, "S2", BL)))
    net = DS.ReactionNetwork(nodes, vcat(f.edges, [edge("EP", "S1"), edge("S2", "T")]), f.set_mappings,
                             f.cofactor_stids, f.containment, nothing, nothing, nothing, closed_table())
    with_env("DS_SILO_BRIDGE_MAX_REACH" => "50") do
        @test !isempty(DS.silo_bridge_edges(net))
        r = solve(net; mode = "closed", p = 2.0)
        @test r.diagnostics["carrier_loop_forms"] == 1
        @test fold(r, "E") ≈ 1.0 rtol = 1e-9
    end
end
end
