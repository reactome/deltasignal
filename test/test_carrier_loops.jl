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

@testset "default is closed_twin; off is byte-identical to a network with no carrier table" begin
    with_env("DS_CARRIER_LOOP_MODE" => nothing) do
        @test DS.carrier_loop_mode() == "closed_twin"
    end
    for p in (2.0, 0.5)
        a = solve(fixture(closed_table()); mode = "off", p = p)
        b = solve(fixture(nothing); mode = "off", p = p)
        @test a.node_activities == b.node_activities
        @test a.diagnostics["carrier_loop_rule"] == "off"
        @test a.diagnostics["carrier_loop_forms"] == 0
        @test a.diagnostics["carrier_loop_edges_dropped"] == 0
        # the default with no table solves as off and says so
        c = solve(fixture(nothing); mode = nothing, p = p)
        @test c.node_activities == b.node_activities
        @test c.diagnostics["carrier_loop_rule"] == "closed_twin: no carrier table"
        # the default with a table is closed_twin
        d = solve(fixture(closed_table()); mode = nothing, p = p)
        e = solve(fixture(closed_table()); mode = "closed_twin", p = p)
        @test d.node_activities == e.node_activities
        @test d.diagnostics["carrier_loop_rule"] == "closed_twin"
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
    # and before it settles, the value is the sweep count's (relaxation alone:
    # since specs/066 a component stopped above tolerance is re-solved to its
    # fixed point by default, which this check turns off)
    early = with_env("DS_NONCONV_FALLBACK" => "off") do
        solve(fixture(nothing); mode = "off", p = 0.5, iters = 10)
    end
    @test fold(early, "E") > 100 * fold(down, "E")
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
    for mode in ("closed", "closed_supply", "closed_twin", "all")
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

@testset "only the listed activating edge is dropped; a negative edge between the pair stays" begin
    f = fixture(closed_table())
    neg = DS.LogicNetworkEdge("EP", "E", false, false, 1.0, "depletion")
    net = DS.ReactionNetwork(f.nodes, vcat(f.edges, [neg]), f.set_mappings, f.cofactor_stids, f.containment,
                             nothing, nothing, nothing, closed_table())
    out, rule, forms, n = DS.apply_carrier_loops(net, "closed")
    @test rule == "closed" && forms == 1 && n == 1
    @test any(e -> e.parent_uuid == "EP" && e.child_uuid == "E" && !e.is_positive, out.edges)
    @test !any(e -> e.parent_uuid == "EP" && e.child_uuid == "E" && e.is_positive, out.edges)
    # a listed pair with no such edge drops nothing and says so
    out, rule, forms, n = DS.apply_carrier_loops(fixture(DS.CarrierLoopEdge[("E", "T", true)]), "closed")
    @test (forms, n) == (0, 0)
end

@testset "closed_supply: a depletion into a closed free form is dropped too (amendment 4)" begin
    # D, a complex holding E elsewhere, depletes E; D's knockout would de-repress
    # E under closed (the GRB2 / PTPN11 case of P10)
    f = fixture(closed_table())
    dep = DS.LogicNetworkEdge("D", "E", false, false, 1.0, "depletion")
    inh = DS.LogicNetworkEdge("I", "E", false, false, 1.0, "regulator")
    mk(table) = DS.ReactionNetwork(merge(f.nodes, Dict(node.(["D", "I"]))), vcat(f.edges, [dep, inh]),
                                   f.set_mappings, f.cofactor_stids, f.containment, nothing, nothing, nothing, table)
    out, rule, forms, n = DS.apply_carrier_loops(mk(closed_table()), "closed_supply")
    @test rule == "closed_supply" && forms == 1 && n == 2
    @test !any(e -> e.parent_uuid == "D", out.edges)
    @test any(e -> e.parent_uuid == "I", out.edges)          # an inhibition is not binding; it stays
    # closed keeps the depletion; a fed loop keeps it under closed_supply
    out, _, _, n = DS.apply_carrier_loops(mk(closed_table()), "closed")
    @test n == 1 && any(e -> e.parent_uuid == "D", out.edges)
    out, _, forms, n = DS.apply_carrier_loops(mk(fed_table()), "closed_supply")
    @test (forms, n) == (0, 0) && any(e -> e.parent_uuid == "D", out.edges)
    # solved: D knocked out leaves E at baseline under closed_supply, de-represses it under closed
    for (mode, expect_base) in (("closed_supply", true), ("closed", false))
        r = solve(mk(closed_table()); mode = mode, p = 2.0,
                  extra = Dict("D" => (0.0, 1.0), "I" => (1.0, 1.0)))
        @test (fold(r, "E") ≈ 1.0) == expect_base
        expect_base && @test fold(r, "EP") ≈ 2.0 rtol = 1e-6
        @test r.diagnostics["carrier_loop_rule"] == mode
    end
end

@testset "closed_twin: a closed free form follows its gene's root copy (amendment 5)" begin
    # G is the gene's root copy of E (same reactome_id, no in-edge); G -> X elsewhere.
    # H is a copy with an in-edge (not a root), so it is not a supply.
    f = fixture(closed_table())
    twin(u) = u => DS.NetworkNode(u, "R-E", "protein", nothing, u, BL)
    mk(table) = DS.ReactionNetwork(merge(f.nodes, Dict(twin("G"), twin("H"), node("X"))),
                                   vcat(f.edges, [edge("G", "X"), edge("X", "H")]),
                                   f.set_mappings, f.cofactor_stids, f.containment, nothing, nothing, nothing, table)
    out, rule, forms, n, twins = DS.apply_carrier_loops(mk(closed_table()), "closed_twin")
    @test rule == "closed_twin" && forms == 1 && n == 1 && twins == 1
    @test count(e -> e.child_uuid == "E" && e.edge_type == "carrier_supply" && e.parent_uuid == "G" &&
                     e.is_positive && !e.is_and, out.edges) == 1
    @test !any(e -> e.parent_uuid == "H" && e.child_uuid == "E", out.edges)
    for (g, e_fold) in ((0.0, 0.0), (3.0, 3.0))
        r = solve(mk(closed_table()); mode = "closed_twin", p = 2.0, extra = Dict("G" => (g, 1.0)))
        @test fold(r, "E") ≈ e_fold atol = 1e-9
        @test fold(r, "EP") ≈ 2.0 * e_fold atol = 1e-6
        @test r.diagnostics["carrier_loop_twin_edges"] == 1
        # closed_supply ignores the gene: E stays at baseline
        r = solve(mk(closed_table()); mode = "closed_supply", p = 2.0, extra = Dict("G" => (g, 1.0)))
        @test fold(r, "E") ≈ 1.0 rtol = 1e-9
        @test r.diagnostics["carrier_loop_twin_edges"] == 0
    end
    # unperturbed, the twin is baseline and E reads baseline as under closed_supply
    r = solve(mk(closed_table()); mode = "closed_twin", p = 2.0)
    @test fold(r, "E") ≈ 1.0 rtol = 1e-9
    # a fed loop gets no twin edge; a listed free form is never a twin
    @test DS.apply_carrier_loops(mk(fed_table()), "closed_twin")[5] == 0
    both = DS.CarrierLoopEdge[("E", "EP", true), ("G", "X", true)]
    @test DS.apply_carrier_loops(mk(both), "closed_twin")[5] == 0
end

@testset "influence scores use the network the solve ran (review of amendment 6)" begin
    f = fixture(closed_table())
    twin(u) = u => DS.NetworkNode(u, "R-E", "protein", nothing, u, BL)
    net = DS.ReactionNetwork(merge(f.nodes, Dict(twin("G"), node("X"))), vcat(f.edges, [edge("G", "X")]),
                             f.set_mappings, f.cofactor_stids, f.containment, nothing, nothing, nothing,
                             closed_table())
    with_env("DS_CARRIER_LOOP_MODE" => nothing) do
        r = solve(net; mode = nothing, p = 2.0, extra = Dict("G" => (5.0, 1.0)))
        @test r.diagnostics["carrier_loop_twin_edges"] == 1
        score(n) = DS.compute_influence_scores(r, DS.convert_to_reaction_network(n); network = n)
        raw, solved = score(net), score(DS.solved_carrier_network(net))
        # G drives E through the twin edge, and EP no longer feeds E
        @test solved["G"] ≈ raw["G"] + 1.0 atol = 1e-6
        @test solved["EP"] ≈ raw["EP"] - 1.0 atol = 1e-6
        # off: the solved network is the network as given
        with_env("DS_CARRIER_LOOP_MODE" => "off") do
            @test DS.solved_carrier_network(net) === net
        end
    end
end

@testset "an empty reactome_id names no entity, so it makes no twin" begin
    blank(u) = u => DS.NetworkNode(u, "", "protein", nothing, u, BL)
    f = fixture(closed_table())
    nodes = merge(f.nodes, Dict(blank("E"), blank("Z")))
    net = DS.ReactionNetwork(nodes, vcat(f.edges, [edge("Z", "T")]), f.set_mappings, f.cofactor_stids,
                             f.containment, nothing, nothing, nothing, closed_table())
    @test DS.apply_carrier_loops(net, "closed_twin")[5] == 0
end

@testset "a typo in the mode is an error" begin
    for bad in ("Closed", "closed ", "on", "inert", "")
        with_env("DS_CARRIER_LOOP_MODE" => bad) do
            @test_throws ArgumentError DS.carrier_loop_mode()
        end
    end
end

@testset "open loops: a fed free form stops reading its release; a scaffold keeps its steps (specs/077)" begin
    # R is made by r0 (from S) and released by r2; r1 binds R and A into RA,
    # which r3 turns into the readout T. The generator lists both R (released by
    # r2) and RA (made by r1) as free forms of one open loop, as it does for the
    # IFN alpha/beta receptor. Producers into an entity are OR edges.
    out(s, t) = DS.LogicNetworkEdge(s, t, false, true, 1.0, "output")
    ids = ["S", "r0", "R", "A", "r1", "RA", "r2", "r3", "T"]
    edges = [edge("S", "r0"), out("r0", "R"), edge("R", "r1"), edge("A", "r1"), out("r1", "RA"),
             edge("RA", "r2"), out("r2", "R"), edge("RA", "r3"), out("r3", "T")]
    tbl = DS.CarrierLoopEdge[("R", "r2", false), ("RA", "r1", false)]
    net = DS.ReactionNetwork(Dict(node.(ids)), edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                             Dict{String, Set{String}}(), nothing, nothing, nothing, tbl)
    run(mode, open; a = 2.0) = with_env("DS_CARRIER_LOOP_MODE" => mode, "DS_CARRIER_LOOP_OPEN" => open) do
        DS.solve_steady_state(net, Dict("S" => (3.0, 1.0), "A" => (a, 1.0)),
                              DS.SteadyStateParams(1.0, 0.1, 500, 1e-6, "penalty"))
    end
    with_env("DS_CARRIER_LOOP_OPEN" => nothing) do
        @test DS.carrier_loop_open() == "fed_single"
    end
    # fed_single is the default since specs/077's adoption, bit for bit; off reports itself
    d = run("closed_twin", nothing)
    o = run("closed_twin", "off")
    @test d.node_activities == run("closed_twin", "fed_single").node_activities
    @test d.diagnostics["carrier_loop_open_rule"] == "fed_single"
    @test o.diagnostics["carrier_loop_open_rule"] == "off"
    @test o.diagnostics["carrier_loop_open_edges"] == 0
    # fed: R reads S only; RA keeps r1, its one producer, so A still reaches T
    f = run("closed_twin", "fed")
    @test f.converged
    @test fold(f, "R") ≈ 3.0 rtol = 1e-6
    @test fold(f, "RA") ≈ 6.0 rtol = 1e-6
    @test fold(f, "T") ≈ 6.0 rtol = 1e-6
    @test fold(run("closed_twin", "fed"; a = 0.5), "T") ≈ 1.5 rtol = 1e-6
    @test f.diagnostics["carrier_loop_open_rule"] == "fed"
    @test f.diagnostics["carrier_loop_open_edges"] == 1
    @test f.diagnostics["carrier_loop_edges_dropped"] == 1
    # all cuts the scaffold: RA has no producer left and T no longer sees A
    a = run("all", "off")
    @test fold(a, "T") ≈ 1.0 rtol = 1e-9
    @test fold(run("all", "off"; a = 0.5), "T") ≈ 1.0 rtol = 1e-9
    # under all and off, fed adds nothing
    @test run("all", "fed").node_activities == a.node_activities
    @test run("off", "fed").node_activities == run("off", "off").node_activities
    @test run("off", "fed").diagnostics["carrier_loop_open_edges"] == 0
    # the influence network is the one the solve ran
    with_env("DS_CARRIER_LOOP_OPEN" => "fed") do
        sn = DS.solved_carrier_network(net)
        @test !any(e -> e.parent_uuid == "r2" && e.child_uuid == "R", sn.edges)
        @test any(e -> e.parent_uuid == "r1" && e.child_uuid == "RA", sn.edges)
    end
    # fed_single (amendment 1): a release with a second non-cofactor input is a
    # displacement and stays; a cofactor beside it does not make it one
    @test fold(run("closed_twin", "fed_single"), "R") ≈ 3.0 rtol = 1e-6
    @test run("closed_twin", "fed_single").diagnostics["carrier_loop_open_edges"] == 1
    for (extra, rid, kept) in (("D", "R-D", true), ("ATP", "R-ALL-113592", false))
        ids2 = [ids; extra]
        nodes2 = Dict(node.(ids2))
        nodes2[extra] = DS.NetworkNode(extra, rid, "protein", nothing, extra, BL)
        net2 = DS.ReactionNetwork(nodes2, [edges; edge(extra, "r2")], Dict{String, DS.SetExpansionMapping}(),
                                  Set{String}(), Dict{String, Set{String}}(), nothing, nothing, nothing, tbl)
        for open in ("fed", "fed_single")
            n2, _, _, _, _, nopen = DS.apply_carrier_loops(net2, "closed_twin"; open = open)
            has = any(e -> e.parent_uuid == "r2" && e.child_uuid == "R", n2.edges)
            @test has == (open == "fed_single" && kept)
            @test nopen == (has ? 0 : 1)
        end
        r2 = with_env("DS_CARRIER_LOOP_MODE" => "closed_twin", "DS_CARRIER_LOOP_OPEN" => "fed_single") do
            DS.solve_steady_state(net2, Dict("S" => (3.0, 1.0), "A" => (2.0, 1.0)),
                                  DS.SteadyStateParams(1.0, 0.1, 500, 1e-6, "penalty"))
        end
        @test r2.diagnostics["carrier_loop_open_rule"] == "fed_single"
        @test r2.diagnostics["carrier_loop_open_edges"] == (kept ? 0 : 1)
    end
    for bad in ("Fed", "on", "all", "", "single")
        with_env("DS_CARRIER_LOOP_OPEN" => bad) do
            @test_throws ArgumentError DS.carrier_loop_open()
        end
    end
end

@testset "a closed complex with no starting copy follows its components' starting copies (specs/084)" begin
    # E is a complex the pathway only regenerates (E + P -> EP -> E + readout T); it
    # contains K1 (which contains S) and K2, each with a starting copy here.
    cont = Dict("R-E" => Set(["R-K1", "R-K2", "R-S"]), "R-K1" => Set(["R-S"]))
    function cnet(; twin = false)
        ids = ["E", "P", "EP", "T", "K1", "K2", "S"]
        twin && push!(ids, "E2")
        nodes = Dict(node.(ids))
        twin && (nodes["E2"] = DS.NetworkNode("E2", "R-E", "protein", nothing, "E2", BL))
        edges = [edge("E", "EP"), edge("P", "EP"), edge("EP", "E"), edge("EP", "T")]
        DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                           cont, nothing, nothing, nothing, closed_table())
    end
    run(net, flag; mode = "closed_twin", k1 = 1.0, s = 1.0) =
        with_env("DS_CARRIER_LOOP_MODE" => mode, "DS_CARRIER_COMPONENT_TWIN" => flag) do
            obs = Dict("P" => (2.0, 1.0), "K1" => (k1, 1.0), "K2" => (1.0, 1.0), "S" => (s, 1.0))
            DS.solve_steady_state(net, obs, DS.SteadyStateParams(1.0, 0.1, 500, 1e-6, "penalty"))
        end
    with_env("DS_CARRIER_COMPONENT_TWIN" => nothing) do
        @test !DS.carrier_component_twin()
    end
    off = run(cnet(), nothing; k1 = 0.0)
    @test fold(off, "E") ≈ 1.0 rtol = 1e-6               # closed and unfed: baseline, K1 unseen
    @test fold(off, "T") ≈ 2.0 rtol = 1e-6
    @test off.diagnostics["carrier_component_rule"] == "off"
    @test off.diagnostics["carrier_component_edges"] == 0
    ko = run(cnet(), "on"; k1 = 0.0)
    @test fold(ko, "E") ≈ 0.0 atol = 1e-9                # a component knocked out: no complex
    @test fold(ko, "T") ≈ 0.0 atol = 1e-9
    @test ko.diagnostics["carrier_component_rule"] == "on"
    @test ko.diagnostics["carrier_component_edges"] == 2  # K1 and K2; S is inside K1
    up = run(cnet(), "on"; k1 = 2.0)
    @test fold(up, "E") ≈ 2.0 rtol = 1e-6
    @test fold(up, "T") ≈ 4.0 rtol = 1e-6
    @test fold(run(cnet(), "on"; s = 3.0), "E") ≈ 1.0 rtol = 1e-6   # a non-maximal component is not read
    # a starting copy of its own entity: closed_twin feeds it, no component edges
    tw = run(cnet(; twin = true), "on"; k1 = 0.0)
    @test tw.diagnostics["carrier_component_edges"] == 0
    @test tw.diagnostics["carrier_loop_twin_edges"] == 1
    # amendment 1: a form something else still feeds (here a bridge B -> E) is left alone
    bnet = let n = cnet()
        nodes = copy(n.nodes); nodes["B"] = DS.NetworkNode("B", "R-B", "protein", nothing, "B", BL)
        DS.ReactionNetwork(nodes, vcat(n.edges, [DS.LogicNetworkEdge("B", "E", false, true, 1.0, "specialization")]),
                           n.set_mappings, n.cofactor_stids, n.containment, nothing, nothing, nothing, n.carrier_loops)
    end
    br = with_env("DS_CARRIER_LOOP_MODE" => "closed_twin", "DS_CARRIER_COMPONENT_TWIN" => "on") do
        DS.solve_steady_state(bnet, Dict("P" => (2.0, 1.0), "K1" => (0.0, 1.0), "K2" => (1.0, 1.0),
                                          "S" => (1.0, 1.0), "B" => (0.25, 1.0)),
                              DS.SteadyStateParams(1.0, 0.1, 500, 1e-6, "penalty"))
    end
    @test br.diagnostics["carrier_component_edges"] == 0
    @test fold(br, "E") ≈ 0.25 rtol = 1e-6               # it keeps reading its bridge
    # only under closed_twin
    cl = run(cnet(), "on"; mode = "closed", k1 = 0.0)
    @test cl.diagnostics["carrier_component_rule"] == "on: inert under DS_CARRIER_LOOP_MODE=closed"
    @test cl.diagnostics["carrier_component_edges"] == 0
    for bad in ("On", "1", "yes")
        with_env("DS_CARRIER_COMPONENT_TWIN" => bad) do
            @test_throws ArgumentError DS.carrier_component_twin()
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
        write(joinpath(dir, "carrier_loops.csv"),
              "loop_id,r,closed,free_uuid,free_stid,loop_source_uuid\nloop1,1,True,u-e,R-HSA-1,\n")
        @test_throws ArgumentError DS.parse_carrier_loops(ln)       # a missing value is an error, not a skip
        write(joinpath(dir, "carrier_loops.csv"), "loop_id,closed,free_uuid\nloop1,True,u-e\n")
        @test_throws ArgumentError DS.parse_carrier_loops(ln)
        write(joinpath(dir, "carrier_loops.csv"),
              "loop_id,r,closed,free_uuid,free_stid,loop_source_uuid\nloop1,1,maybe,u-e,R-HSA-1,u-ep\n")
        @test_throws ArgumentError DS.parse_carrier_loops(ln)
    end
    @test DS.carrier_loops_json(fixture(nothing)) === nothing
    @test DS.carrier_loops_from_json(nothing) === nothing
    @test_throws ArgumentError DS.carrier_loops_from_json("E")
    @test_throws ArgumentError DS.carrier_loops_from_json(DS.JSON3.read("[1]"))
    @test_throws ArgumentError DS.carrier_loops_from_json(
        DS.JSON3.read("[{\"free_uuid\": 1, \"loop_source_uuid\": \"EP\", \"closed\": true}]"))
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
