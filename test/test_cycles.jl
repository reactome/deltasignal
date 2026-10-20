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

function fixture(; pools = :default, rl = identity, extra = [])
    ids = ["G", "P", "A", "B", "EF", "EB", "F", "R", "W", "D", "Z", "Y", "F2", "K", "Q", "M", "So", "Sr", "Sq"]
    nodes = Dict(rl(k) => DS.NetworkNode(rl(k), "R-" * k, "protein", nothing, rl(k), BL) for k in ids)
    E = [("G", "P", "input", true), ("P", "A", "output", false),
         ("A", "F", "input", true), ("EF", "F", "catalyst", true), ("F", "B", "output", false),
         ("B", "R", "input", true), ("EB", "R", "catalyst", true), ("R", "A", "output", false),
         ("B", "W", "input", true), ("W", "D", "output", false)]
    edges = [edge(rl(s), rl(t), et; and = a) for (s, t, et, a) in E]
    for (s, t, et, pos) in extra
        push!(edges, DS.LogicNetworkEdge(rl(s), rl(t), true, pos, 1.0, et))
    end
    p = pools === :default ?
        [DS.CyclePool("pool1", sort([rl("A"), rl("B")]), rl("A"),
                      [(rl("A"), rl("B"), rl("F")), (rl("B"), rl("A"), rl("R"))])] : pools
    DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                       Dict{String, Set{String}}(), nothing, p)
end

function solve(net; mode = "balance", phi = "0.1", obs = Dict{String, Float64}(), rl = identity,
               method = nothing)
    with_env("DS_CYCLE_MODE" => mode, "DS_CYCLE_PHI" => phi, "DS_SELF_INHIBITOR_WEIGHT" => "off",
             "DS_SCC_METHOD" => method) do
        o = Dict(rl(k) => (1.0, 1.0) for k in ("G", "EF", "EB", "Z"))
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

@testset "default balance, and a missing table is reported" begin
    with_env("DS_CYCLE_MODE" => nothing) do
        @test DS.cycle_mode() == "balance"
    end
    r = solve(fixture(); mode = "off")
    @test r.diagnostics["cycle_rule"] == "off"
    @test r.diagnostics["cycle_pools_solved"] == 0
    r = solve(fixture(; pools = nothing))
    @test r.diagnostics["cycle_rule"] == "balance: no pool table"
    @test r.diagnostics["cycle_pools_solved"] == 0
    # only the SCC fixed point solves pools; any other method must not say "balance"
    for m in ("minimize", "pool", "pool_parity", "pool_all")
        r = solve(fixture(); method = m)
        @test r.diagnostics["cycle_rule"] == "balance: inert under DS_SCC_METHOD=$m"
        @test r.diagnostics["cycle_pools_solved"] == 0
    end
    r = solve(fixture(); method = "fixed_point")
    @test r.diagnostics["cycle_rule"] == "balance"
    @test r.diagnostics["cycle_pools_solved"] == 1
end

@testset "a pinned form sets the pool" begin
    r = solve(fixture(); obs = Dict("B" => 5.0))
    @test fold(r, "B") ≈ 5.0 atol = 1e-9
    @test fold(r, "A") ≈ 5.0 atol = 1e-4              # unperturbed split: supply 5x
end

@testset "review fixes (specs/039 amendment 1)" begin
    # finding 5: every drive at zero leaves the split at pi0, not wherever a
    # constant added to Q puts it
    r = solve(fixture(); obs = Dict("EF" => 0.0, "EB" => 0.0))
    @test r.converged
    @test fold(r, "A") ≈ 1.0 atol = 1e-6
    @test fold(r, "B") ≈ 1.0 atol = 1e-6
    # finding 6: a producer of a form inside the pool's own component is not
    # supply. Y: B -> A is such a producer; the supply stays P's.
    net = fixture(; extra = [("B", "Y", "input", true), ("Y", "A", "output", true)])
    r = solve(net; obs = Dict("G" => 80.0))
    @test r.converged
    @test fold(r, "A") ≈ 80.0 atol = 1e-4
    @test fold(r, "B") ≈ 80.0 atol = 1e-4
    # finding 4: a depletion edge into a form from outside the pool still acts
    net = fixture(; extra = [("Z", "B", "depletion", false)])
    r0 = solve(net)
    @test fold(r0, "B") ≈ 1.0 atol = 1e-6               # Z at baseline: no effect
    r = solve(net; obs = Dict("Z" => 10.0))
    @test r.converged
    @test fold(r, "B") < 0.99
    @test fold(r, "D") ≈ fold(r, "B") atol = 1e-6         # and the readout sees it
    # ... but the pool's own catalyst ⊣ source-form depletion is not applied
    # again (pre-registered: the balance already contains it)
    net = fixture(; extra = [("EF", "A", "depletion", false)])
    r = solve(net; obs = Dict("EF" => 80.0))
    @test fold(r, "A") ≈ 0.11236 atol = 1e-4
    @test fold(r, "B") ≈ 8.98876 atol = 1e-4
    # Fable finding 5: an uncatalysed step beside a catalysed one for the same
    # pair carries a negligible baseline share, so GEF knockout empties the
    # modified form instead of halving it
    xtra = [("A", "F2", "input", true), ("F2", "B", "output", true)]
    pools2(uncat) = [DS.CyclePool("pool1", ["A", "B"], "A",
                                  sort([("A", "B", "F"), ("A", "B", "F2"), ("B", "A", "R")]), uncat)]
    r = solve(fixture(; extra = xtra, pools = pools2(Set(["F2"]))); obs = Dict("EF" => 0.0))
    @test r.converged
    @test fold(r, "B") < 0.01
    r = solve(fixture(; extra = xtra, pools = pools2(Set{String}())); obs = Dict("EF" => 0.0))
    @test 0.4 < fold(r, "B") < 0.7                       # without the flag: the old equal split
    r = solve(fixture(; extra = xtra, pools = pools2(Set(["F2"]))))
    @test fold(r, "B") ≈ 1.0 atol = 1e-6                 # baseline is still exact
    # Fable finding 2: a producer that shares the pool's component only through
    # a pinned cofactor (K, released by R and consumed by M) is still supply
    net = fixture(; extra = [("R", "K", "output", true), ("K", "M", "input", true),
                             ("Q", "M", "input", true), ("M", "A", "output", true)])
    r = solve(net; obs = Dict("K" => 1.0, "Q" => 80.0))
    @test r.converged
    @test fold(r, "A") ≈ 40.5 atol = 1e-3                # mean of P (1x) and M (80x)
    @test fold(r, "B") ≈ 40.5 atol = 1e-3
    # review of amendment 2, finding 1: a depletion edge from the pool's own
    # downstream onto the supply reaction does not make the supply "fed by the
    # pool" -- only mass flow does
    net = fixture(; extra = [("D", "P", "depletion", false)])
    r = solve(net; obs = Dict("G" => 80.0))
    @test r.converged
    @test fold(r, "A") ≈ fold(r, "P") atol = 1e-6
    @test fold(r, "A") > 5.0
    # amendment 3: a derived depletion edge from a pool-fed node onto the
    # supply chain (B -> W; W -| P, where G -> P -> A is the supply) is not
    # applied: the pool does not suppress its own supply
    net = fixture(; extra = [("W", "P", "depletion", false)])
    r = solve(net; obs = Dict("G" => 80.0))
    @test r.converged
    @test fold(r, "A") ≈ 80.0 atol = 1e-3
    @test fold(r, "B") ≈ 80.0 atol = 1e-3
    @test r.diagnostics["cycle_supply_depletions_held"] == 1
    r = solve(net; obs = Dict("EF" => 80.0))
    @test fold(r, "B") ≈ 8.98876 atol = 1e-3           # the fixture, unmoved by W's rise
    # ... and off, the pool table changes nothing: a pool-bearing network
    # solves to exactly the values of the same network with no table
    for obs in (Dict("G" => 80.0), Dict("EF" => 80.0), Dict("EB" => 0.0))
        ra = solve(net; mode = "off", obs = obs)
        rb = solve(fixture(; extra = [("W", "P", "depletion", false)], pools = nothing); mode = "off", obs = obs)
        @test ra.diagnostics["cycle_rule"] == "off"
        @test all(k -> ra.node_activities[k] == rb.node_activities[k], keys(rb.node_activities))
    end
    # amendment 5 (post hoc): an inhibitor built from the pool's own base state
    # (SOCS-bound receptor: A + So -> Sq, Sq -| F) is re-evaluated in the
    # drives with the pool at baseline, so it carries SOCS's fold but not A's:
    # gene 80x moves both forms 80x instead of the inhibitor tracking A and
    # cancelling the forward drive exactly (IFN alpha/beta) ...
    socs = [("A", "Sr", "input", true), ("So", "Sr", "input", true), ("Sr", "Sq", "output", true),
            ("Sq", "F", "regulator", false)]
    net = fixture(; extra = socs)
    r = solve(net; obs = Dict("G" => 80.0, "So" => 1.0))
    @test r.converged
    @test fold(r, "A") ≈ 80.0 atol = 1e-3
    @test fold(r, "B") ≈ 80.0 atol = 1e-3
    @test r.diagnostics["cycle_self_fed_inputs"] >= 1
    r = solve(net; obs = Dict("So" => 1.0))
    @test fold(r, "B") ≈ 1.0 atol = 1e-9                  # baseline exact
    # ... and SOCS itself still acts (amendment 4's hold silenced it): SOCS KO
    # de-represses the forward step, SOCS 80x represses it
    r = solve(net; obs = Dict("So" => 0.0))
    @test r.converged
    @test fold(r, "B") > 1.15
    r = solve(net; obs = Dict("So" => 80.0))
    @test r.converged
    @test fold(r, "B") < 0.5
    @test fold(r, "B") ≈ (1 / 80) / (0.1 * (1 / 80) + 0.9) atol = 1e-3   # u = 1/SOCS fold
    # an OUTSIDE inhibitor of the same step still acts (Z, not fed by the pool)
    net = fixture(; extra = [("Z", "F", "regulator", false)])
    r = solve(net; obs = Dict("Z" => 10.0))
    @test fold(r, "B") < 0.5
    # finding 9: a component method that bypasses the sweep manages no pool,
    # and says so
    for m in ("pool", "minimize")
        r = solve(fixture(); method = m)
        @test r.diagnostics["cycle_pools_solved"] == 0
        @test r.diagnostics["cycle_pools_unmanaged"] == 1
    end
end

# A kinase/phosphatase cycle curated in steps (problem2.md, amendment 2):
#   S + K -> SE -> SsE -> Ss + K        (bind, modify, release)
#   Ss + Ph -> SsP -> SP -> S + Ph
# G -> P -> S supplies the protein; Ss -> W -> D is the readout. K has its own
# supply KG -> KP -> K besides the release step (a carrier).
function ring(; rl = identity, carriers = true, copies = false, c1cat = false, extra_carriers = [])
    ids = ["G", "P", "S", "SE", "SsE", "Ss", "SsP", "SP", "K", "Ph", "b1", "c1", "r1", "b2", "c2", "r2",
           "W", "D", "KG", "KP", "K2", "b1x", "C1"]
    nodes = Dict(rl(k) => DS.NetworkNode(rl(k), "R-" * k, "protein", nothing, rl(k), BL) for k in ids)
    E = [("G", "P", "input"), ("P", "S", "output"),
         ("S", "b1", "input"), ("K", "b1", "input"), ("b1", "SE", "output"),
         ("SE", "c1", "input"), ("c1", "SsE", "output"),
         ("SsE", "r1", "input"), ("r1", "Ss", "output"), ("r1", "K", "output"),
         ("Ss", "b2", "input"), ("Ph", "b2", "input"), ("b2", "SsP", "output"),
         ("SsP", "c2", "input"), ("c2", "SP", "output"),
         ("SP", "r2", "input"), ("r2", "S", "output"), ("r2", "Ph", "output"),
         ("Ss", "W", "input"), ("W", "D", "output"),
         ("KG", "KP", "input"), ("KP", "K", "output")]
    copies && append!(E, [("S", "b1x", "input"), ("K2", "b1x", "input"), ("b1x", "SE", "output")])
    c1cat && push!(E, ("C1", "c1", "input"))
    edges = [edge(rl(s), rl(t), et; and = et == "input") for (s, t, et) in E]
    st(a, b, rs...) = (rl(a), rl(b), sort([rl(r) for r in rs]))
    bind1 = copies ? st("S", "SE", "b1", "b1x") : st("S", "SE", "b1")
    paths = [DS.PoolPath(rl("S"), rl("Ss"), [bind1, st("SE", "SsE", "c1"), st("SsE", "Ss", "r1")], true),
             DS.PoolPath(rl("Ss"), rl("S"), [st("Ss", "SsP", "b2"), st("SsP", "SP", "c2"), st("SP", "S", "r2")], true)]
    cs = carriers ? [(rl("K"), rl("r1")), (rl("Ph"), rl("r2"))] : Tuple{String, String}[]
    append!(cs, [(rl(a), rl(b)) for (a, b) in extra_carriers])
    p = [DS.CyclePool("pool1", sort([rl("S"), rl("Ss")]), rl("S"), paths,
                      sort([rl(k) for k in ("SE", "SsE", "SsP", "SP")]), sort(cs))]
    DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                       Dict{String, Set{String}}(), nothing, p)
end

function solve_ring(net; obs = Dict{String, Float64}(), phi = "0.1", carriers = "on", rl = identity,
                    pinned = ("G", "K", "Ph"))
    with_env("DS_CYCLE_MODE" => "balance", "DS_CYCLE_PHI" => phi, "DS_CYCLE_CARRIERS" => carriers,
             "DS_SELF_INHIBITOR_WEIGHT" => "off") do
        o = Dict(rl(k) => (1.0, 1.0) for k in pinned)
        for (k, v) in obs; o[rl(k)] = (v, 1.0); end
        DS.solve_steady_state(net, o, DS.SteadyStateParams(1.0, 0.1, 500, 1e-9, "penalty"))
    end
end

@testset "a carrier shared by two pools leaves out both pools' release steps (specs/049 amendment 7)" begin
    # two copies of the ring share the kinase K (with its supply KG -> KP -> K):
    # K is a carrier of both pools and must read only KP, not the other pool's
    # release flux (the RAS GAPs in RAF were written from it, and the solve never
    # converged).
    shared = ("K", "KG", "KP")
    ra(k) = k in shared ? k : "a_" * k
    rb(k) = k in shared ? k : "b_" * k
    na, nb = ring(; rl = ra), ring(; rl = rb)
    seen = Set{Tuple{String, String, String}}()
    edges = DS.LogicNetworkEdge[]
    for e in vcat(na.edges, nb.edges)
        key = (e.parent_uuid, e.child_uuid, e.edge_type)
        key in seen && continue
        push!(seen, key); push!(edges, e)
    end
    net = DS.ReactionNetwork(merge(na.nodes, nb.nodes), edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                             Dict{String, Set{String}}(), nothing, vcat(na.pools, nb.pools))
    r = with_env("DS_CYCLE_MODE" => "balance", "DS_CYCLE_PHI" => "0.1", "DS_CYCLE_CARRIERS" => "on",
                 "DS_SELF_INHIBITOR_WEIGHT" => "off") do
        o = Dict(k => (1.0, 1.0) for k in ("a_G", "a_Ph", "b_G", "b_Ph"))
        o["KG"] = (80.0, 1.0)
        DS.solve_steady_state(net, o, DS.SteadyStateParams(1.0, 0.1, 500, 1e-9, "penalty"))
    end
    @test r.converged
    @test fold(r, "K") ≈ 80.0 atol = 1e-6
    @test fold(r, "a_Ss") ≈ fold(r, "b_Ss") atol = 1e-9
    @test r.diagnostics["cycle_carriers_held"] == 1      # one distinct carrier
    # specs/072 amendment 1: the other pool's listed release (the same reaction,
    # R-r1) is already left out, so `reaction` adds nothing and moves nothing
    r2 = with_env("DS_CYCLE_MODE" => "balance", "DS_CYCLE_PHI" => "0.1", "DS_CYCLE_CARRIERS" => "on",
                  "DS_SELF_INHIBITOR_WEIGHT" => "off", "DS_CYCLE_CARRIER_RELEASE" => "reaction") do
        o = Dict(k => (1.0, 1.0) for k in ("a_G", "a_Ph", "b_G", "b_Ph"))
        o["KG"] = (80.0, 1.0)
        DS.solve_steady_state(net, o, DS.SteadyStateParams(1.0, 0.1, 500, 1e-9, "penalty"))
    end
    @test net.nodes["a_r1"].reactome_id == net.nodes["b_r1"].reactome_id
    @test r2.diagnostics["cycle_carrier_release_added"] == 0
    for (k, v) in r.node_activities
        @test r2.node_activities[k] ≈ v atol = 1e-12
    end
end

@testset "multi-step cycles (amendment 2)" begin
    net = ring()
    # Fable's fixtures: the lumped ring reduces to the two-form rule
    for (obs, s, ss, cx) in ((Dict{String,Float64}(), 1.0, 1.0, 1.0),
                             (Dict("K" => 0.0), 1.11111, 0.0, 0.0),
                             (Dict("Ph" => 0.0), 0.0, 10.0, 0.0),
                             (Dict("K" => 80.0), 0.11236, 8.98876, 8.98876),
                             (Dict("K" => 80.0, "Ph" => 80.0), 1.0, 1.0, 80.0),
                             (Dict("G" => 80.0), 80.0, 80.0, 80.0))
        r = solve_ring(net; obs = obs)
        @test r.converged
        @test fold(r, "S") ≈ s atol = 1e-3
        @test fold(r, "Ss") ≈ ss atol = 1e-3
        @test fold(r, "D") ≈ ss atol = 1e-3             # the readout follows the modified state
        @test fold(r, "SE") ≈ cx atol = 1e-3            # an intermediate reads its path's flux
        @test fold(r, "r1") ≈ cx atol = 1e-3            # and so does each step's reaction node
        @test r.diagnostics["cycle_pools_solved"] == 1
        @test r.diagnostics["cycle_pools_multistep"] == 1
    end
    # carriers: with the kinase not pinned, its free form reads its OTHER
    # producer (KP), not the release step that closes its own loop
    r = solve_ring(net; obs = Dict("KG" => 80.0), pinned = ("G", "Ph"))
    @test r.converged
    @test fold(r, "K") ≈ 80.0 atol = 1e-6
    @test fold(r, "Ss") ≈ 8.98876 atol = 1e-3
    @test r.diagnostics["cycle_carriers"] == "on" && r.diagnostics["cycle_carriers_held"] == 1
    r0 = solve_ring(net; pinned = ("G", "Ph"))
    @test fold(r0, "K") ≈ 1.0 atol = 1e-6              # baseline exact with the carrier held
    roff = solve_ring(net; obs = Dict("KG" => 80.0), pinned = ("G", "Ph"), carriers = "off")
    @test roff.diagnostics["cycle_carriers"] == "off" && roff.diagnostics["cycle_carriers_held"] == 0
    @test !(fold(roff, "K") ≈ 80.0)                     # off: the release step still feeds K
    # copies of one step (a set-member kinase K2 beside K): the step's drive is
    # their mean, so K2 KO halves it, and each copy node reads its share
    rc = solve_ring(ring(; copies = true); obs = Dict("K2" => 0.0), pinned = ("G", "K", "Ph", "K2"))
    @test rc.converged
    @test fold(rc, "Ss") ≈ 0.5 / (0.1 * 0.5 + 0.9) atol = 1e-3
    @test fold(rc, "b1x") ≈ 0.0 atol = 1e-6
    @test fold(rc, "b1") ≈ 2 * fold(rc, "SE") atol = 1e-6
    rc0 = solve_ring(ring(; copies = true); pinned = ("G", "K", "Ph", "K2"))
    @test fold(rc0, "b1") ≈ 1.0 atol = 1e-9
    @test fold(rc0, "b1x") ≈ 1.0 atol = 1e-9
    # review of amendment 2, finding 4: an intermediate is v / k_cat. When the
    # step that EXITS it is blocked (C1, the modification step's co-input,
    # knocked out) the complex upstream accumulates; the one downstream empties
    rb = solve_ring(ring(; c1cat = true); obs = Dict("C1" => 0.0), pinned = ("G", "K", "Ph", "C1"))
    @test rb.converged
    @test fold(rb, "Ss") ≈ 0.0 atol = 1e-3
    @test fold(rb, "SE") ≈ fold(rb, "S") atol = 1e-6          # source fold x the other drives (1)
    @test fold(rb, "SE") > 1.0
    @test fold(rb, "SsE") ≈ 0.0 atol = 1e-6
    rb0 = solve_ring(ring(; c1cat = true); pinned = ("G", "K", "Ph", "C1"))
    @test fold(rb0, "SE") ≈ 1.0 atol = 1e-9
    # finding 2: a declared carrier that a pool writes is not a carrier
    r = solve_ring(ring(; extra_carriers = [("SE", "r1")]); obs = Dict("K" => 80.0))
    @test r.converged
    @test r.diagnostics["cycle_carrier_conflicts"] == 1
    @test fold(r, "Ss") ≈ 8.98876 atol = 1e-3
    # finding 5: a pinned intermediate makes the pool fall back, counted
    r = solve_ring(ring(); obs = Dict("SE" => 1.0))
    @test r.diagnostics["cycle_pools_solved"] == 0
    @test r.diagnostics["cycle_pools_pinned_fallback"] == 1
    # label independence on the ring
    rl(u) = "zz_" * u
    r1 = solve_ring(ring(); obs = Dict("K" => 3.0))
    r2 = solve_ring(ring(; rl = rl); obs = Dict("K" => 3.0), rl = rl)
    for k in ("S", "Ss", "SE", "SsP", "D")
        @test r1.node_activities[k] ≈ r2.node_activities[rl(k)] atol = 1e-12
    end
end

# specs/049: two pools that pass through one complex (MAPK1 and MAPK3 through a
# capped RAF:MAP2K:MAPKs node). S + T + K -> X (b), X -> Ss + Ts + K (r);
# Ss -> S (d1), Ts -> T (d2). X, b and r belong to both pools.
function twin(; rl = identity)
    ids = ["GS", "PS", "GT", "PT", "S", "T", "Ss", "Ts", "X", "K", "b", "r", "d1", "d2", "W", "D"]
    nodes = Dict(rl(k) => DS.NetworkNode(rl(k), "R-" * k, "protein", nothing, rl(k), BL) for k in ids)
    E = [("GS", "PS", "input"), ("PS", "S", "output"), ("GT", "PT", "input"), ("PT", "T", "output"),
         ("S", "b", "input"), ("T", "b", "input"), ("K", "b", "input"), ("b", "X", "output"),
         ("X", "r", "input"), ("r", "Ss", "output"), ("r", "Ts", "output"), ("r", "K", "output"),
         ("Ss", "d1", "input"), ("d1", "S", "output"), ("Ts", "d2", "input"), ("d2", "T", "output"),
         ("Ss", "W", "input"), ("Ts", "W", "input"), ("W", "D", "output")]
    edges = [edge(rl(s), rl(t), et; and = et == "input") for (s, t, et) in E]
    st(a, b, rs...) = (rl(a), rl(b), sort([rl(x) for x in rs]))
    pool(id, s, ss, d) = DS.CyclePool(id, sort([rl(s), rl(ss)]), rl(s),
        [DS.PoolPath(rl(s), rl(ss), [st(s, "X", "b"), st("X", ss, "r")], true),
         DS.PoolPath(rl(ss), rl(s), [st(ss, s, d)], false)], [rl("X")], Tuple{String, String}[])
    DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                       Dict{String, Set{String}}(), nothing, [pool("pool1", "S", "Ss", "d1"), pool("pool2", "T", "Ts", "d2")])
end

@testset "pools that share an intermediate write it once (specs/049)" begin
    pinned = ("GS", "GT", "K")
    r = solve_ring(twin(); pinned = pinned)
    @test r.converged
    @test r.diagnostics["cycle_pools_solved"] == 2
    @test r.diagnostics["cycle_shared_nodes"] == 3          # X, b and r
    for u in ("S", "T", "Ss", "Ts", "X", "b", "r", "D")
        @test fold(r, u) ≈ 1.0 atol = 1e-9                  # exactly baseline at rest
    end
    # the kinase knocked out: both modified forms and the shared complex go to zero
    r = solve_ring(twin(); obs = Dict("K" => 0.0), pinned = pinned)
    @test r.converged
    for u in ("Ss", "Ts", "X", "D")
        @test fold(r, u) ≈ 0.0 atol = 1e-6
    end
    # the kinase up: both rise, symmetrically, and the readout with them
    r = solve_ring(twin(); obs = Dict("K" => 80.0), pinned = pinned)
    @test r.converged
    @test fold(r, "Ss") > 2.0
    @test fold(r, "Ss") ≈ fold(r, "Ts") rtol = 1e-9
    @test fold(r, "D") > 2.0
    # one substrate's supply up: its own pool grows; the other pool's modified form
    # rises only through the shared binding step (T reads its partner S there)
    r = solve_ring(twin(); obs = Dict("GS" => 80.0), pinned = pinned)
    @test r.converged
    @test fold(r, "S") > 2.0 && fold(r, "Ss") > 2.0
    @test fold(r, "Ss") > fold(r, "Ts")
    # label-independent
    rl = u -> "q" * reverse(u)
    a = solve_ring(twin(); obs = Dict("K" => 80.0), pinned = pinned)
    b = solve_ring(twin(; rl = rl); obs = Dict("K" => 80.0), pinned = pinned, rl = rl)
    @test all(isapprox(a.node_activities[u], b.node_activities[rl(u)]; atol = 1e-12) for u in keys(a.node_activities))
end

@testset "an irreversible three-state ring is exactly baseline at rest" begin
    # A -> B -> C -> A, one step each, no reverse: sqrt(pi0_j/pi0_i) rates would
    # not make pi0 stationary; the baseline flow does
    ids = ["G", "P", "A", "B", "C", "E1", "E2", "E3", "t1", "t2", "t3"]
    nodes = Dict(k => DS.NetworkNode(k, "R-" * k, "protein", nothing, k, BL) for k in ids)
    E = [("G", "P", "input", true), ("P", "A", "output", false),
         ("A", "t1", "input", true), ("E1", "t1", "catalyst", true), ("t1", "B", "output", false),
         ("B", "t2", "input", true), ("E2", "t2", "catalyst", true), ("t2", "C", "output", false),
         ("C", "t3", "input", true), ("E3", "t3", "catalyst", true), ("t3", "A", "output", false)]
    edges = [edge(s, t, et; and = a) for (s, t, et, a) in E]
    p = [DS.CyclePool("pool1", ["A", "B", "C"], "A", [("A", "B", "t1"), ("B", "C", "t2"), ("C", "A", "t3")])]
    net = DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                             Dict{String, Set{String}}(), nothing, p)
    r = with_env("DS_CYCLE_MODE" => "balance", "DS_SELF_INHIBITOR_WEIGHT" => "off") do
        DS.solve_steady_state(net, Dict(k => (1.0, 1.0) for k in ("G", "E1", "E2", "E3")),
                              DS.SteadyStateParams(1.0, 0.1, 500, 1e-9, "penalty"))
    end
    @test r.converged
    for k in ("A", "B", "C", "t1", "t2", "t3")
        @test fold(r, k) ≈ 1.0 atol = 1e-9
    end
end

@testset "several base states: a protein resting in two copies of its starting complex (specs/081)" begin
    # p130 rests in two copies of its E2F complex, C1 and C2, assembled from roots S1, S2.
    # The kinase K moves it to M (F1, F2); PP2A (PP) returns it by R, which the pool
    # lists for both return paths (as the generator does for p130). D reads C1.
    ids = ["S1", "S2", "C1", "C2", "M", "K", "PP", "F1", "F2", "R", "X", "W", "D"]
    nodes = Dict(k => DS.NetworkNode(k, "R-" * k, "protein", nothing, k, BL) for k in ids)
    E = [("S1", "C1", "assembly", true), ("S2", "C2", "assembly", true),
         ("C1", "F1", "input", true), ("K", "F1", "catalyst", true), ("F1", "M", "output", false),
         ("C2", "F2", "input", true), ("K", "F2", "catalyst", true), ("F2", "M", "output", false),
         ("M", "R", "input", true), ("PP", "R", "catalyst", true), ("R", "X", "output", false),
         ("C1", "W", "input", true), ("W", "D", "output", false),
         # the recycle edges close the cycle structurally: free X is the starting copies' entity
         ("X", "S1", "recycle", false), ("X", "S2", "recycle", false)]
    edges = [edge(s, t, et; and = a) for (s, t, et, a) in E]
    paths = [DS.PoolPath("C1", "M", [("C1", "M", ["F1"])], true), DS.PoolPath("C2", "M", [("C2", "M", ["F2"])], true),
             DS.PoolPath("M", "C1", [("M", "C1", ["R"])], true), DS.PoolPath("M", "C2", [("M", "C2", ["R"])], true)]
    pool = DS.CyclePool("recycle1", ["C1", "C2", "M"], "C1", paths, String[], Tuple{String, String}[], ["C1", "C2"])
    net = DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                             Dict{String, Set{String}}(), nothing, [pool])
    go(obs) = with_env("DS_CYCLE_MODE" => "balance", "DS_CYCLE_PHI" => "0.1", "DS_SELF_INHIBITOR_WEIGHT" => "off") do
        o = Dict(k => (1.0, 1.0) for k in ("S1", "S2", "K", "PP"))
        for (k, v) in obs; o[k] = (v, 1.0); end
        DS.solve_steady_state(net, o, DS.SteadyStateParams(1.0, 0.1, 500, 1e-9, "penalty"))
    end
    rho = 0.1 / 0.9
    r0 = go(Dict{String, Float64}())
    @test r0.converged && r0.diagnostics["cycle_pools_solved"] == 1
    @test r0.diagnostics["recycle_edges"] == 2
    # a recycle edge is never an input: the starting copies have no update of their own
    @test !any(r -> r.target_uuid in ("S1", "S2"), DS.convert_to_reaction_network(net))
    # without the recycle edges the pool's nodes are not one component, so it is not solved
    norec = DS.ReactionNetwork(nodes, [e for e in edges if e.edge_type != "recycle"],
                               Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                               Dict{String, Set{String}}(), nothing, [pool])
    nr = with_env("DS_CYCLE_MODE" => "balance", "DS_SELF_INHIBITOR_WEIGHT" => "off") do
        DS.solve_steady_state(norec, Dict(k => (1.0, 1.0) for k in ("S1", "S2", "K", "PP")),
                              DS.SteadyStateParams(1.0, 0.1, 500, 1e-9, "penalty"))
    end
    @test nr.diagnostics["cycle_pools_solved"] == 0 && nr.diagnostics["recycle_edges"] == 0
    # unpinned starting copies stay at their own (baseline) value whatever X reads
    up2 = with_env("DS_CYCLE_MODE" => "balance", "DS_SELF_INHIBITOR_WEIGHT" => "off") do
        DS.solve_steady_state(net, Dict("K" => (80.0, 1.0), "PP" => (1.0, 1.0)),
                              DS.SteadyStateParams(1.0, 0.1, 500, 1e-9, "penalty"))
    end
    @test fold(up2, "S1") ≈ 1.0 atol = 1e-12
    @test fold(up2, "S2") ≈ 1.0 atol = 1e-12
    for k in ("C1", "C2", "M", "D")
        @test fold(r0, k) ≈ 1.0 atol = 1e-9          # both resting copies start at full weight
    end
    up = go(Dict("K" => 80.0))                       # more kinase: p130 leaves both complexes
    @test fold(up, "C1") ≈ fold(up, "C2") atol = 1e-9
    @test fold(up, "C1") < 0.25 && fold(up, "M") > 1.5 && fold(up, "D") < 0.25
    @test fold(go(Dict("PP" => 0.0)), "C1") < 0.01   # no phosphatase: none returns
    ko = go(Dict("K" => 0.0))                        # no kinase: all of it rests, (2 + ρ)/2
    @test fold(ko, "C1") ≈ (2 + rho) / 2 atol = 1e-6
    @test fold(ko, "M") < 0.01
    # the table reads several bases, and they survive a JSON round trip
    mktempdir() do dir
        write(joinpath(dir, "logic_network.csv"), "source_id,target_id,pos_neg,and_or,edge_type,stoichiometry\n")
        write(joinpath(dir, "pools.csv"), "pool_id,node_uuid,stable_id,role,is_base\n" *
              "recycle1,C1,R-C1,state,True\nrecycle1,C2,R-C2,state,True\nrecycle1,M,R-M,state,False\n")
        write(joinpath(dir, "pool_transitions.csv"),
              "pool_id,path_id,step,source_uuid,target_uuid,reaction_uuid,reaction_stid,enzyme_driven\n" *
              "recycle1,f1,1,C1,M,F1,R-F,True\nrecycle1,r1,1,M,C1,R,R-R,True\n")
        p = DS.parse_pools(joinpath(dir, "logic_network.csv"))
        @test length(p) == 1 && p[1].base == "C1" && p[1].bases == ["C1", "C2"]
        write(joinpath(dir, "pools.csv"), "pool_id,node_uuid,stable_id,role,is_base\n" *
              "recycle1,C1,R-C1,state,True\nrecycle1,Z,R-Z,intermediate,True\nrecycle1,M,R-M,state,False\n")
        @test_throws ArgumentError DS.parse_pools(joinpath(dir, "logic_network.csv"))
    end
    back = DS.pools_from_json(DS.JSON3.read(DS.JSON3.write(DS.pools_json(net))))
    @test back[1].bases == ["C1", "C2"] && back[1].base == "C1"
    old = DS.JSON3.read(DS.JSON3.write([Dict("id" => "p", "forms" => ["A", "B"], "base" => "A",
        "intermediates" => String[], "carriers" => [], "paths" => [Dict("from" => "A", "to" => "B", "enzyme" => true,
        "steps" => [["A", "B", ["F"]]])])]))
    @test DS.pools_from_json(old)[1].bases == ["A"]   # a table written before specs/081
end

@testset "configuration: DS_CYCLE_CARRIERS" begin
    for bad in ("On", "yes", "")
        with_env("DS_CYCLE_CARRIERS" => bad) do
            @test_throws ArgumentError DS.cycle_carriers()
        end
    end
    with_env("DS_CYCLE_CARRIERS" => nothing) do
        @test DS.cycle_carriers()
    end
end

@testset "the pool table survives a JSON round trip (code review 2026-10-02)" begin
    net = ring()
    raw = DS.JSON3.read(DS.JSON3.write(DS.pools_json(net)))
    back = DS.pools_from_json(raw)
    @test length(back) == length(net.pools)
    for (a, b) in zip(net.pools, back)
        @test (a.id, a.forms, a.base, a.intermediates, a.carriers) == (b.id, b.forms, b.base, b.intermediates, b.carriers)
        @test [(q.from, q.to, q.steps, q.enzyme) for q in a.paths] == [(q.from, q.to, q.steps, q.enzyme) for q in b.paths]
    end
    # the round-tripped network solves exactly as the original
    net2 = DS.ReactionNetwork(net.nodes, net.edges, net.set_mappings, net.cofactor_stids,
                              net.containment, net.drug_stids, back)
    for obs in (Dict("K" => 80.0), Dict("Ph" => 0.0))
        r1 = solve_ring(net; obs = obs); r2 = solve_ring(net2; obs = obs)
        @test r1.node_activities == r2.node_activities
        @test r2.diagnostics["cycle_rule"] == "balance"
    end
    # "no table" stays "no table"; malformed input is an error, not a smaller table
    @test DS.pools_json(fixture(; pools = nothing)) === nothing
    @test DS.pools_from_json(nothing) === nothing
    @test_throws ArgumentError DS.pools_from_json("x")
    @test_throws ArgumentError DS.pools_from_json(DS.JSON3.read("""[{"id":1,"forms":["a"],"base":"a","intermediates":[],"carriers":[],"paths":[]}]"""))
    @test_throws ArgumentError DS.pools_from_json(DS.JSON3.read("""[{"id":"p","forms":["a"],"base":"a","intermediates":[],"carriers":[["c"]],"paths":[]}]"""))
    @test_throws ArgumentError DS.pools_from_json(DS.JSON3.read("""[{"id":"p","forms":["a"],"base":"a","intermediates":[],"carriers":[]}]"""))
end

@testset "supply reached only by regulation stays a supply (specs/071)" begin
    # B (the modified form) induces P, the reaction that makes the protein: the
    # pool reaches its own supply, but only through a regulator edge. Under the
    # default the review-finding-6 test drops P from the supply, so knocking
    # out the gene cannot move the pool (TP53 inducing MDM2, specs/070 §5).
    net = fixture(extra = [("B", "P", "regulator", true)])
    solve_r(obs, reach) = with_env("DS_POOL_SUPPLY_REACH" => reach) do
        solve(net; obs = obs)
    end
    for reach in (nothing, "all", "material")
        r = solve_r(Dict{String,Float64}(), reach)
        @test r.converged
        @test fold(r, "A") ≈ 1.0 atol = 1e-6
        @test fold(r, "B") ≈ 1.0 atol = 1e-6
        @test r.diagnostics["pool_supply_reach"] == something(reach, "all")
    end
    ko_all = solve_r(Dict("G" => 0.0), "all")
    @test fold(ko_all, "A") ≈ 1.0 atol = 1e-6       # the gene cannot reach its own pool
    @test fold(ko_all, "B") ≈ 1.0 atol = 1e-6
    @test ko_all.diagnostics["cycle_supply_regulated_kept"] == 0
    ko_mat = solve_r(Dict("G" => 0.0), "material")
    @test ko_mat.converged
    @test fold(ko_mat, "A") ≈ 0.0 atol = 1e-6       # no protein made, none in either form
    @test fold(ko_mat, "B") ≈ 0.0 atol = 1e-6
    @test fold(ko_mat, "D") ≈ 0.0 atol = 1e-6
    @test ko_mat.diagnostics["cycle_supply_regulated_kept"] == 1
    # a supply the pool feeds with MATERIAL is still dropped under `material`
    netm = fixture(extra = [("B", "P", "input", true)])
    r = with_env("DS_POOL_SUPPLY_REACH" => "material") do
        solve(netm; obs = Dict("G" => 0.0))
    end
    @test r.diagnostics["cycle_supply_regulated_kept"] == 0
    @test fold(r, "A") ≈ fold(with_env(() -> solve(netm; obs = Dict("G" => 0.0)), "DS_POOL_SUPPLY_REACH" => "all"), "A") atol = 1e-9
    # the edge classifier: a pair that is also an input carries material
    n2 = fixture(extra = [("B", "P", "regulator", true), ("EF", "F", "input", true), ("EB", "W", "regulator", false)])
    sp = DS.pool_signal_edge_pairs(n2)
    @test ("B", "P") in sp
    @test ("EB", "R") in sp                         # catalyst only
    @test !(("EF", "F") in sp)                      # catalyst AND input
    @test !(("EB", "W") in sp)                      # negative regulators are not followed at all
    @test !(("A", "F") in sp)
    for bad in ("", "Material", "regulator")
        with_env("DS_POOL_SUPPLY_REACH" => bad) do
            @test_throws ArgumentError solve(net)
        end
    end
end

@testset "every copy of a listed release reaction is a release step (specs/072)" begin
    # r1x is a second copy of the release reaction r1 (same stable id): it
    # releases K from an intermediate SsEx that is not in the pool (USP7 from a
    # second DAXX:PolyUb-MDM2:USP7 variant, specs/071). Only r1 is listed.
    function ring_copy()
        n = ring()
        nodes = copy(n.nodes)
        for k in ("V", "SsEx"); nodes[k] = DS.NetworkNode(k, "R-" * k, "protein", nothing, k, BL); end
        nodes["r1x"] = DS.NetworkNode("r1x", "R-r1", "reaction", nothing, "r1x", BL)
        edges = vcat(n.edges, [edge("V", "SsEx", "output"; and = false), edge("SsEx", "r1x", "input"),
                               edge("r1x", "K", "output"; and = false)])
        DS.ReactionNetwork(nodes, edges, n.set_mappings, n.cofactor_stids, n.containment, nothing, n.pools)
    end
    net = ring_copy()
    rr(obs, mode) = with_env("DS_CYCLE_CARRIER_RELEASE" => mode) do
        solve_ring(net; obs = obs, pinned = ("G", "Ph", "V"))
    end
    for mode in (nothing, "listed", "reaction")
        r = rr(Dict{String,Float64}(), mode)
        @test r.converged
        @test fold(r, "K") ≈ 1.0 atol = 1e-6
        @test r.diagnostics["cycle_carrier_release"] == something(mode, "listed")
    end
    listed = rr(Dict("V" => 80.0), "listed")
    @test listed.diagnostics["cycle_carrier_release_added"] == 0
    @test fold(listed, "K") > 2.0                    # the unlisted copy feeds the carrier
    react = rr(Dict("V" => 80.0), "reaction")
    @test react.converged
    @test react.diagnostics["cycle_carrier_release_added"] == 1
    @test fold(react, "K") ≈ 1.0 atol = 1e-6         # both copies left out: K reads KP only
    @test fold(react, "Ss") ≈ 1.0 atol = 1e-6
    # a producer with another stable id is still read
    kg = rr(Dict("V" => 80.0, "KG" => 80.0), "reaction")
    @test fold(kg, "K") ≈ 80.0 atol = 1e-6
    # carriers off: nothing is added and the mode is reported as listed
    off = with_env("DS_CYCLE_CARRIER_RELEASE" => "reaction") do
        solve_ring(net; obs = Dict("V" => 80.0), pinned = ("G", "Ph", "V"), carriers = "off")
    end
    @test off.diagnostics["cycle_carrier_release"] == "listed"
    @test off.diagnostics["cycle_carrier_release_added"] == 0
    for bad in ("", "Reaction", "stid")
        with_env("DS_CYCLE_CARRIER_RELEASE" => bad) do
            @test_throws ArgumentError DS.cycle_carrier_release()
            @test_throws ArgumentError solve_ring(net)
        end
    end
end

@testset "a non-converged pool component is re-solved jointly (specs/073)" begin
    # The iteration is stopped after `iters` sweeps, so the pool component ends
    # above tolerance; the converged solve (500 sweeps) is the reference.
    function capped(net, obs, iters, pools_mode; fallback = "lm", pinned = ("G", "K", "Ph"))
        with_env("DS_CYCLE_MODE" => "balance", "DS_CYCLE_PHI" => "0.1", "DS_CYCLE_CARRIERS" => "on",
                 "DS_SELF_INHIBITOR_WEIGHT" => "off", "DS_NONCONV_FALLBACK" => fallback,
                 "DS_NONCONV_POOLS" => pools_mode) do
            o = Dict(k => (1.0, 1.0) for k in pinned)
            for (k, v) in obs; o[k] = (v, 1.0); end
            DS.solve_steady_state(net, o, DS.SteadyStateParams(1.0, 0.1, iters, 1e-9, "penalty"))
        end
    end
    with_env("DS_NONCONV_POOLS" => nothing) do
        @test DS.nonconv_pools_mode() == "joint"     # default since 2026-10-18
    end
    # A pool alone is solved exactly in one sweep; a loop from its readout D
    # back into one of its steps makes the component iterate.
    function looped(n, s, t)
        DS.ReactionNetwork(n.nodes, vcat(n.edges, [edge(s, t, "input")]), n.set_mappings, n.cofactor_stids,
                           n.containment, nothing, n.pools)
    end
    rg = looped(ring(), "D", "c1")
    cases = [(rg, Dict("G" => 3.0), ("G", "K", "Ph")),
             (rg, Dict("G" => 80.0), ("G", "K", "Ph")),
             (rg, Dict("K" => 0.2), ("G", "K", "Ph")),
             (looped(ring(; copies = true), "D", "c1"), Dict("K2" => 5.0), ("G", "K", "Ph", "K2")),
             (looped(fixture(), "D", "F"), Dict("EF" => 3.0), ("G", "EF", "EB", "Z")),
             (looped(fixture(), "D", "F"), Dict("EB" => 0.2), ("G", "EF", "EB", "Z"))]
    for (net, obs, pinned) in cases
        ref = capped(net, obs, 500, "skip"; pinned = pinned)
        @test ref.converged
        skip = capped(net, obs, 3, "skip"; pinned = pinned)
        @test !skip.converged
        @test skip.diagnostics["nonconv_pools"] == "skip"
        @test skip.diagnostics["nonconv_skipped_pools"] == 1
        @test skip.diagnostics["nonconv_pool_components"] == 0
        joint = capped(net, obs, 3, "joint"; pinned = pinned)
        @test joint.diagnostics["nonconv_pools"] == "joint"
        @test joint.diagnostics["nonconv_skipped_pools"] == 0
        @test joint.diagnostics["nonconv_pool_components"] == 1
        @test joint.diagnostics["nonconv_pool_improved"] == 1
        @test joint.diagnostics["nonconv_pool_rescued"] == 1
        @test joint.converged
        @test maximum(abs(joint.node_activities[k] - ref.node_activities[k]) for k in keys(ref.node_activities)) < 1e-7
        # converged solves are untouched, bit for bit
        @test capped(net, obs, 500, "joint"; pinned = pinned).node_activities == ref.node_activities
    end
    # with the fallback off, the pool rule is reported as skip and nothing runs
    off = capped(rg, Dict("G" => 3.0), 3, "joint"; fallback = "off")
    @test !off.converged
    @test off.diagnostics["nonconv_pools"] == "skip"
    @test off.diagnostics["nonconv_pool_components"] == 0
    for bad in ("", "Joint", "lm")
        with_env("DS_NONCONV_POOLS" => bad) do
            @test_throws ArgumentError DS.nonconv_pools_mode()
            @test_throws ArgumentError solve_ring(ring())
        end
    end
end

@testset "a carrier re-evaluated in another pool's drives reads by the carrier rule (specs/075)" begin
    # pool2 (A <-> B) returns B to A through R2x, catalysed by K, the carrier of
    # the ring's pool. B feeds the ring's binding step, so K is pool-fed for
    # pool2 and its drives re-evaluate it (specs/039 amendment 5). K's whole
    # reaction reads its release step r1, which carries the ring's flux; the
    # carrier rule leaves r1 out (USP7 in the TP53 pool, read 50.5x).
    function two()
        n = ring()
        nodes = copy(n.nodes)
        for k in ("A", "B", "EF2", "F2x", "R2x")
            nodes[k] = DS.NetworkNode(k, "R-" * k, k in ("F2x", "R2x") ? "reaction" : "protein", nothing, k, BL)
        end
        E = [edge("A", "F2x", "input"), edge("EF2", "F2x", "catalyst"), edge("F2x", "B", "output"; and = false),
             edge("B", "R2x", "input"), edge("K", "R2x", "catalyst"), edge("R2x", "A", "output"; and = false),
             edge("B", "b1", "input")]
        p2 = DS.CyclePool("pool2", sort(["A", "B"]), "A", [("A", "B", "F2x"), ("B", "A", "R2x")])
        DS.ReactionNetwork(nodes, vcat(n.edges, E), n.set_mappings, n.cofactor_stids, n.containment, nothing,
                           vcat(n.pools, [p2]))
    end
    net = two()
    rr(g, mode) = with_env("DS_CYCLE_REEVAL_CARRIERS" => mode) do
        solve_ring(net; obs = Dict("G" => g), pinned = ("G", "Ph", "EF2"))
    end
    for mode in (nothing, "full", "carrier")
        r = rr(1.0, mode)
        @test r.converged
        @test r.diagnostics["cycle_reeval_carrier_rule"] == something(mode, "carrier")   # default since 2026-10-18
        @test r.diagnostics["cycle_reeval_carriers"] == 1        # counted under both modes
        for k in ("A", "B", "K", "Ss"); @test fold(r, k) ≈ 1.0 atol = 1e-6; end
    end
    for g in (80.0, 0.2)
        full = rr(g, "full")
        car = rr(g, "carrier")
        @test full.converged && car.converged
        @test fold(full, "K") ≈ 1.0 atol = 1e-6
        @test fold(car, "K") ≈ 1.0 atol = 1e-6
        # full: pool2's drive reads K's release flux, so B moves with G
        @test abs(log(fold(full, "B"))) > 0.3
        @test abs(log(fold(full, "Ss") / g)) > 0.3
        # carrier: K is 1 wherever it is read, so pool2 stays at rest and the
        # ring follows its supply
        @test fold(car, "A") ≈ 1.0 atol = 1e-6
        @test fold(car, "B") ≈ 1.0 atol = 1e-6
        @test fold(car, "Ss") ≈ g rtol = 1e-6
    end
    # the ring alone has no other pool to re-evaluate its carriers
    @test solve_ring(ring()).diagnostics["cycle_reeval_carriers"] == 0
    # carriers off: nothing is counted and the rule is reported as off
    off = with_env("DS_CYCLE_REEVAL_CARRIERS" => "carrier") do
        solve_ring(net; obs = Dict("G" => 80.0), pinned = ("G", "Ph", "EF2"), carriers = "off")
    end
    @test off.diagnostics["cycle_reeval_carrier_rule"] == "off"
    @test off.diagnostics["cycle_reeval_carriers"] == 0
    for bad in ("", "Carrier", "rule")
        with_env("DS_CYCLE_REEVAL_CARRIERS" => bad) do
            @test_throws ArgumentError DS.cycle_reeval_carriers()
            @test_throws ArgumentError solve_ring(net)
        end
    end
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
        @test [(q.from, q.to, q.steps) for q in p[1].paths] ==
              [("u-a", "u-b", [("u-a", "u-b", ["u-f"])]), ("u-b", "u-a", [("u-b", "u-a", ["u-r"])])]
        @test all(q -> q.enzyme, p[1].paths) && isempty(p[1].intermediates)
        write(joinpath(dir, "pools.csv"), "pool_id,node_uuid,stable_id,is_base\npool1,u-a,R-A,False\n")
        @test_throws ArgumentError DS.parse_pools(ln)   # a pool must have a base form
        # ... a malformed table is an error when the rule uses it (balance, the
        # default), and dropped with a warning only when the rule is off
        write(ln, "source_id,target_id,pos_neg,and_or,edge_type\n")
        write(joinpath(dir, "stid_to_uuid_mapping.csv"), "stable_id,uuid\n")
        with_env("DS_CYCLE_MODE" => "off") do
            net = DS.parse_complete_network(ln, joinpath(dir, "stid_to_uuid_mapping.csv"))
            @test net.pools === nothing
        end
        with_env("DS_CYCLE_MODE" => nothing) do
            @test_throws ArgumentError DS.parse_complete_network(ln, joinpath(dir, "stid_to_uuid_mapping.csv"))
        end
        with_env("DS_CYCLE_MODE" => "balance") do
            @test_throws ArgumentError DS.parse_complete_network(ln, joinpath(dir, "stid_to_uuid_mapping.csv"))
        end
        write(ln, "")
        # amendment 2: roles, one row per step, carriers
        write(joinpath(dir, "pools.csv"), "pool_id,node_uuid,stable_id,role,is_base\n" *
              "pool1,u-s,R-S,state,True\npool1,u-t,R-T,state,False\npool1,u-c,R-C,intermediate,False\n")
        write(joinpath(dir, "pool_transitions.csv"),
              "pool_id,path_id,step,source_uuid,target_uuid,reaction_uuid,reaction_stid,enzyme_driven\n" *
              "pool1,p2,1,u-t,u-s,u-h,R-H,False\npool1,p1,2,u-c,u-t,u-r,R-R,True\npool1,p1,1,u-s,u-c,u-b,R-B,True\n")
        write(joinpath(dir, "pool_carriers.csv"), "pool_id,carrier_uuid,release_reaction_uuid\npool1,u-e,u-r\n")
        p = DS.parse_pools(ln)
        @test p[1].forms == ["u-s", "u-t"] && p[1].intermediates == ["u-c"]
        @test [(q.from, q.to, length(q.steps), q.enzyme) for q in p[1].paths] ==
              [("u-s", "u-t", 2, true), ("u-t", "u-s", 1, false)]
        @test p[1].paths[1].steps == [("u-s", "u-c", ["u-b"]), ("u-c", "u-t", ["u-r"])]
        # copies: several rows for one (path, step) are one step with its copies
        write(joinpath(dir, "pool_transitions.csv"),
              "pool_id,path_id,step,source_uuid,target_uuid,reaction_uuid,reaction_stid,enzyme_driven\n" *
              "pool1,p1,1,u-s,u-c,u-b2,R-B,True\npool1,p1,1,u-s,u-c,u-b,R-B,True\npool1,p1,2,u-c,u-t,u-r,R-R,True\n" *
              "pool1,p2,1,u-t,u-s,u-h,R-H,False\n")
        p = DS.parse_pools(ln)
        @test p[1].paths[1].steps == [("u-s", "u-c", ["u-b", "u-b2"]), ("u-c", "u-t", ["u-r"])]
        write(joinpath(dir, "pool_transitions.csv"),
              "pool_id,path_id,step,source_uuid,target_uuid,reaction_uuid,reaction_stid,enzyme_driven\n" *
              "pool1,p1,1,u-s,u-c,u-b,R-B,True\npool1,p1,1,u-s,u-t,u-b2,R-B,True\npool1,p1,2,u-c,u-t,u-r,R-R,True\n")
        @test_throws ArgumentError DS.parse_pools(ln)    # copies of one step must join the same nodes
        @test p[1].carriers == [("u-e", "u-r")]
        # a path whose steps do not chain is an error, not a silent pool
        write(joinpath(dir, "pool_transitions.csv"),
              "pool_id,path_id,step,source_uuid,target_uuid,reaction_uuid,reaction_stid,enzyme_driven\n" *
              "pool1,p1,1,u-s,u-c,u-b,R-B,True\npool1,p1,2,u-x,u-t,u-r,R-R,True\n")
        @test_throws ArgumentError DS.parse_pools(ln)
    end
end

end
