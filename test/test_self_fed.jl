# specs/040: a node's own downstream multiplied back into it.
#
# Rule A (DS_SELF_FED_MODE=entry): an input of an iterated component that is
# fed by nothing except what it feeds reads its component-entry value, i.e.
# baseline. Rule B (DS_SELF_INHIBITOR_LEAVES=1): specs/022's self-contained
# inhibitor test widened to a shared leaf whose carrier reaches both.
#
# Fixtures are the three RAF amplifiers in miniature:
#   A1  a recycled set-pool member:   X, pool P -> R1 -> D -> R2 -> M -(set_member)-> P
#   A2  bundle leaves fed downstream: X, L1, L2 -> R1 -> C -> R2 -> L1, L2
#   B   a leaf-sharing inhibitor:     X -> A (input of T); X -> A2; A2 + Q -> I -| T
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
edge(s, t, pos, et; and = true) = DS.LogicNetworkEdge(s, t, and, pos, 1.0, et)
net(nodes, edges; cont = Dict{String, Set{String}}(), cof = Set{String}()) =
    DS.ReactionNetwork(Dict(node.(nodes)), edges, Dict{String, DS.SetExpansionMapping}(), cof, cont)

"A1: M is produced only by R2, downstream of R1, and feeds R1 back through pool P."
fixture_pool() = net(["X", "M0", "M", "P", "R1", "D", "R2"],
    [edge("X", "R1", true, "input"), edge("P", "R1", true, "catalyst"),
     edge("R1", "D", true, "output"; and = false),
     edge("D", "R2", true, "input"),
     edge("R2", "M", true, "output"; and = false),
     edge("M", "P", true, "set_member"; and = false), edge("M0", "P", true, "set_member"; and = false)])

"A2: the leaves L1, L2 are AND inputs of R1 and are released only by R2. With `outside`, Z also produces L1."
fixture_leaves(; outside = false) = net(["X", "L1", "L2", "R1", "C", "R2", "Z"],
    vcat([edge("X", "R1", true, "input"), edge("L1", "R1", true, "input"), edge("L2", "R1", true, "input"),
          edge("R1", "C", true, "output"; and = false),
          edge("C", "R2", true, "input"),
          edge("R2", "L1", true, "output"; and = false), edge("R2", "L2", true, "output"; and = false)],
         outside ? [edge("Z", "L1", true, "output"; and = false)] : DS.LogicNetworkEdge[]))

"B: T = A, inhibited by I = A2 + Q, where A and A2 are sibling products of X and share the leaf a."
const CONT_B = Dict("R-I" => Set(["R-A2", "R-a", "R-q"]), "R-A2" => Set(["R-a"]),
                    "R-A" => Set(["R-a"]), "R-X" => Set(["R-a"]), "R-Q" => Set(["R-q"]))
fixture_inh(; cont = CONT_B, cof = Set{String}()) = net(["X", "Q", "A", "A2", "I", "T"],
    [edge("X", "A", true, "output"; and = false), edge("X", "A2", true, "output"; and = false),
     edge("A2", "I", true, "assembly"), edge("Q", "I", true, "assembly"),
     edge("A", "T", true, "input"), edge("I", "T", false, "regulator"; and = false)];
    cont = cont, cof = cof)

const PARAMS = DS.SteadyStateParams(1.0, 0.1, 500, 1e-9, "penalty")
function solve(network, obs; env...)
    pairs = [String(k) => v for (k, v) in env]
    with_env(pairs...) do
        DS.solve_steady_state(network, obs, PARAMS)
    end
end
fold(r, u) = r.node_activities[u] / BL
obs(pairs...) = Dict{String, Tuple{Float64, Float64}}(u => (v, 1.0) for (u, v) in pairs)

"A2 plus N: a node with NO activator input (only a depletion from C) inside the component."
fixture_noact() = net(["X", "L1", "L2", "R1", "C", "R2", "N"],
    [edge("X", "R1", true, "input"), edge("L1", "R1", true, "input"), edge("L2", "R1", true, "input"),
     edge("N", "R1", true, "input"),
     edge("R1", "C", true, "output"; and = false),
     edge("C", "R2", true, "input"),
     edge("R2", "L1", true, "output"; and = false), edge("R2", "L2", true, "output"; and = false),
     edge("C", "N", false, "depletion")])

@testset "self-fed inputs and leaf-shared inhibitors (specs/040)" begin

@testset "rule A: entries are only what the pre-registration lists (review of specs/040)" begin
    # A node with no activator input is not an entry by that fact alone: at rest
    # nothing is off baseline, so nothing is self-fed.
    r = solve(fixture_noact(), obs("X" => 1.0); DS_SELF_FED_MODE = "entry")
    @test r.diagnostics["self_fed_nodes"] == 0
    @test r.diagnostics["self_fed_edges_held"] == 0
    @test all(v -> isapprox(v, BL; atol = 1e-12), values(r.node_activities))
end


@testset "mode readers reject typos, defaults are off" begin
    with_env("DS_SELF_FED_MODE" => nothing) do
        @test DS.self_fed_mode() == "off"
    end
    with_env("DS_SELF_FED_MODE" => "entry") do
        @test DS.self_fed_mode() == "entry"
    end
    for bad in ("entri", "on", "1", "")
        with_env("DS_SELF_FED_MODE" => bad) do
            @test_throws ArgumentError DS.self_fed_mode()
            @test_throws ArgumentError solve(fixture_pool(), obs("X" => 2.0))
        end
    end
    with_env("DS_SELF_INHIBITOR_LEAVES" => nothing) do
        @test DS.self_inhibitor_leaves() == false
    end
    for bad in ("yes", "on", "2", "")
        with_env("DS_SELF_INHIBITOR_LEAVES" => bad) do
            @test_throws ArgumentError DS.self_inhibitor_leaves()
            @test_throws ArgumentError solve(fixture_inh(), obs("X" => 2.0))
        end
    end
end

@testset "A1: a recycled set-pool member reads baseline" begin
    r_on = solve(fixture_pool(), obs("X" => 2.0); DS_SELF_FED_MODE = "entry")
    # P = product(M0 = 1, M held at 1) = 1, so the step carries X alone.
    @test fold(r_on, "D") ≈ 2.0 rtol = 1e-6
    @test fold(r_on, "M") ≈ 2.0 rtol = 1e-6
    @test fold(r_on, "P") ≈ 1.0 rtol = 1e-6
    @test r_on.diagnostics["self_fed_rule"] == "entry"
    @test r_on.diagnostics["self_fed_nodes"] == 1
    @test r_on.diagnostics["self_fed_edges_held"] == 1
    # Off: the loop R1 -> D -> R2 -> M -> P -> R1 has gain 1 and a push rails it.
    r_off = solve(fixture_pool(), obs("X" => 2.0); DS_SELF_FED_MODE = "off")
    @test !(fold(r_off, "D") ≈ 2.0)
    @test r_off.diagnostics["self_fed_rule"] == "off"
    @test r_off.diagnostics["self_fed_edges_held"] == 0
    # A knockdown reads DOWN, not to the all-zero root through the loop.
    r_ko = solve(fixture_pool(), obs("X" => 0.5); DS_SELF_FED_MODE = "entry")
    @test fold(r_ko, "D") ≈ 0.5 rtol = 1e-6
end

@testset "A1: a pinned member is an entry and reads live" begin
    # M pinned at 0.5x: the signal enters the loop at M, so M is not self-fed.
    r = solve(fixture_pool(), obs("X" => 2.0, "M" => 0.5); DS_SELF_FED_MODE = "entry")
    @test fold(r, "P") ≈ 0.5 rtol = 1e-6
    @test fold(r, "D") ≈ 1.0 rtol = 1e-6
    @test r.diagnostics["self_fed_edges_held"] == 0
end

@testset "A2: bundle leaves fed only by the step's own downstream" begin
    r_on = solve(fixture_leaves(), obs("X" => 2.0); DS_SELF_FED_MODE = "entry")
    @test fold(r_on, "C") ≈ 2.0 rtol = 1e-6
    @test fold(r_on, "L1") ≈ 2.0 rtol = 1e-6
    @test r_on.diagnostics["self_fed_nodes"] == 2
    @test r_on.diagnostics["self_fed_edges_held"] == 2
    r_off = solve(fixture_leaves(), obs("X" => 2.0); DS_SELF_FED_MODE = "off")
    @test !(fold(r_off, "C") ≈ 2.0)        # gain-2 loop rails
    # An outside producer AT baseline carries no signal: L1 is still self-fed.
    r_z = solve(fixture_leaves(outside = true), obs("X" => 2.0, "Z" => 1.0); DS_SELF_FED_MODE = "entry")
    @test r_z.diagnostics["self_fed_edges_held"] == 2
    @test fold(r_z, "C") ≈ 2.0 rtol = 1e-6
    # An outside producer OFF baseline makes L1 the (only) entry. L2 is held
    # as before, and so is R2: its supply path starts at L1, its own product,
    # so R2's edges into its products read baseline and L1 reads its
    # exogenous part, mean(1, 3) = 2 -- amendment 5's re-evaluation.
    r_z2 = solve(fixture_leaves(outside = true), obs("X" => 1.0, "Z" => 3.0); DS_SELF_FED_MODE = "entry")
    @test r_z2.diagnostics["self_fed_edges_held"] == 3
    @test fold(r_z2, "L1") ≈ 2.0 rtol = 1e-6
    @test fold(r_z2, "C") ≈ 2.0 rtol = 1e-6
end

@testset "held edges: the final residual reads them as the iteration did" begin
    # Before the fix the verification pass read a held edge live and reported
    # a residual of up to 1 on a converged solve (3 of 8 RAF cases).
    for f in (fixture_pool(), fixture_leaves())
        r = solve(f, obs("X" => 2.0); DS_SELF_FED_MODE = "entry")
        @test r.converged
        @test r.final_residual < PARAMS.tolerance
    end
end

@testset "off is byte-identical to unset; baseline is exact" begin
    for f in (fixture_pool(), fixture_leaves(), fixture_inh())
        a = solve(f, obs("X" => 2.0); DS_SELF_FED_MODE = nothing, DS_SELF_INHIBITOR_LEAVES = nothing)
        b = solve(f, obs("X" => 2.0); DS_SELF_FED_MODE = "off", DS_SELF_INHIBITOR_LEAVES = "0")
        @test a.node_activities == b.node_activities
        # At rest every held edge reads fold 1: nothing moves, on or off.
        for mode in ("off", "entry"), lv in ("0", "1")
            r = solve(f, obs(); DS_SELF_FED_MODE = mode, DS_SELF_INHIBITOR_LEAVES = lv)
            @test all(abs(v - BL) < 1e-12 for v in values(r.node_activities))
        end
    end
end

@testset "the self-fed set does not depend on node labels" begin
    f = fixture_pool()
    relabel = Dict(u => "zz_" * reverse(u) for u in keys(f.nodes))
    nodes2 = Dict(relabel[u] => DS.NetworkNode(relabel[u], n.reactome_id, n.entity_type,
                                             n.original_set_id, n.display_name, n.baseline)
                  for (u, n) in f.nodes)
    edges2 = [DS.LogicNetworkEdge(relabel[e.parent_uuid], relabel[e.child_uuid], e.is_and,
                                  e.is_positive, e.stoichiometry, e.edge_type) for e in f.edges]
    f2 = DS.ReactionNetwork(nodes2, edges2, Dict{String, DS.SetExpansionMapping}())
    r1 = solve(f, obs("X" => 2.0); DS_SELF_FED_MODE = "entry")
    r2 = solve(f2, obs(relabel["X"] => 2.0); DS_SELF_FED_MODE = "entry")
    @test r2.diagnostics["self_fed_edges_held"] == r1.diagnostics["self_fed_edges_held"]
    for u in keys(f.nodes)
        @test r2.node_activities[relabel[u]] ≈ r1.node_activities[u] rtol = 1e-9
    end
end

@testset "B: a leaf-sharing inhibitor built from the input's sibling" begin
    # Strict (specs/022): I contains A2, not A, and A does not reach I: no pair.
    @test isempty(DS.self_contained_inhibitor_map(fixture_inh(); leaves = false))
    # Leaf clause: A and I share the leaf a; X carries a and reaches both.
    m = DS.self_contained_inhibitor_map(fixture_inh(); leaves = true)
    @test m == Dict(("I", "T") => ["A"])
    with_env("DS_SELF_INHIBITOR_LEAVES" => "1") do
        @test DS.self_inhibitor_setup(fixture_inh())[3] == 1
    end
    # Off: I tracks A exactly and cancels it (2 * 1/2 = 1).
    r_off = solve(fixture_inh(), obs("X" => 2.0); DS_SELF_INHIBITOR_LEAVES = "0", DS_SELF_INHIBITOR_WEIGHT = "0.1")
    @test fold(r_off, "I") ≈ 2.0 rtol = 1e-6
    @test fold(r_off, "T") ≈ 1.0 rtol = 1e-6
    @test r_off.diagnostics["self_inhibitor_leaves"] == "off"
    @test r_off.diagnostics["self_inhibitor_leaf_pairs"] == 0
    # On, w = 0.1: the tracking part is kept at weight w, so T reads 2^0.9.
    r_on = solve(fixture_inh(), obs("X" => 2.0); DS_SELF_INHIBITOR_LEAVES = "1", DS_SELF_INHIBITOR_WEIGHT = "0.1")
    @test fold(r_on, "T") ≈ 2.0^0.9 rtol = 1e-6
    @test r_on.diagnostics["self_inhibitor_leaves"] == "on"
    @test r_on.diagnostics["self_inhibitor_leaf_pairs"] == 1
    @test r_on.diagnostics["self_inhibitors"] == 1
    # A knockdown of X reads down at T, damped, not inverted.
    r_kd = solve(fixture_inh(), obs("X" => 0.5); DS_SELF_INHIBITOR_LEAVES = "1", DS_SELF_INHIBITOR_WEIGHT = "0.1")
    @test fold(r_kd, "T") ≈ 0.5^0.9 rtol = 1e-6
end

@testset "B: what the leaf clause does not flag" begin
    # No carrier of the shared leaf reaches both (X does not carry a): not built
    # from the same species, so containment alone is not enough.
    cont = Dict("R-I" => Set(["R-A2", "R-a", "R-q"]), "R-A2" => Set(["R-a"]), "R-A" => Set(["R-a"]))
    @test isempty(DS.self_contained_inhibitor_map(fixture_inh(cont = cont); leaves = true))
    # The shared leaf is a cofactor: not a shared species.
    @test isempty(DS.self_contained_inhibitor_map(fixture_inh(cof = Set(["R-a"])); leaves = true))
    # Without a containment table the rule is inert and says so.
    bare = fixture_inh(cont = Dict{String, Set{String}}())
    with_env("DS_SELF_INHIBITOR_LEAVES" => "1") do
        m, status, n = DS.self_inhibitor_setup(bare)
        @test isempty(m) && startswith(status, "inert") && n == 0
    end
    # With the weight off, the leaf clause has nothing to act on.
    with_env("DS_SELF_INHIBITOR_LEAVES" => "1", "DS_SELF_INHIBITOR_WEIGHT" => "off") do
        m, status, n = DS.self_inhibitor_setup(fixture_inh())
        @test isempty(m) && status == "off" && n == 0
    end
end

@testset "self_fed_nodes, by hand" begin
    # ring e -> t -> a -> b -> t with one entry t: only b (Out = {t}) is self-fed.
    act = Dict(1 => [2], 2 => [3], 3 => [1])
    @test DS.self_fed_nodes([1, 2, 3], act, Set{Int}(), Set{Int}(), Set([1]), Set{Int}()) == Set([3])
    # a second entry at a: b is reached from a without passing t: nothing is self-fed.
    @test isempty(DS.self_fed_nodes([1, 2, 3], act, Set{Int}(), Set{Int}(), Set([1, 2]), Set{Int}()))
    # no entry: nothing is decided.
    @test isempty(DS.self_fed_nodes([1, 2, 3], act, Set{Int}(), Set{Int}(), Set{Int}(), Set{Int}()))
    # pool transparency: 3 feeds pool 4 which feeds t; Out(3) = {1}.
    act2 = Dict(1 => [2], 2 => [3], 3 => [4], 4 => [1])
    @test DS.self_fed_nodes([1, 2, 3, 4], act2, Set([4]), Set([4]), Set([1]), Set{Int}()) == Set([3])
    # 5 also feeds 3. As an entry (pinned off baseline) it reaches 3 without
    # passing t, so nothing is self-fed; pinned AT baseline it is a constant:
    # not a candidate, not passable, and 3 is self-fed again.
    act3 = Dict(1 => [2], 2 => [3], 3 => [1], 5 => [3])
    @test isempty(DS.self_fed_nodes([1, 2, 3, 5], act3, Set{Int}(), Set{Int}(), Set([1, 5]), Set{Int}()))
    @test DS.self_fed_nodes([1, 2, 3, 5], act3, Set{Int}(), Set{Int}(), Set([1]), Set([5])) == Set([3])
    # an entry that is one of u's own products still counts as passing Out(u):
    # ring with the entry at a (Out(b) = {t}, Out(t) = {a}): t is self-fed, b is not.
    @test DS.self_fed_nodes([1, 2, 3], act, Set{Int}(), Set{Int}(), Set([2]), Set{Int}()) == Set([1])
end

end
