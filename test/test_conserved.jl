# specs/035: loop participants released unchanged by their loop (a conserved
# pool, like a cofactor) held at baseline under DS_CONSERVED_MODE=inert.
#
# Fixture = the RAF shape. Loop: F (e.g. F-actin) and the signal S enter R1
# (AND) -> complex C (contains F) -> R2 (a release: consumes C) -> F and S2 ->
# back. F is released unchanged: CONSERVED. S2 is produced by R2 but C does not
# contain it: transformed, NOT conserved. G is a loop input with an outside
# producer (Z -> G): NOT conserved.
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

function fixture()
    nodes = Dict(node.(["X", "F", "S2", "G", "Z", "R1", "C", "R2", "T"]))
    edges = [edge("X", "R1", "input"), edge("F", "R1", "input"), edge("S2", "R1", "input"),
             edge("G", "R1", "input"),
             edge("R1", "C", "output"; and = false),
             edge("C", "R2", "input"),
             edge("R2", "F", "output"; and = false), edge("R2", "S2", "output"; and = false),
             edge("R2", "G", "output"; and = false), edge("Z", "G", "output"; and = false),
             edge("C", "T", "input")]
    cont = Dict("R-C" => Set(["R-F", "R-G"]))        # C holds F and G, not S2
    DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(), cont)
end

solve(mode; obs = Dict("X" => (2.0, 1.0), "Z" => (1.0, 1.0))) =
    with_env("DS_CONSERVED_MODE" => mode, "DS_SELF_INHIBITOR_WEIGHT" => "off") do
        DS.solve_steady_state(fixture(), obs, DS.SteadyStateParams(1.0, 0.1, 500, 1e-6, "penalty"))
    end

@testset "conserved participants (specs/035)" begin

@testset "detection: released unchanged, all supply in the loop" begin
    @test DS.conserved_uuids(fixture()) == Set(["F"])
    # without the containment table nothing is known to be released unchanged
    f = fixture()
    bare = DS.ReactionNetwork(f.nodes, f.edges, f.set_mappings)
    @test isempty(DS.conserved_uuids(bare))
end

@testset "detection is label-independent" begin
    f = fixture()
    rl(u) = "q_" * u
    nodes = Dict(rl(k) => DS.NetworkNode(rl(k), v.reactome_id, v.entity_type, nothing, rl(k), v.baseline)
                 for (k, v) in f.nodes)
    edges = [DS.LogicNetworkEdge(rl(e.parent_uuid), rl(e.child_uuid), e.is_and, e.is_positive,
                                 e.stoichiometry, e.edge_type) for e in f.edges]
    net = DS.ReactionNetwork(nodes, edges, f.set_mappings, Set{String}(), f.containment)
    @test DS.conserved_uuids(net) == Set(["q_F"])
end

@testset "default off, byte-identical and reported" begin
    with_env("DS_CONSERVED_MODE" => nothing) do
        @test DS.conserved_mode() == "off"
    end
    r = solve(nothing)
    @test r.diagnostics["conserved_rule"] == "off"
    @test r.diagnostics["conserved_held"] == 0
end

@testset "inert holds F at baseline; an observation of F wins" begin
    r = solve("inert")
    @test r.diagnostics["conserved_rule"] == "inert"
    @test r.diagnostics["conserved_held"] == 1
    @test r.node_activities["F"] / BL ≈ 1.0 rtol = 1e-9
    r = solve("inert"; obs = Dict("X" => (2.0, 1.0), "Z" => (1.0, 1.0), "F" => (0.0, 1.0)))
    @test r.node_activities["F"] / BL ≈ 0.0 atol = 1e-12
    @test r.diagnostics["conserved_held"] == 0
    r = solve("inert"; obs = Dict("X" => (2.0, 1.0), "Z" => (1.0, 1.0), "F" => (0.0, 0.0)))
    @test r.node_activities["F"] / BL ≈ 1.0 rtol = 1e-9       # a gated-out row does not unpin
end

@testset "a typo is an error" begin
    for bad in ("Inert", "on", "")
        with_env("DS_CONSERVED_MODE" => bad) do
            @test_throws ArgumentError DS.conserved_mode()
        end
    end
end

end
