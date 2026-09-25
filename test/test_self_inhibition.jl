# specs/022: an inhibitor that contains its own reaction's input.
#
# Fixture = the WNT shape traced in specs/021: target T has AND inputs L (the
# ligand) and F (a partner at baseline), and OR-flagged inhibitors C1, C2 whose
# containment includes L. C1 and C2 are pinned to L's fold, which is what the
# benchmark's set pinning produced. C3 is an unrelated inhibitor.
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

function fixture(; containment = Dict("R-C1" => Set(["R-L", "R-X"]), "R-C2" => Set(["R-L"])),
                 extra_inhibitor = false)
    nodes = Dict(node.(["T", "L", "F", "C1", "C2", "C3"]))
    edges = [DS.LogicNetworkEdge("L", "T", true, true, 1.0, "input"),
             DS.LogicNetworkEdge("F", "T", true, true, 1.0, "input"),
             DS.LogicNetworkEdge("C1", "T", false, false, 1.0, "regulator"),
             DS.LogicNetworkEdge("C2", "T", false, false, 1.0, "regulator")]
    extra_inhibitor && push!(edges, DS.LogicNetworkEdge("C3", "T", false, false, 1.0, "regulator"))
    DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(), containment)
end

"Fold of T with L, C1, C2 pinned at `x` (UI = fold) and C3 at `c3`."
function t_fold(net; x, c1 = x, c2 = x, c3 = 1.0, w = "off")
    with_env("DS_SELF_INHIBITOR_WEIGHT" => w) do
        obs = Dict("L" => (x, 1.0), "F" => (1.0, 1.0), "C1" => (c1, 1.0),
                   "C2" => (c2, 1.0), "C3" => (c3, 1.0))
        r = DS.solve_steady_state(net, obs, DS.SteadyStateParams(1.0, 0.1, 500, 1e-6, "penalty"))
        (r.node_activities["T"] / BL, get(r.diagnostics, "self_inhibitors", -1))
    end
end

@testset "self-contained inhibition (specs/022)" begin

@testset "default is 0.1 (specs/023); off reproduces the double count exactly" begin
    with_env("DS_SELF_INHIBITOR_WEIGHT" => nothing) do
        @test DS.resolve_reaction_eval_config().self_inhibitor_weight == 0.1
    end
    with_env("DS_SELF_INHIBITOR_WEIGHT" => "off") do
        @test DS.resolve_reaction_eval_config().self_inhibitor_weight == -1.0
    end
    @test t_fold(fixture(); x = 0.5, w = nothing)[1] ≈ 0.5^0.9 rtol = 1e-6   # unset = the default
    net = fixture()
    for (x, want) in ((0.5, 2.0), (0.2, 2.0), (0.05, 0.5), (2.0, 0.5))
        got, n = t_fold(net; x = x)
        @test got ≈ want rtol = 1e-6        # x * min(x^-2, 10): the specs/021 inversion
        @test n == 0                         # nothing damped when off
    end
end

@testset "with w: the input sets the direction, the inhibitors only damp it" begin
    net = fixture()
    for x in (0.05, 0.2, 0.5, 2.0, 5.0, 20.0)
        got, n = t_fold(net; x = x, w = "0.1")
        @test got ≈ x^0.9 rtol = 1e-6        # x^(1-w), however many inhibitors
        @test n == 2
    end
    ups = [t_fold(net; x = x, w = "0.1")[1] for x in (1.0, 2.0, 5.0, 20.0)]
    downs = [t_fold(net; x = x, w = "0.1")[1] for x in (1.0, 0.5, 0.2, 0.05)]
    @test issorted(ups) && issorted(downs; rev = true)       # monotone both ways
    @test t_fold(net; x = 0.5, w = "0")[1] ≈ 0.5 rtol = 1e-6   # w = 0: input alone
    @test t_fold(net; x = 0.5, w = "1")[1] ≈ 1.0 rtol = 1e-6   # w = 1: the old single-count cancel
    @test t_fold(net; x = 1.0, w = "0.1")[1] ≈ 1.0 rtol = 1e-9 # baseline untouched
end

@testset "the inhibitor's independent part keeps full strength" begin
    net = fixture()
    # C1 doubles for its own reason (WIF1 overexpressed) while L is at baseline.
    @test t_fold(net; x = 1.0, c1 = 2.0, c2 = 1.0, w = "0.1")[1] ≈ 0.5 rtol = 1e-6
    @test t_fold(net; x = 1.0, c1 = 2.0, c2 = 1.0)[1] ≈ 0.5 rtol = 1e-6
    # An inhibitor that contains nothing of T's inputs is not touched.
    net3 = fixture(extra_inhibitor = true)
    @test t_fold(net3; x = 1.0, c1 = 1.0, c2 = 1.0, c3 = 2.0, w = "0.1")[1] ≈ 0.5 rtol = 1e-6
    @test t_fold(net3; x = 1.0, c3 = 2.0, w = "0.1")[2] == 2
end

@testset "inert without a containment table, and says so" begin
    net = fixture(containment = Dict{String, Set{String}}())
    got, n = @test_logs (:warn, r"no containment table") match_mode = :any t_fold(net; x = 0.5, w = "0.1")
    @test got ≈ 2.0 rtol = 1e-6
    @test n == 0
end

@testset "invalid weights are startup errors" begin
    for bad in ("1.5", "-0.1", "abc", "NaN")
        with_env("DS_SELF_INHIBITOR_WEIGHT" => bad) do
            @test_throws ArgumentError DS.resolve_reaction_eval_config()
        end
    end
    with_env("DS_SELF_INHIBITOR_WEIGHT" => "0.1") do
        @test DS.resolve_reaction_eval_config().self_inhibitor_weight == 0.1
    end
end

@testset "containment.csv is read beside the logic network, self rows dropped" begin
    mktempdir() do d
        write(joinpath(d, "containment.csv"),
              "stable_id,contains_stable_id,reactome_release\nR-C1,R-C1,97\nR-C1,R-L,97\nR-C2,R-L,97\n")
        c = DS.parse_containment(joinpath(d, "logic_network.csv"))
        @test c == Dict("R-C1" => Set(["R-L"]), "R-C2" => Set(["R-L"]))
        @test isempty(DS.parse_containment(joinpath(d, "sub", "logic_network.csv")))
    end
end

@testset "the map matches per target node, never a sibling's input" begin
    nodes = Dict(node.(["T1", "T2", "L1", "L2", "C"]))
    edges = [DS.LogicNetworkEdge("L1", "T1", true, true, 1.0, "input"),
             DS.LogicNetworkEdge("L2", "T2", true, true, 1.0, "input"),
             DS.LogicNetworkEdge("C", "T1", false, false, 1.0, "regulator"),
             DS.LogicNetworkEdge("C", "T2", false, false, 1.0, "regulator"),
             DS.LogicNetworkEdge("C", "L1", false, false, 1.0, "depletion")]
    net = DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                             Dict("R-C" => Set(["R-L1"])))
    @test DS.self_contained_inhibitor_map(net) == Dict(("C", "T1") => ["L1"])
end

end
