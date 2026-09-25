#!/usr/bin/env julia
#
# Regression tests for configuration guard rails and the API's silent
# network-substitution fallbacks.
#
# Uses `Test` and real `@test` assertions, so a regression makes this file exit
# non-zero. Most other files in this directory wrap everything in try/catch,
# print a failure and return `false` to a caller that discards it, so they
# cannot fail — see test/README or the suite notes.

using Test
using JSON3

include(joinpath(@__DIR__, "..", "src", "DeltaSignal.jl"))
using .DeltaSignal
using .DeltaSignal: resolve_reaction_eval_config, DS_VALID_MODES,
                    _bool_env, _float_env, _mode_env

"""Run `f` with `ENV[name] = value`, restoring the previous state after."""
function with_env(f, name::String, value)
    had = haskey(ENV, name)
    old = had ? ENV[name] : ""
    value === nothing ? delete!(ENV, name) : (ENV[name] = value)
    try
        f()
    finally
        had ? (ENV[name] = old) : delete!(ENV, name)
    end
end

# A top-level `@testset` THROWS when it finishes with any failure, which
# aborts the file and silently skips every testset below it. That is how a
# config mismatch in the dev container hid 61 assertions here -- the cofactor
# tests, the silo-bridge tests and the observation-membership test all looked
# green because they never ran at all.
#
# Nesting them inside one outer testset makes failures accumulate: every
# testset executes, and the outer one reports the full tally at the end. A
# failing run still exits non-zero, it just stops lying about coverage.
@testset "configuration guard rails" begin

@testset "DS_* mode validation" begin
    # Every value the dispatch actually implements must still be accepted.
    # This is the half that matters: an over-tight allowlist breaks working
    # configurations, which is worse than the bug it replaces.
    @testset "all implemented modes accepted" begin
        for (name, values) in DS_VALID_MODES, value in values
            with_env(name, value) do
                @test resolve_reaction_eval_config() isa DeltaSignal.ReactionEvalConfig
            end
        end
    end

    # An unrecognised mode used to fall through to a *different model*. For
    # DS_INHIBITION_MODE the fallback was "spec" with default beta 0, i.e.
    # inhibition silently switched off entirely.
    @testset "typos rejected, not silently remapped" begin
        for bad in ("Divide", "divide ", "dividee", "", "DIVIDE")
            with_env("DS_INHIBITION_MODE", bad) do
                @test_throws ArgumentError resolve_reaction_eval_config()
            end
        end
        for bad in ("Hill_Log", "hill-log", "bogus")
            with_env("DS_AND_MODE", bad) do
                @test_throws ArgumentError resolve_reaction_eval_config()
            end
        end
        with_env("DS_OR_COMBINE", "Gate") do
            @test_throws ArgumentError resolve_reaction_eval_config()
        end
    end

    @testset "error names the variable and lists valid values" begin
        with_env("DS_INHIBITION_MODE", "nope") do
            err = try resolve_reaction_eval_config() catch e; e end
            @test err isa ArgumentError
            @test occursin("DS_INHIBITION_MODE", err.msg)
            @test occursin("divide", err.msg)
        end
    end
end

@testset "DS_* boolean parsing" begin
    # Flags were tested two incompatible ways: some `== "1"`, others `!= "0"`.
    # DS_ASSEMBLY_LIMITING="true" therefore turned the feature OFF while
    # DS_SCC_SOLVE="false" left it ON. Both spellings must now mean the same.
    @testset "true spellings" begin
        for value in ("1", "true", "TRUE", "yes", "on", " true ")
            @test _bool_env("DS_TEST_FLAG", false) == false  # unset -> default
            with_env("DS_TEST_FLAG", value) do
                @test _bool_env("DS_TEST_FLAG", false) == true
            end
        end
    end

    @testset "false spellings" begin
        for value in ("0", "false", "FALSE", "no", "off")
            with_env("DS_TEST_FLAG", value) do
                @test _bool_env("DS_TEST_FLAG", true) == false
            end
        end
    end

    @testset "non-boolean rejected" begin
        for value in ("bogus", "2", "")
            with_env("DS_TEST_FLAG", value) do
                @test_throws ArgumentError _bool_env("DS_TEST_FLAG", false)
            end
        end
    end

    @testset "assembly limiting honours word spellings" begin
        with_env("DS_ASSEMBLY_LIMITING", "true") do
            @test resolve_reaction_eval_config().assembly_limiting == true
        end
        with_env("DS_ASSEMBLY_LIMITING", "false") do
            @test resolve_reaction_eval_config().assembly_limiting == false
        end
    end
end

@testset "DS_* numeric parsing" begin
    @test _float_env("DS_TEST_NUM", 2.5) == 2.5           # unset -> default
    with_env("DS_TEST_NUM", "3.5") do
        @test _float_env("DS_TEST_NUM", 2.5) == 3.5
    end
    for bad in ("abc", "", "NaN", "Inf")
        with_env("DS_TEST_NUM", bad) do
            @test_throws ArgumentError _float_env("DS_TEST_NUM", 2.5)
        end
    end
    with_env("DS_HILL_LOG_ZMAX", "not-a-number") do
        err = try resolve_reaction_eval_config() catch e; e end
        @test err isa ArgumentError
        @test occursin("DS_HILL_LOG_ZMAX", err.msg)
    end
end

@testset "defaults are the validated winning config" begin
    for name in ("DS_INHIBITION_MODE", "DS_AND_MODE", "DS_OR_MODE",
                 "DS_ASSEMBLY_LIMITING", "DS_OR_COMBINE", "DS_INHIBITOR_OR",
                 "DS_SELF_INHIBITOR_WEIGHT",
                 "DS_HILL_SAT_EPS")
        haskey(ENV, name) && delete!(ENV, name)
    end
    config = resolve_reaction_eval_config()
    @test config.inhibition_mode == "divide"
    # hill_sat / clamp-off / eps 1e-5 as of specs/002-upregulation-propagation.
    # test/test_and_curves.jl covers why; this only pins that the defaults are
    # what that feature measured.
    # Re-measured 2026-09-16 on 23,022 wide-curator cases at Release97, one
    # catalog build, each arm changing ONE variable and paired against the
    # others. specs/002 set `hill_sat` and `assembly_limiting=false` on +31
    # over 564 experimental cases; on the curator set those same two choices
    # cost -90 and -196 respectively. The wide set decides (FR-008), so they
    # are reverted here. Attribution and numbers: specs/009-solver-defaults.
    @test config.and_mode == "hill_sat"
    @test config.or_mode == "mean"
    # specs/023 (2026-09-25): off under the root-pinning protocol of record.
    @test config.assembly_limiting == false
    # specs/022, adopted in specs/023.
    @test config.self_inhibitor_weight == 0.1
    # Only consulted when and_mode is hill_sat, so inert at the current
    # default. Kept sized correctly so switching AND mode does not also
    # silently re-introduce a 10%-of-baseline epsilon.
    #
    # 1e-9, lowered from 1e-5. This is the smooth-max floor against zero, so
    # it bounds how far BELOW baseline a value can travel, and 1e-5 against a
    # baseline of 0.01 distorted the low end badly: a fold of 0.001 read 20.7%
    # high, a fold of 0.0001 read 452% high. The sizing rule that was applied
    # to DS_INHIBITOR_EPS applies here for the same reason -- baseline is 0.01,
    # so the epsilon must be orders of magnitude below it, and 1e-5 is only
    # three. At 1e-9 every AND product is exact across the range and the 100x
    # ceiling still holds.
    @test config.hill_sat_eps == 1e-9
    # A divide-by-zero guard, and nothing else. Baseline is 0.01, so this
    # must stay orders of magnitude below it — at the old 1e-3 it was 10% of
    # baseline and silently set the de-repression ceiling, compressed the
    # response curve and shifted maximum suppression.
    # Swept 2026-09-16: 1e-6 and 1e-12 give BIT-IDENTICAL predictions (same
    # 269 changes, same +109/-98), so at or below 1e-6 this is a pure
    # divide-by-zero guard -- exactly what specs/006 FR-002 required and
    # nobody had confirmed. 1e-4 scores 3 cases better in 23,022 but still
    # changes results, i.e. still acts as a model parameter; taking it would
    # be tuning a guard against the evaluation set.
    @test config.inhibitor_eps == 1e-12
    @test config.inhibitor_eps < 0.01 / 1000
    @test config.or_combine == "max"
    @test config.inhibitor_or == false
    @test config.or_redundancy == 1.0
end

@testset "the inhibition epsilon is a guard, not a model parameter" begin
    # The defect this pins: at 1e-3 the epsilon was 10% of baseline and did
    # three things nobody chose. At a true guard size it does one.
    bl = 0.01
    network = DeltaSignal.ReactionNetwork(
        Dict(
            "T" => DeltaSignal.NetworkNode("T", "T", "protein", nothing, "T", bl),
            "I" => DeltaSignal.NetworkNode("I", "I", "protein", nothing, "I", bl),
        ),
        [DeltaSignal.LogicNetworkEdge("I", "T", true, false, 1.0, "regulator")],
        Dict{String, DeltaSignal.SetExpansionMapping}(),
    )
    function target_fold(eps, inhibitor_ui)
        prev = get(ENV, "DS_INHIBITOR_EPS", nothing)
        ENV["DS_INHIBITOR_EPS"] = string(eps)
        try
            params = DeltaSignal.SteadyStateParams(1.0, 0.1, 500, 1e-6, "penalty")
            result = DeltaSignal.solve_steady_state(network, Dict("I" => (inhibitor_ui, 1.0)), params)
            return result.node_activities["T"] / bl
        finally
            prev === nothing ? delete!(ENV, "DS_INHIBITOR_EPS") : (ENV["DS_INHIBITOR_EPS"] = prev)
        end
    end

    # Halving the inhibitor should double its target: h = bl/x exactly.
    # The old default returned 1.833 here — an 8% error from a "guard".
    @test target_fold(1e-12, 0.5) ≈ 2.0 atol=0.01
    @test target_fold(1e-3, 0.5) < 1.9        # the defect, pinned

    # A full knockout is bounded by the per-reaction clamp, not by epsilon,
    # so shrinking epsilon by nine orders of magnitude must not change it.
    @test target_fold(1e-12, 0.0) ≈ target_fold(1e-9, 0.0) atol=1e-6
    @test target_fold(1e-12, 0.0) <= 10.0 + 1e-6
end

@testset "damping is range-checked" begin
    # nv = (1-λ)x + λF(x) is written back unclamped, so λ outside [0,1] left the
    # model's [0,1] domain and threw an uncaught DomainError from inside
    # hill_log; λ = 0 silently never updated the component.
    network = DeltaSignal.ReactionNetwork(
        Dict(
            "A" => DeltaSignal.NetworkNode("A", "A", "protein", nothing, "A", 0.01),
            "B" => DeltaSignal.NetworkNode("B", "B", "protein", nothing, "B", 0.01),
        ),
        [DeltaSignal.LogicNetworkEdge("A", "B", true, true, 1.0, "input")],
        Dict{String, DeltaSignal.SetExpansionMapping}(),
    )
    observations = Dict("A" => (50.0, 1.0))
    params = DeltaSignal.SteadyStateParams(1.0, 0.1, 100, 1e-6, "penalty")

    with_env("DS_SCC_DAMPING", "0.5") do
        @test DeltaSignal.solve_steady_state(network, observations, params) !== nothing
    end
    for bad in ("1.5", "-0.5", "0.0", "3.0")
        with_env("DS_SCC_DAMPING", bad) do
            @test_throws ArgumentError DeltaSignal.solve_steady_state(
                network, observations, params)
        end
    end
    for bad in ("1.0", "-0.1")
        with_env("DS_DAMPING", bad) do
            @test_throws ArgumentError DeltaSignal.solve_steady_state(
                network, observations, params)
        end
    end
end

@testset "DS_SCC_SWEEP" begin
    # The sweep scheme inside a cyclic component. `gauss_seidel` (the default,
    # historical) writes x[t] back mid-sweep, so the visit order -- which comes
    # from Julia Dict/UUID hashing -- selects which fixed point a multi-root
    # component lands in. Relabelling every UUID in the 92-pathway catalog,
    # verified isomorphic, moved 14 of 23,908 curator predictions. `jacobi`
    # evaluates against the state at the sweep's start and commits together,
    # which is order-free. See test/test_solver_determinism.jl.
    #
    # A typo here must throw rather than silently pick a scheme, which is the
    # same guard-rail rule the other DS_* mode knobs follow.
    network = DeltaSignal.ReactionNetwork(
        Dict(
            "A" => DeltaSignal.NetworkNode("A", "A", "protein", nothing, "A", 0.01),
            "B" => DeltaSignal.NetworkNode("B", "B", "protein", nothing, "B", 0.01),
        ),
        [DeltaSignal.LogicNetworkEdge("A", "B", true, true, 1.0, "input")],
        Dict{String, DeltaSignal.SetExpansionMapping}(),
    )
    observations = Dict("A" => (50.0, 1.0))
    params = DeltaSignal.SteadyStateParams(1.0, 0.1, 100, 1e-6, "penalty")

    # Unset must behave exactly as the explicit default.
    ref = with_env("DS_SCC_SWEEP", nothing) do
        DeltaSignal.solve_steady_state(network, observations, params)
    end
    @test ref !== nothing
    gs = with_env("DS_SCC_SWEEP", "gauss_seidel") do
        DeltaSignal.solve_steady_state(network, observations, params)
    end
    @test gs.node_activities == ref.node_activities

    with_env("DS_SCC_SWEEP", "jacobi") do
        @test DeltaSignal.solve_steady_state(network, observations, params) !== nothing
    end

    # Typos, case variants, near-misses and empty all throw.
    for bad in ("gauss-seidel", "Jacobi", "GAUSS_SEIDEL", "jacobi ", "", "1", "true", "seidel")
        with_env("DS_SCC_SWEEP", bad) do
            @test_throws ArgumentError DeltaSignal.solve_steady_state(
                network, observations, params)
        end
    end
end

@testset "DS_COFACTOR_MODE" begin
    # `inert` is the default: +37 cases on 21,450 curator cases. Deleting the
    # nodes instead cost 84, so there is no "drop" mode to select.
    with_env("DS_COFACTOR_MODE", nothing) do
        @test DeltaSignal.cofactor_mode() == "inert"
    end
    for mode in ("propagate", "inert")
        with_env("DS_COFACTOR_MODE", mode) do
            @test DeltaSignal.cofactor_mode() == mode
        end
    end
    # A typo must fail loudly rather than silently selecting a different
    # model — the defect DS #13 was filed for. "drop" is in this list on
    # purpose: it was a real mode that did nothing, and must not come back
    # as a silent no-op.
    for bad in ("Inert", "drop", "none", "off", "")
        with_env("DS_COFACTOR_MODE", bad) do
            @test_throws ArgumentError DeltaSignal.cofactor_mode()
        end
    end

    # The list is keyed off Reactome stable ids, so it must survive a network
    # whose nodes carry none, and must not match on a prefix.
    @testset "cofactor_uuids" begin
        nodes = Dict(
            "atp" => DeltaSignal.NetworkNode("atp", "R-ALL-113592", "unknown",
                                             nothing, "ATP", 0.01),
            "gene" => DeltaSignal.NetworkNode("gene", "R-HSA-69541", "unknown",
                                              nothing, "TP53", 0.01),
            "none" => DeltaSignal.NetworkNode("none", nothing, "unknown",
                                              nothing, "none", 0.01),
        )
        network = DeltaSignal.ReactionNetwork(
            nodes, DeltaSignal.LogicNetworkEdge[],
            Dict{String, DeltaSignal.SetExpansionMapping}())
        found = DeltaSignal.cofactor_uuids(network)
        @test found == Set(["atp"])
    end

    # The list is the feature; pin its shape and its size so an edit that
    # silently widens it fails here rather than in a benchmark six steps later.
    @test "R-ALL-113592" in DeltaSignal.COFACTOR_STIDS      # ATP [cytosol]
    @test length(DeltaSignal.COFACTOR_STIDS) == 253         # stated in the docstring
    @test all(id -> occursin(r"^R-(ALL|HSA)-\d+$", id), DeltaSignal.COFACTOR_STIDS)
end

@testset "an explicit observation beats the cofactor pin" begin
    # Someone measuring an ATP depletion is not making the modelling
    # assumption `inert` encodes. Overwriting their input would be the silent
    # substitution the rest of this file exists to prevent.
    nodes = Dict(
        "atp" => DeltaSignal.NetworkNode("atp", "R-ALL-113592", "unknown",
                                         nothing, "ATP", 0.01),
        "out" => DeltaSignal.NetworkNode("out", "R-HSA-69541", "unknown",
                                         nothing, "OUT", 0.01),
    )
    edges = [DeltaSignal.LogicNetworkEdge("atp", "out", false, true, 1.0, "input")]
    network = DeltaSignal.ReactionNetwork(
        nodes, edges, Dict{String, DeltaSignal.SetExpansionMapping}())
    params = DeltaSignal.SteadyStateParams(1.0, 0.1, 100, 1e-6, "penalty")

    with_env("DS_COFACTOR_MODE", "inert") do
        # Unobserved: pinned at baseline (UI 1.0 -> internal 0.01).
        quiet = DeltaSignal.solve_steady_state(
            network, Dict{String, Tuple{Float64, Float64}}(), params)
        @test isapprox(quiet.node_activities["atp"], 0.01, atol = 1e-9)

        # Observed: the observation stands, pin or no pin.
        depleted = DeltaSignal.solve_steady_state(
            network, Dict("atp" => (0.0, 1.0)), params)
        @test isapprox(depleted.node_activities["atp"], 0.0, atol = 1e-9)
    end
end

@testset "the cofactor list that ships with a bundle wins" begin
    # The generator derives the list from the same release it generated the
    # network from, so it cannot drift from it. The built-in list can, and did.
    mktempdir() do dir
        logic = joinpath(dir, "logic_network.csv")
        write(logic, "source_id,target_id,pos_neg,and_or,edge_type,stoichiometry\n" *
                     "u-atp,u-out,pos,and,input,1\n")
        write(joinpath(dir, "stid_to_uuid_mapping.csv"),
              "uuid,stable_id\nu-atp,R-ALL-113592\nu-out,R-HSA-69541\n")

        # No cofactors.csv: fall back to the built-in list, do NOT treat the
        # bundle as cofactor-free.
        bare = DeltaSignal.parse_complete_network(
            logic, joinpath(dir, "stid_to_uuid_mapping.csv"))
        @test isempty(bare.cofactor_stids)
        @test DeltaSignal.cofactor_uuids(bare) == Set(["u-atp"])

        # A bundle that declares its own list is believed, including when it
        # declares something the built-in list has never heard of.
        write(joinpath(dir, "cofactors.csv"),
              "stable_id,molecule,chebi_id,name,in_network,reactome_release\n" *
              "R-HSA-69541,Invented,1,Invented [cytosol],1,97\n" *
              "R-ALL-113592,ATP,30616,ATP [cytosol],0,97\n")
        bundled = DeltaSignal.parse_complete_network(
            logic, joinpath(dir, "stid_to_uuid_mapping.csv"))
        # in_network=0 rows are listed for completeness but cannot match here.
        @test bundled.cofactor_stids == Set(["R-HSA-69541"])
        @test DeltaSignal.cofactor_uuids(bundled) == Set(["u-out"])
    end

    # A file that is not a cofactor list must fail loudly, not read as empty.
    mktempdir() do dir
        bad = joinpath(dir, "cofactors.csv")
        write(bad, "something_else\n1\n")
        @test_throws ArgumentError DeltaSignal.parse_cofactor_list(bad)
    end
    @test DeltaSignal.parse_cofactor_list(nothing) == Set{String}()
end

@testset "a corrupt cofactor bundle is never silent" begin
    # CSV.jl types the in_network column from its contents, so the same
    # generated file arrives as Int, Float64, Bool or String depending on what
    # else is in it, and a round-trip through another tool can quote it.
    for truthy in (1, 1.0, true, "1", "1.0", "true", "yes", "", missing)
        @test DeltaSignal._in_network_flag(truthy, "x.csv") === true
    end
    for falsy in (0, 0.0, false, "0", "0.0", "false", "no")
        @test DeltaSignal._in_network_flag(falsy, "x.csv") === false
    end
    # An unreadable value used to die with a bare MethodError naming neither
    # the file nor the column.
    err = try
        DeltaSignal._in_network_flag("maybe", "bundle.csv"); nothing
    catch e; e end
    @test err isa ArgumentError
    @test occursin("bundle.csv", err.msg) && occursin("in_network", err.msg)

    # A truncated or partly-written bundle narrows the model silently: the
    # bundle is authoritative, so declaring one cofactor where the network
    # holds several quietly stops treating the rest as cofactors. It must warn.
    mktempdir() do dir
        logic = joinpath(dir, "logic_network.csv")
        write(logic, "source_id,target_id,pos_neg,and_or,edge_type,stoichiometry\n" *
                     "u-atp,u-rxn,pos,and,input,1\nu-h2o,u-rxn,pos,and,input,1\n")
        write(joinpath(dir, "stid_to_uuid_mapping.csv"),
              "uuid,stable_id\nu-atp,R-ALL-113592\nu-h2o,R-ALL-29356\nu-rxn,R-HSA-1\n")
        # Declares ATP but not H2O, which the built-in list does carry.
        write(joinpath(dir, "cofactors.csv"),
              "stable_id,molecule,chebi_id,name,in_network,reactome_release\n" *
              "R-ALL-113592,ATP,30616,ATP [cytosol],1,97\n")
        net = @test_logs (:warn,) match_mode = :any DeltaSignal.parse_complete_network(
            logic, joinpath(dir, "stid_to_uuid_mapping.csv"))
        # The bundle still wins — this is a warning, not an override.
        @test net.cofactor_stids == Set(["R-ALL-113592"])
        @test DeltaSignal.cofactor_uuids(net) == Set(["u-atp"])
    end
end

@testset "US1: a perturbation does not travel through a cofactor" begin
    # The feature's whole purpose, previously untested. The existing fixtures
    # all make the cofactor a ROOT input, which is the one topology where
    # pinning is a guaranteed no-op — so they would pass with the cofactor
    # block deleted. Here ATP sits BETWEEN the perturbation and the readout.
    #
    #     up --> atp --> out
    #
    # Under `propagate` the knockout must reach `out`; under `inert` it must
    # not, because ATP cannot carry it.
    nodes = Dict(
        "up"  => DeltaSignal.NetworkNode("up", "R-HSA-111", "unknown", nothing, "UP", 0.01),
        "atp" => DeltaSignal.NetworkNode("atp", "R-ALL-113592", "unknown", nothing, "ATP", 0.01),
        "out" => DeltaSignal.NetworkNode("out", "R-HSA-222", "unknown", nothing, "OUT", 0.01),
    )
    edges = [
        DeltaSignal.LogicNetworkEdge("up", "atp", false, true, 1.0, "input"),
        DeltaSignal.LogicNetworkEdge("atp", "out", false, true, 1.0, "input"),
    ]
    network = DeltaSignal.ReactionNetwork(
        nodes, edges, Dict{String, DeltaSignal.SetExpansionMapping}())
    params = DeltaSignal.SteadyStateParams(1.0, 0.1, 100, 1e-6, "penalty")
    knockout = Dict("up" => (0.0, 1.0))

    local propagated, pinned
    with_env("DS_COFACTOR_MODE", "propagate") do
        propagated = DeltaSignal.solve_steady_state(network, knockout, params)
    end
    with_env("DS_COFACTOR_MODE", "inert") do
        pinned = DeltaSignal.solve_steady_state(network, knockout, params)
    end

    # The knockout reaches ATP, and through it the readout, when propagating.
    @test propagated.node_activities["atp"] < 0.01
    @test propagated.node_activities["out"] < 0.01

    # Pinned, ATP holds at baseline and the readout never moves.
    @test isapprox(pinned.node_activities["atp"], 0.01, atol = 1e-9)
    @test isapprox(pinned.node_activities["out"], 0.01, atol = 1e-6)

    # And the two modes genuinely disagree — the guard against a fixture where
    # both arms happen to give the same answer.
    @test pinned.node_activities["out"] > propagated.node_activities["out"]
end

@testset "DS_SILO_BRIDGE_MAX_REACH" begin
    with_env("DS_SILO_BRIDGE_MAX_REACH", nothing) do
        @test DeltaSignal.silo_bridge_max_reach() == 0   # off by default
    end
    for good in ("0", "50", "1000")
        with_env("DS_SILO_BRIDGE_MAX_REACH", good) do
            @test DeltaSignal.silo_bridge_max_reach() == parse(Int, good)
        end
    end
    for bad in ("-1", "50.5", "lots", "")
        with_env("DS_SILO_BRIDGE_MAX_REACH", bad) do
            @test_throws ArgumentError DeltaSignal.silo_bridge_max_reach()
        end
    end

    # sink --(no route)--> source, both the same curated entity.
    #   up -> sink        (sink receives signal, has no outgoing edge)
    #   source -> out     (source feeds a reaction, has no incoming edge)
    nodes = Dict(
        "up"     => DeltaSignal.NetworkNode("up", "R-HSA-1", "unknown", nothing, "UP", 0.01),
        "sink"   => DeltaSignal.NetworkNode("sink", "R-HSA-SPLIT", "unknown", nothing, "S", 0.01),
        "source" => DeltaSignal.NetworkNode("source", "R-HSA-SPLIT", "unknown", nothing, "S", 0.01),
        "out"    => DeltaSignal.NetworkNode("out", "R-HSA-2", "unknown", nothing, "OUT", 0.01),
    )
    edges = [
        DeltaSignal.LogicNetworkEdge("up", "sink", false, true, 1.0, "input"),
        DeltaSignal.LogicNetworkEdge("source", "out", false, true, 1.0, "input"),
    ]
    net = DeltaSignal.ReactionNetwork(
        nodes, edges, Dict{String, DeltaSignal.SetExpansionMapping}())

    @test isempty(DeltaSignal.silo_bridge_edges(net, 0))      # off
    made = DeltaSignal.silo_bridge_edges(net, 50)
    @test length(made) == 1
    @test made[1].parent_uuid == "sink" && made[1].child_uuid == "source"
    @test made[1].edge_type == "silo_bridge"
    @test !made[1].is_and    # must not impose AND-completeness on the target

    # The cap is the whole idea: a bridge reaching more than the cap is refused.
    # This one reaches `source` and `out`, a gain of 2.
    @test length(DeltaSignal.silo_bridge_edges(net, 2)) == 1
    @test isempty(DeltaSignal.silo_bridge_edges(net, 1))

    # An entity split into two uuids that ALREADY have a route is not bridged.
    linked = DeltaSignal.ReactionNetwork(
        nodes,
        vcat(edges, DeltaSignal.LogicNetworkEdge("sink", "source", false, true, 1.0, "input")),
        Dict{String, DeltaSignal.SetExpansionMapping}())
    @test isempty(DeltaSignal.silo_bridge_edges(linked, 50))
end

@testset "the solver silently drops observations for nodes it does not have" begin
    # This is WHY src/api/server.jl checks observation membership. The solver's
    # pinning loop does `haskey(uuid_to_idx, uuid) || continue`, so an
    # observation naming a uuid the network does not contain is discarded
    # without an error: the caller asked for a perturbation, got HTTP 200, and
    # received the unperturbed baseline. The usual cause is solving against a
    # different catalog build than the uuids came from — uuids are regenerated
    # per build, and two builds of one pathway share none.
    #
    # Pinned here so that if the solver ever starts rejecting unknown
    # observations itself, whoever makes that change learns the API-layer guard
    # has become redundant rather than leaving two checks to drift.
    N(id) = DeltaSignal.NetworkNode(id, id, "protein", nothing, id, 0.01)
    net = DeltaSignal.ReactionNetwork(
        Dict(id => N(id) for id in ("A", "r", "B")),
        [DeltaSignal.LogicNetworkEdge("A", "r", true, true, 1.0, "input"),
         DeltaSignal.LogicNetworkEdge("r", "B", false, true, 1.0, "output")],
        Dict{String, DeltaSignal.SetExpansionMapping}())
    params = DeltaSignal.SteadyStateParams(1.0, 0.1, 500, 1e-8, "penalty")

    base = DeltaSignal.solve_steady_state(
        net, Dict{String, Tuple{Float64, Float64}}(), params)
    # A real knockout moves the readout; this is the control.
    real_ko = DeltaSignal.solve_steady_state(net, Dict("A" => (0.0, 1.0)), params)
    @test real_ko.node_activities["B"] < base.node_activities["B"] / 10

    # An observation for a uuid that is not in the network: no error, no effect.
    ghost = DeltaSignal.solve_steady_state(
        net, Dict("not-a-node-in-this-network" => (0.0, 1.0)), params)
    @test ghost.node_activities["B"] == base.node_activities["B"]
    @test ghost.converged == base.converged
end



@testset "a full network payload round-trips through /api/solve" begin
    # The bug this pins: reaction_network_from_json takes `data`, but the
    # cofactor branch referenced an undefined `network_json`. Julia only catches
    # that at runtime, so /api/solve returned 500 for EVERY request that sent a
    # network payload instead of a cached network_id — the path an upload uses,
    # and the path the benchmark's DS_SKIP_EDGE_TYPES diagnostic uses. Both were
    # silently broken; the diagnostic reported 0/0 scored cases, which is how it
    # was found.
    #
    # Reached directly rather than over HTTP: this repo has no API test harness,
    # and the defect is in the JSON→network conversion, not in the routing.
    payload = """
    {"nodes": [{"uuid": "A", "name": "A", "reactome_id": "R-HSA-1",
                "entity_type": "protein", "baseline": 0.01, "set_id": null},
               {"uuid": "B", "name": "B", "reactome_id": "R-HSA-2",
                "entity_type": "protein", "baseline": 0.01, "set_id": null}],
     "edges": [{"parent_uuid": "A", "child_uuid": "B", "is_and": true,
                "is_positive": true, "stoichiometry": 1.0, "edge_type": "input"}]}
    """
    net = DeltaSignal.reaction_network_from_json(JSON3.read(payload))
    @test length(net.nodes) == 2
    @test length(net.edges) == 1
    @test isempty(net.cofactor_stids)

    # And the round-trip the code comment promises: a bundle's cofactor list
    # must survive parse → solve rather than being silently dropped.
    with_cof = replace(payload, "\"edges\":" =>
        "\"cofactor_stids\": [\"R-ALL-29372\", \"R-ALL-113582\"], \"edges\":")
    net2 = DeltaSignal.reaction_network_from_json(JSON3.read(with_cof))
    @test net2.cofactor_stids == Set(["R-ALL-29372", "R-ALL-113582"])
    @test length(net2.nodes) == 2
end

@testset "edge sign and logic values fail loudly" begin
    # `is_positive = pos_raw == "pos"` made every UNRECOGNISED sign an
    # inhibitor, silently inverting the edge. That is the same class as the
    # DS_* mode typos guarded above: a wrong value selected different
    # behaviour with no error. Only pos/neg and and/or/"" occur in the catalog
    # today -- all 92 networks and 334,844 edges parse unchanged -- so this
    # guards a future typo rather than fixing current data.
    #
    # Goes through `parse_logic_network`, the public entry point, so the
    # schema detection is exercised too.
    function gen_edge(pn, ao)
        path = tempname() * ".csv"
        write(path, "source_id,target_id,pos_neg,and_or,stoichiometry,edge_type\n" *
                    "A,B,$(pn),$(ao),1.0,input\n")
        try
            return DeltaSignal.parse_logic_network(path)[1]
        finally
            rm(path; force=true)
        end
    end
    function sample_edge(a, p)
        path = tempname() * ".csv"
        write(path, "parent,child,is_and,is_positive,stoichiometry\nA,B,$(a),$(p),1.0\n")
        try
            return DeltaSignal.parse_logic_network(path)[1]
        finally
            rm(path; force=true)
        end
    end

    # Valid values are unchanged, including case and the documented empty
    # and_or meaning OR.
    @test gen_edge("pos", "and").is_positive
    @test gen_edge("pos", "and").is_and
    @test !gen_edge("neg", "or").is_positive
    @test !gen_edge("neg", "or").is_and
    @test gen_edge("pos", "").is_positive
    @test !gen_edge("pos", "").is_and
    @test gen_edge("POS", "AND").is_positive
    @test gen_edge("POS", "AND").is_and

    # A typo throws instead of inverting the edge.
    @test_throws ArgumentError gen_edge("positive", "and")
    @test_throws ArgumentError gen_edge("negative", "and")
    @test_throws ArgumentError gen_edge("pos", "xor")

    # The message names the offending value so the row can be found.
    err = try
        gen_edge("positive", "and")
        ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("positive", err)
    @test occursin("pos_neg", err)

    # Sample format carries the same guard: anything but 0/1 is a mistake.
    @test sample_edge(1, 1).is_and
    @test sample_edge(1, 1).is_positive
    @test !sample_edge(0, 0).is_and
    @test !sample_edge(0, 0).is_positive
    @test_throws ArgumentError sample_edge(2, 1)
    @test_throws ArgumentError sample_edge(1, -1)
end

end  # outer testset
