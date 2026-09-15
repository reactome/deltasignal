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
                 "DS_HILL_SAT_EPS")
        haskey(ENV, name) && delete!(ENV, name)
    end
    config = resolve_reaction_eval_config()
    @test config.inhibition_mode == "divide"
    # hill_sat / clamp-off / eps 1e-5 as of specs/002-upregulation-propagation.
    # test/test_and_curves.jl covers why; this only pins that the defaults are
    # what that feature measured.
    @test config.and_mode == "hill_sat"
    @test config.or_mode == "mean"
    @test config.assembly_limiting == false
    @test config.hill_sat_eps == 1e-5
    # A divide-by-zero guard, and nothing else. Baseline is 0.01, so this
    # must stay orders of magnitude below it — at the old 1e-3 it was 10% of
    # baseline and silently set the de-repression ceiling, compressed the
    # response curve and shifted maximum suppression.
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
