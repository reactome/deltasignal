#!/usr/bin/env julia

"""Benchmark DeltaSignal solver invariants across historical and current configurations.

This is a solver-only benchmark: every configuration receives the same synthetic
ReactionNetwork objects. It separates propagation-math changes from SCC loop
handling and complements, rather than replaces, empirical perturbation benchmarks.
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV
using DataFrames
using JSON3
using DeltaSignal

const BASELINE_UI = 1.0

const CONFIGS = [
    (
        "legacy_math_flat",
        Dict(
            "DS_SCC_SOLVE" => "0",
            "DS_INHIBITION_MODE" => "spec",
            "DS_AND_MODE" => "geomean",
            "DS_OR_MODE" => "max",
            "DS_ASSEMBLY_LIMITING" => "0",
            "DS_HILL_LOG_ZMAX" => "4.605170185988092",
        ),
    ),
    (
        "legacy_math_scc",
        Dict(
            "DS_SCC_SOLVE" => "1",
            "DS_INHIBITION_MODE" => "spec",
            "DS_AND_MODE" => "geomean",
            "DS_OR_MODE" => "max",
            "DS_ASSEMBLY_LIMITING" => "0",
            "DS_HILL_LOG_ZMAX" => "4.605170185988092",
        ),
    ),
    (
        "current_math_flat",
        Dict(
            "DS_SCC_SOLVE" => "0",
            "DS_INHIBITION_MODE" => "divide",
            "DS_AND_MODE" => "hill_log",
            "DS_OR_MODE" => "mean",
            "DS_ASSEMBLY_LIMITING" => "1",
            "DS_HILL_LOG_ZMAX" => "10.0",
        ),
    ),
    (
        "current_math_scc",
        Dict(
            "DS_SCC_SOLVE" => "1",
            "DS_INHIBITION_MODE" => "divide",
            "DS_AND_MODE" => "hill_log",
            "DS_OR_MODE" => "mean",
            "DS_ASSEMBLY_LIMITING" => "1",
            "DS_HILL_LOG_ZMAX" => "10.0",
        ),
    ),
]


function parse_args(args)
    output_dir = joinpath(@__DIR__, "results", "solver_invariants")
    i = 1
    while i <= length(args)
        args[i] == "--output-dir" || error("Unknown argument: $(args[i])")
        i == length(args) && error("--output-dir requires a value")
        output_dir = args[i + 1]
        i += 2
    end
    return abspath(output_dir)
end


node(id::String) = NetworkNode(id, nothing, "synthetic", nothing, id, 0.01)


function network(node_ids::Vector{String}, edge_specs)
    nodes = Dict(id => node(id) for id in node_ids)
    edges = LogicNetworkEdge[
        LogicNetworkEdge(source, target, is_and, is_positive, stoich, edge_type)
        for (source, target, is_and, is_positive, stoich, edge_type) in edge_specs
    ]
    return ReactionNetwork(nodes, edges, Dict{String, SetExpansionMapping}())
end


function solve_ui(net::ReactionNetwork, observations; max_iters=500, tolerance=1e-8)
    params = SteadyStateParams(1.0, 0.1, max_iters, tolerance, "penalty")
    result = solve_steady_state(net, observations, params)
    activities = Dict(id => 100.0 * value for (id, value) in result.node_activities)
    return result, activities
end


function config_pairs(config::Dict{String, String})
    controlled = Dict{String, Union{String, Nothing}}(
        "DS_DAMPING" => "0.0",
        "DS_DROP_PASSTHROUGH" => "0",
        "DS_INHIBITOR_BETA" => nothing,
        "DS_INHIBITOR_FLOOR" => "0.0",
        "DS_INHIBITOR_FLOOR_SCOPE" => "loops",
        "DS_SCC_BREAK_CATALYST" => "0",
        "DS_SCC_DAMPING" => "0.5",
        "DS_SCC_NEG_MODE" => "converge",
    )
    merge!(controlled, config)
    return Pair{String, Union{String, Nothing}}[key => value for (key, value) in controlled]
end


function add_result!(rows, config, test, category, passed, value, expected, details)
    push!(rows, (
        config=config,
        test=test,
        category=category,
        passed=passed,
        value=Float64(value),
        expected=expected,
        details=details,
    ))
end


function run_config(config_name::String, config::Dict{String, String})
    rows = NamedTuple[]
    withenv(config_pairs(config)...) do
        chain = network(
            ["A", "B", "C"],
            [
                ("A", "B", false, true, 1.0, "output"),
                ("B", "C", false, true, 1.0, "output"),
            ],
        )
        up_result, up = solve_ui(chain, Dict("A" => (80.0, 1.0)))
        down_result, down = solve_ui(chain, Dict("A" => (0.0, 1.0)))
        add_result!(rows, config_name, "positive_chain_up", "causal_polarity",
                    up["C"] > BASELINE_UI, up["C"], "> 1", "A=80 should increase C")
        add_result!(rows, config_name, "positive_chain_down", "causal_polarity",
                    down["C"] < BASELINE_UI, down["C"], "< 1", "A=0 should decrease C")
        add_result!(rows, config_name, "positive_chain_converges", "numerics",
                    up_result.converged && down_result.converged, max(up_result.final_residual, down_result.final_residual),
                    "converged", "Acyclic chain must solve exactly")

        inhibition = network(
            ["I", "B"],
            [("I", "B", false, false, 1.0, "regulator")],
        )
        _, inhibitor_down = solve_ui(inhibition, Dict("I" => (0.0, 1.0)))
        _, inhibitor_up = solve_ui(inhibition, Dict("I" => (80.0, 1.0)))
        add_result!(rows, config_name, "inhibitor_knockout_derepresses", "causal_polarity",
                    inhibitor_down["B"] > BASELINE_UI, inhibitor_down["B"], "> 1",
                    "Removing an inhibitor should increase its target")
        add_result!(rows, config_name, "inhibitor_up_suppresses", "causal_polarity",
                    inhibitor_up["B"] < BASELINE_UI, inhibitor_up["B"], "< 1",
                    "Increasing an inhibitor should decrease its target")

        assembly = network(
            ["S1", "S2", "COMPLEX"],
            [
                ("S1", "COMPLEX", true, true, 1.0, "assembly"),
                ("S2", "COMPLEX", true, true, 1.0, "assembly"),
            ],
        )
        _, one_subunit = solve_ui(
            assembly,
            Dict("S1" => (80.0, 1.0), "S2" => (BASELINE_UI, 1.0)),
        )
        _, both_subunits = solve_ui(
            assembly,
            Dict("S1" => (80.0, 1.0), "S2" => (80.0, 1.0)),
        )
        add_result!(rows, config_name, "assembly_single_subunit_capped", "complex_assembly",
                    one_subunit["COMPLEX"] <= 1.05, one_subunit["COMPLEX"], "<= 1.05",
                    "One abundant subunit cannot create an abundant complex")
        add_result!(rows, config_name, "assembly_all_subunits_propagate", "complex_assembly",
                    both_subunits["COMPLEX"] > 10.0, both_subunits["COMPLEX"], "> 10",
                    "All abundant subunits should increase the complex")

        disconnected = network(
            ["A", "B", "X", "Y"],
            [
                ("A", "B", false, true, 1.0, "output"),
                ("X", "Y", false, true, 1.0, "output"),
            ],
        )
        _, locality_outputs = solve_ui(disconnected, Dict("A" => (80.0, 1.0)))
        add_result!(rows, config_name, "disconnected_component_locality", "causal_locality",
                    isapprox(locality_outputs["Y"], BASELINE_UI; atol=1e-8), locality_outputs["Y"], "1.0",
                    "Perturbing A must not alter disconnected Y")

        permuted = network(
            ["A", "B", "C"],
            [
                ("B", "C", false, true, 1.0, "output"),
                ("A", "B", false, true, 1.0, "output"),
            ],
        )
        _, permuted_up = solve_ui(permuted, Dict("A" => (80.0, 1.0)))
        add_result!(rows, config_name, "edge_order_invariance", "determinism",
                    isapprox(permuted_up["C"], up["C"]; atol=1e-10),
                    abs(permuted_up["C"] - up["C"]), "absolute delta <= 1e-10",
                    "Reordering edges must not change a solve")

        negative_loop = network(
            ["S", "A", "B"],
            [
                ("S", "A", false, true, 1.0, "input"),
                ("B", "A", false, false, 1.0, "regulator"),
                ("A", "B", false, true, 1.0, "output"),
            ],
        )
        negative_result, negative = solve_ui(negative_loop, Dict("S" => (80.0, 1.0)))
        add_result!(rows, config_name, "negative_feedback_converges", "feedback_loop",
                    negative_result.converged, negative_result.final_residual, "residual < 1e-8",
                    "A<->B negative-feedback SCC should converge")
        add_result!(rows, config_name, "negative_feedback_bounded", "feedback_loop",
                    0.0 <= negative["B"] <= 100.0, negative["B"], "0..100",
                    "Negative feedback output must remain bounded")

        positive_loop = network(
            ["S", "A", "B"],
            [
                ("S", "A", false, true, 1.0, "input"),
                ("B", "A", false, true, 1.0, "input"),
                ("A", "B", false, true, 1.0, "output"),
            ],
        )
        positive_result, positive = solve_ui(positive_loop, Dict("S" => (80.0, 1.0)))
        add_result!(rows, config_name, "positive_feedback_converges", "feedback_loop",
                    positive_result.converged, positive_result.final_residual, "residual < 1e-8",
                    "A<->B positive-feedback SCC should converge")
        add_result!(rows, config_name, "positive_feedback_bounded", "feedback_loop",
                    0.0 <= positive["B"] <= 100.0, positive["B"], "0..100",
                    "Positive feedback output must remain bounded")

        pinned_result, pinned = solve_ui(chain, Dict("A" => (37.0, 1.0), "B" => (12.0, 1.0)))
        pin_error = max(abs(pinned["A"] - 37.0), abs(pinned["B"] - 12.0))
        add_result!(rows, config_name, "observations_are_hard_pinned", "observation_contract",
                    pin_error <= 1e-10 && pinned_result.converged, pin_error, "max error <= 1e-10",
                    "Observed nodes must equal supplied values")
    end
    return rows
end


function main()
    output_dir = parse_args(ARGS)
    mkpath(output_dir)
    rows = NamedTuple[]
    for (name, config) in CONFIGS
        append!(rows, run_config(name, config))
    end
    table = DataFrame(rows)
    CSV.write(joinpath(output_dir, "solver_invariant_results.tsv"), table; delim='\t')

    summary = combine(
        groupby(table, :config),
        nrow => :tests,
        :passed => sum => :passed,
    )
    summary.failed = summary.tests .- summary.passed
    CSV.write(joinpath(output_dir, "solver_invariant_summary.tsv"), summary; delim='\t')

    manifest = Dict(
        "deltasignal_commit" => try
            readchomp(`git -C $(joinpath(@__DIR__, "..")) rev-parse HEAD`)
        catch
            nothing
        end,
        "configurations" => Dict(name => config for (name, config) in CONFIGS),
        "summary" => [Dict(String(k) => v for (k, v) in pairs(row)) for row in eachrow(summary)],
    )
    open(joinpath(output_dir, "solver_invariant_manifest.json"), "w") do handle
        JSON3.pretty(handle, manifest)
    end

    println(summary)
    any(.!table.passed) && println("See solver_invariant_results.tsv for failed invariants.")
end


main()
