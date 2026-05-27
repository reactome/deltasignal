#!/usr/bin/env julia

# Comprehensive biological realism validation tests
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Random
using Statistics
using CSV
using DataFrames

# Add the src directory to the load path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
include("../src/io/tsv_parser.jl")
include("../src/core/sensitivity.jl")
include("../src/core/aggregators.jl") 
include("../src/core/hill_functions.jl")
include("../src/core/reaction_model.jl")
include("../src/core/feedback_enhancements.jl")
include("../src/solvers/steady_state.jl")
include("../src/solvers/enhanced_steady_state.jl")
include("../src/solvers/time_dynamic.jl")

function test_biological_realism()
    println("🧬 Testing Biological Realism & Expected Behaviors")
    println("=" ^ 60)
    
    try
        # Load the test network
        logic_path = joinpath(@__DIR__, "..", "examples", "propagation_test_network.tsv")
        uuid_path = joinpath(@__DIR__, "..", "examples", "propagation_test_uuid_mapping.tsv")
        
        network = parse_complete_network(logic_path, uuid_path, nothing)
        
        println("📁 Network loaded: $(length(network.nodes)) nodes, $(length(network.edges)) edges")
        
        # Biological validation tests
        test_results = Dict{String, Bool}()
        
        # Test 1: Dose-Response Curves (Monotonicity)
        println("\\n" * "="^60)
        println("TEST 1: Dose-Response Monotonicity")
        println("="^60)
        
        test_results["dose_response"] = test_dose_response_monotonicity(network)
        
        # Test 2: Signal Attenuation with Distance
        println("\\n" * "="^60)
        println("TEST 2: Signal Attenuation with Distance")
        println("="^60)
        
        test_results["signal_attenuation"] = test_signal_attenuation(network)
        
        # Test 3: Competitive Inhibition
        println("\\n" * "="^60)
        println("TEST 3: Competitive Inhibition Behavior")
        println("="^60)
        
        test_results["competitive_inhibition"] = test_competitive_inhibition(network)
        
        # Test 4: Feedback Loop Behavior
        println("\\n" * "="^60)
        println("TEST 4: Feedback Loop Biological Behavior")
        println("="^60)
        
        test_results["feedback_loops"] = test_feedback_loop_behavior(network)
        
        # Test 5: Dynamic Response Patterns
        println("\\n" * "="^60)
        println("TEST 5: Time-Dynamic Response Patterns")
        println("="^60)
        
        test_results["dynamic_patterns"] = test_dynamic_response_patterns(network)
        
        # Test 6: Pathway Crosstalk
        println("\\n" * "="^60)
        println("TEST 6: Pathway Crosstalk & Independence")
        println("="^60)
        
        test_results["pathway_crosstalk"] = test_pathway_crosstalk(network)
        
        # Test 7: Saturation & Hill Kinetics
        println("\\n" * "="^60)
        println("TEST 7: Saturation & Hill Kinetics Behavior")
        println("="^60)
        
        test_results["hill_kinetics"] = test_hill_kinetics_behavior(network)
        
        # Test 8: Biological Range Validation
        println("\\n" * "="^60)
        println("TEST 8: Biological Range & Bounds Validation")
        println("="^60)
        
        test_results["biological_ranges"] = test_biological_ranges(network)
        
        # Summary of biological realism assessment
        println("\\n" * "="^70)
        println("🧬 BIOLOGICAL REALISM ASSESSMENT SUMMARY")
        println("="^70)
        
        passed_tests = sum(values(test_results))
        total_tests = length(test_results)
        success_rate = (passed_tests / total_tests) * 100
        
        println("\\nTest Results:")
        for (test_name, passed) in test_results
            status = passed ? "✅ PASS" : "❌ FAIL"
            println("  $(rpad(test_name, 25)): $status")
        end
        
        println("\\n📊 Overall Biological Realism Score: $(round(success_rate, digits=1))% ($passed_tests/$total_tests)")
        
        # Biological realism categories
        if success_rate >= 90.0
            println("🎉 EXCELLENT biological realism - System behaves as expected")
            println("   ✅ Ready for biological pathway analysis")
        elseif success_rate >= 75.0
            println("✅ GOOD biological realism - Minor issues identified")
            println("   ⚠️  Some biological behaviors may need refinement")
        elseif success_rate >= 50.0
            println("⚠️  MODERATE biological realism - Several issues found")
            println("   🔧 Significant biological tuning recommended")
        else
            println("❌ POOR biological realism - Major issues detected")
            println("   🚨 System needs substantial biological validation")
        end
        
        return success_rate >= 75.0
        
    catch e
        println("❌ Biological realism test failed: $e")
        println("Stack trace:")
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
        return false
    end
end

"""
Test 1: Dose-response curves should be monotonic (higher input → higher output for activators)
"""
function test_dose_response_monotonicity(network::ReactionNetwork)::Bool
    println("🧪 Testing dose-response monotonicity...")
    
    # Test ROOT1 → A1 → TERM1 pathway (should be monotonically increasing)
    dose_levels = [20.0, 40.0, 60.0, 80.0]
    term1_responses = Float64[]
    
    ss_params = SteadyStateParams(1.0, 0.05, 150, 1e-5, "fixed_point")
    
    for dose in dose_levels
        observations = Dict("ROOT1" => (dose, 0.9))
        result = solve_steady_state(network, observations, ss_params)
        
        if result.converged
            term1_activity = get(result.node_activities, "TERM1", 0.2) * 100
            push!(term1_responses, term1_activity)
        else
            println("  ❌ Failed to converge at dose $dose")
            return false
        end
    end
    
    println("  📈 ROOT1 doses: $(dose_levels)")
    println("  📊 TERM1 responses: $([round(r, digits=1) for r in term1_responses])")
    
    # Check monotonicity
    is_monotonic = true
    for i in 2:length(term1_responses)
        if term1_responses[i] <= term1_responses[i-1]
            is_monotonic = false
            break
        end
    end
    
    # Check reasonable response range
    response_range = maximum(term1_responses) - minimum(term1_responses)
    has_good_range = response_range >= 10.0  # At least 10% range
    
    if is_monotonic && has_good_range
        println("  ✅ Dose-response is monotonic with good range ($(round(response_range, digits=1))%)")
        return true
    elseif is_monotonic
        println("  ⚠️  Dose-response is monotonic but range is limited ($(round(response_range, digits=1))%)")
        return true  # Still acceptable
    else
        println("  ❌ Dose-response is NOT monotonic - violates biological expectation")
        return false
    end
end

"""
Test 2: Signal should attenuate with distance from source
"""
function test_signal_attenuation(network::ReactionNetwork)::Bool
    println("🛡️  Testing signal attenuation with distance...")
    
    # Test ROOT7 → I1 → J1 → K1 → L1 → TERM6 (5-step cascade)
    observations = Dict("ROOT7" => (80.0, 0.9))
    
    ss_params = SteadyStateParams(1.0, 0.05, 150, 1e-5, "fixed_point")
    result = solve_steady_state(network, observations, ss_params)
    
    if !result.converged
        println("  ❌ Failed to converge")
        return false
    end
    
    # Measure activities along the cascade
    cascade_nodes = ["ROOT7", "I1", "J1", "K1", "L1", "TERM6"]
    activities = Float64[]
    
    for node in cascade_nodes
        activity = get(result.node_activities, node, 0.2) * 100
        push!(activities, activity)
    end
    
    println("  🔗 Cascade path: $(join(cascade_nodes, " → "))")
    println("  📊 Activities: $([round(a, digits=1) for a in activities])")
    
    # Check that signal generally decreases (allowing some fluctuation)
    source_activity = activities[1]  # ROOT7
    terminal_activity = activities[end]  # TERM6
    
    # Signal should be attenuated but still significant
    attenuation_reasonable = (terminal_activity < source_activity * 1.1) && (terminal_activity > source_activity * 0.3)
    
    # No major amplification in intermediate steps
    no_major_amplification = true
    for i in 2:length(activities)
        if activities[i] > activities[1] * 1.5  # Allow some reasonable amplification
            no_major_amplification = false
            break
        end
    end
    
    signal_transmitted = terminal_activity > 30.0  # Signal should reach the end
    
    if attenuation_reasonable && no_major_amplification && signal_transmitted
        println("  ✅ Signal attenuation is biologically reasonable")
        println("    Source: $(round(source_activity, digits=1))% → Terminal: $(round(terminal_activity, digits=1))%")
        return true
    else
        println("  ❌ Signal attenuation is not biologically reasonable")
        println("    Attenuation OK: $attenuation_reasonable, No amplification: $no_major_amplification, Signal transmitted: $signal_transmitted")
        return false
    end
end

"""
Test 3: Competitive inhibition - inhibitor should reduce activator effect
"""
function test_competitive_inhibition(network::ReactionNetwork)::Bool
    println("🥊 Testing competitive inhibition...")
    
    # Test ROOT5(+) and ROOT6(-) competing at G1
    ss_params = SteadyStateParams(1.0, 0.05, 150, 1e-5, "fixed_point")
    
    # Scenario 1: Activator only
    obs1 = Dict("ROOT5" => (80.0, 0.9))
    result1 = solve_steady_state(network, obs1, ss_params)
    
    # Scenario 2: Inhibitor only  
    obs2 = Dict("ROOT6" => (80.0, 0.9))
    result2 = solve_steady_state(network, obs2, ss_params)
    
    # Scenario 3: Both competing
    obs3 = Dict("ROOT5" => (80.0, 0.9), "ROOT6" => (80.0, 0.9))
    result3 = solve_steady_state(network, obs3, ss_params)
    
    if !(result1.converged && result2.converged && result3.converged)
        println("  ❌ Some simulations failed to converge")
        return false
    end
    
    # Get terminal output (TERM5)
    term5_activator = get(result1.node_activities, "TERM5", 0.2) * 100
    term5_inhibitor = get(result2.node_activities, "TERM5", 0.2) * 100
    term5_competing = get(result3.node_activities, "TERM5", 0.2) * 100
    
    println("  📊 TERM5 activities:")
    println("    Activator only (ROOT5): $(round(term5_activator, digits=1))%")
    println("    Inhibitor only (ROOT6): $(round(term5_inhibitor, digits=1))%")
    println("    Both competing:        $(round(term5_competing, digits=1))%")
    
    # Biological expectations:
    # 1. Activator only should be higher than baseline (~20%)
    activator_effective = term5_activator > 30.0
    
    # 2. Inhibitor only should be lower than activator only
    inhibitor_effective = term5_inhibitor < term5_activator
    
    # 3. Competition should result in intermediate value
    competition_intermediate = (term5_competing > term5_inhibitor) && (term5_competing < term5_activator)
    
    # 4. Competition should reduce activator effect significantly
    significant_competition = term5_competing < term5_activator * 0.8
    
    all_tests_pass = activator_effective && inhibitor_effective && competition_intermediate && significant_competition
    
    if all_tests_pass
        println("  ✅ Competitive inhibition behaves biologically correctly")
        competition_strength = ((term5_activator - term5_competing) / term5_activator) * 100
        println("    Competition reduces activator effect by $(round(competition_strength, digits=1))%")
        return true
    else
        println("  ❌ Competitive inhibition issues detected:")
        println("    Activator effective: $activator_effective")
        println("    Inhibitor effective: $inhibitor_effective") 
        println("    Competition intermediate: $competition_intermediate")
        println("    Significant competition: $significant_competition")
        return false
    end
end

"""
Test 4: Feedback loops should behave as expected biologically
"""
function test_feedback_loop_behavior(network::ReactionNetwork)::Bool
    println("🔄 Testing feedback loop biological behavior...")
    
    ss_params_std = SteadyStateParams(1.0, 0.05, 150, 1e-5, "fixed_point")
    ss_params_enh = SteadyStateParams(1.0, 0.1, 200, 1e-5, "feedback_aware")
    
    # Test negative feedback: ROOT10 → Q1 → R1 → S1 → Q1(-)
    println("\\n  🔻 Testing NEGATIVE feedback (ROOT10 → Q1 → R1 → S1 → Q1(-)):")
    
    upregulation = Dict("ROOT10" => (80.0, 0.9))
    downregulation = Dict("ROOT10" => (40.0, 0.9))
    
    # Standard solver
    result_up_std = solve_steady_state(network, upregulation, ss_params_std)
    result_down_std = solve_steady_state(network, downregulation, ss_params_std)
    
    # Enhanced solver
    result_up_enh = solve_steady_state_enhanced(network, upregulation, ss_params_enh)
    result_down_enh = solve_steady_state_enhanced(network, downregulation, ss_params_enh)
    
    if !(result_up_std.converged && result_down_std.converged && result_up_enh.converged && result_down_enh.converged)
        println("    ❌ Some negative feedback simulations failed to converge")
        return false
    end
    
    # Analyze negative feedback dampening
    baseline = 20.0  # Expected baseline
    
    term8_up_std = get(result_up_std.node_activities, "TERM8", 0.2) * 100
    term8_up_enh = get(result_up_enh.node_activities, "TERM8", 0.2) * 100
    
    input_deviation = abs(80.0 - baseline)
    output_deviation_std = abs(term8_up_std - baseline)
    output_deviation_enh = abs(term8_up_enh - baseline)
    
    dampening_std = output_deviation_std / input_deviation
    dampening_enh = abs(output_deviation_enh) / input_deviation
    
    println("    📊 Negative feedback analysis (ROOT10: 80%):")
    println("      Standard solver - TERM8: $(round(term8_up_std, digits=1))%, dampening: $(round(dampening_std, digits=2))x")
    println("      Enhanced solver - TERM8: $(round(term8_up_enh, digits=1))%, dampening: $(round(dampening_enh, digits=2))x")
    
    # Test positive feedback: ROOT4 → D1 → E1 → F1 → D1(+)
    println("\\n  🔺 Testing POSITIVE feedback (ROOT4 → D1 → E1 → F1 → D1(+)):")
    
    pos_fb_obs = Dict("ROOT4" => (60.0, 0.9))
    result_pos_std = solve_steady_state(network, pos_fb_obs, ss_params_std)
    result_pos_enh = solve_steady_state_enhanced(network, pos_fb_obs, ss_params_enh)
    
    if !(result_pos_std.converged && result_pos_enh.converged)
        println("    ❌ Positive feedback simulations failed to converge")
        return false
    end
    
    term4_pos_std = get(result_pos_std.node_activities, "TERM4", 0.2) * 100
    term4_pos_enh = get(result_pos_enh.node_activities, "TERM4", 0.2) * 100
    
    println("    📊 Positive feedback analysis (ROOT4: 60%):")
    println("      Standard solver - TERM4: $(round(term4_pos_std, digits=1))%")
    println("      Enhanced solver - TERM4: $(round(term4_pos_enh, digits=1))%")
    
    # Biological validation criteria
    negative_fb_working = dampening_enh < 0.8  # Enhanced solver should dampen
    positive_fb_amplifying = term4_pos_enh > 40.0  # Should show some amplification
    enhanced_better = dampening_enh < dampening_std  # Enhanced should be better
    
    if negative_fb_working && positive_fb_amplifying && enhanced_better
        println("  ✅ Feedback loops behave biologically correctly")
        improvement = ((dampening_std - dampening_enh) / dampening_std) * 100
        println("    Enhanced solver improves negative feedback by $(round(improvement, digits=1))%")
        return true
    else
        println("  ❌ Feedback loop issues detected:")
        println("    Negative FB dampening: $negative_fb_working ($(round(dampening_enh, digits=2))x)")
        println("    Positive FB amplifying: $positive_fb_amplifying ($(round(term4_pos_enh, digits=1))%)")
        println("    Enhanced solver better: $enhanced_better")
        return false
    end
end

"""
Test 5: Time-dynamic patterns should show realistic kinetics
"""
function test_dynamic_response_patterns(network::ReactionNetwork)::Bool
    println("⏰ Testing time-dynamic response patterns...")
    
    # Test step response
    observations = Dict("ROOT1" => (70.0, 0.9))
    td_params = TimeDynamicParams(0.0, 5.0, 0.1)
    
    result = rollout_time_dynamic(network, observations, td_params)
    
    if !result.converged
        println("  ❌ Time-dynamic simulation failed to converge")
        return false
    end
    
    # Analyze TERM1 response
    if !haskey(result.node_activities, "TERM1")
        println("  ❌ TERM1 not found in time-dynamic results")
        return false
    end
    
    term1_timeseries = result.node_activities["TERM1"] .* 100
    time_points = result.time_points
    
    # Extract key metrics
    initial_value = term1_timeseries[1]
    final_value = term1_timeseries[end]
    max_value = maximum(term1_timeseries)
    
    # Find rise time (time to reach 63% of final value)
    target_value = initial_value + 0.63 * (final_value - initial_value)
    rise_time_idx = findfirst(x -> x >= target_value, term1_timeseries)
    rise_time = rise_time_idx !== nothing ? time_points[rise_time_idx] : NaN
    
    println("  📊 TERM1 step response analysis:")
    println("    Initial: $(round(initial_value, digits=1))%")
    println("    Final: $(round(final_value, digits=1))%")
    println("    Max: $(round(max_value, digits=1))%")
    println("    Rise time (63%): $(isnan(rise_time) ? "N/A" : round(rise_time, digits=2))s")
    
    # Biological validation criteria
    shows_response = abs(final_value - initial_value) > 5.0  # At least 5% change
    no_overshoot = max_value <= final_value * 1.1  # Minimal overshoot
    reasonable_kinetics = !isnan(rise_time) && rise_time > 0.1 && rise_time < 3.0
    stable_final = abs(term1_timeseries[end] - term1_timeseries[end-5]) < 1.0  # Stable in last 0.5s
    
    if shows_response && no_overshoot && reasonable_kinetics && stable_final
        println("  ✅ Time-dynamic patterns are biologically realistic")
        return true
    else
        println("  ❌ Time-dynamic pattern issues:")
        println("    Shows response: $shows_response")
        println("    No overshoot: $no_overshoot") 
        println("    Reasonable kinetics: $reasonable_kinetics")
        println("    Stable final: $stable_final")
        return false
    end
end

"""
Test 6: Independent pathways should not interfere significantly
"""
function test_pathway_crosstalk(network::ReactionNetwork)::Bool
    println("🔀 Testing pathway crosstalk and independence...")
    
    ss_params = SteadyStateParams(1.0, 0.05, 150, 1e-5, "fixed_point")
    
    # Test independent pathways: ROOT7→TERM6 and ROOT8→TERM7
    
    # Pathway 1 alone
    obs1 = Dict("ROOT7" => (70.0, 0.9))
    result1 = solve_steady_state(network, obs1, ss_params)
    
    # Pathway 2 alone  
    obs2 = Dict("ROOT8" => (60.0, 0.9))
    result2 = solve_steady_state(network, obs2, ss_params)
    
    # Both pathways together
    obs_both = Dict("ROOT7" => (70.0, 0.9), "ROOT8" => (60.0, 0.9))
    result_both = solve_steady_state(network, obs_both, ss_params)
    
    if !(result1.converged && result2.converged && result_both.converged)
        println("  ❌ Some crosstalk simulations failed to converge")
        return false
    end
    
    # Compare outputs
    term6_alone = get(result1.node_activities, "TERM6", 0.2) * 100
    term6_together = get(result_both.node_activities, "TERM6", 0.2) * 100
    
    term7_alone = get(result2.node_activities, "TERM7", 0.2) * 100
    term7_together = get(result_both.node_activities, "TERM7", 0.2) * 100
    
    println("  📊 Pathway independence analysis:")
    println("    TERM6: alone $(round(term6_alone, digits=1))% vs together $(round(term6_together, digits=1))%")
    println("    TERM7: alone $(round(term7_alone, digits=1))% vs together $(round(term7_together, digits=1))%")
    
    # Calculate interference
    term6_interference = abs(term6_together - term6_alone) / max(term6_alone, 1.0)
    term7_interference = abs(term7_together - term7_alone) / max(term7_alone, 1.0)
    
    println("    Interference: TERM6 $(round(term6_interference * 100, digits=1))%, TERM7 $(round(term7_interference * 100, digits=1))%")
    
    # Biological expectation: independent pathways should have minimal interference
    low_crosstalk = (term6_interference < 0.2) && (term7_interference < 0.2)  # <20% interference
    
    if low_crosstalk
        println("  ✅ Pathway independence is maintained (low crosstalk)")
        return true
    else
        println("  ❌ Excessive pathway crosstalk detected")
        return false
    end
end

"""
Test 7: Hill kinetics should show proper saturation behavior
"""
function test_hill_kinetics_behavior(network::ReactionNetwork)::Bool
    println("📈 Testing Hill kinetics and saturation behavior...")
    
    # Test saturation by increasing ROOT1 to very high levels
    dose_levels = [20.0, 40.0, 60.0, 80.0, 90.0, 95.0]
    responses = Float64[]
    
    ss_params = SteadyStateParams(1.0, 0.05, 150, 1e-5, "fixed_point")
    
    for dose in dose_levels
        observations = Dict("ROOT1" => (dose, 0.9))
        result = solve_steady_state(network, observations, ss_params)
        
        if result.converged
            term1_activity = get(result.node_activities, "TERM1", 0.2) * 100
            push!(responses, term1_activity)
        else
            println("  ❌ Failed to converge at dose $dose")
            return false
        end
    end
    
    println("  📊 Saturation analysis:")
    println("    Doses: $dose_levels")
    println("    Responses: $([round(r, digits=1) for r in responses])")
    
    # Calculate response increments
    increments = Float64[]
    for i in 2:length(responses)
        increment = responses[i] - responses[i-1]
        push!(increments, increment)
    end
    
    println("    Increments: $([round(inc, digits=1) for inc in increments])")
    
    # Biological validation: should show diminishing returns (saturation)
    shows_saturation = length(increments) >= 3 && increments[end] < increments[1] * 0.7
    reaches_reasonable_max = maximum(responses) > 40.0  # Should reach reasonable activity
    no_negative_responses = all(r -> r >= 0.0, responses)
    bounded_responses = all(r -> r <= 100.0, responses)
    
    if shows_saturation && reaches_reasonable_max && no_negative_responses && bounded_responses
        println("  ✅ Hill kinetics show proper saturation behavior")
        saturation_ratio = increments[end] / increments[1]
        println("    Saturation ratio: $(round(saturation_ratio, digits=2)) (lower = more saturation)")
        return true
    else
        println("  ❌ Hill kinetics issues detected:")
        println("    Shows saturation: $shows_saturation")
        println("    Reaches reasonable max: $reaches_reasonable_max")
        println("    No negative responses: $no_negative_responses")
        println("    Bounded responses: $bounded_responses")
        return false
    end
end

"""
Test 8: Biological ranges and bounds validation
"""
function test_biological_ranges(network::ReactionNetwork)::Bool
    println("📏 Testing biological ranges and bounds...")
    
    ss_params = SteadyStateParams(1.0, 0.05, 150, 1e-5, "fixed_point")
    
    # Test extreme scenarios to check bounds
    extreme_scenarios = [
        ("No input", Dict{String, Tuple{Float64, Float64}}()),
        ("Maximum input", Dict("ROOT1" => (100.0, 0.9), "ROOT2" => (100.0, 0.9))),
        ("Minimum input", Dict("ROOT1" => (0.0, 0.9), "ROOT2" => (0.0, 0.9))),
        ("Mixed extreme", Dict("ROOT1" => (100.0, 0.9), "ROOT6" => (100.0, 0.9)))  # Strong competition
    ]
    
    all_valid = true
    
    for (scenario_name, observations) in extreme_scenarios
        result = solve_steady_state(network, observations, ss_params)
        
        if !result.converged
            println("  ❌ Failed to converge for scenario: $scenario_name")
            all_valid = false
            continue
        end
        
        # Check all activities are in valid range [0, 1]
        activities = collect(values(result.node_activities))
        
        min_activity = minimum(activities)
        max_activity = maximum(activities)
        
        valid_range = (min_activity >= 0.0) && (max_activity <= 1.0)
        reasonable_diversity = (max_activity - min_activity) > 0.05  # Some diversity in responses
        
        println("  📊 $scenario_name: range [$(round(min_activity, digits=3)), $(round(max_activity, digits=3))]")
        
        if !valid_range
            println("    ❌ Activities outside valid range [0,1]")
            all_valid = false
        end
        
        # Check for NaN or Inf values
        has_nan_inf = any(x -> isnan(x) || isinf(x), activities)
        if has_nan_inf
            println("    ❌ Contains NaN or Inf values")
            all_valid = false
        end
    end
    
    # Test baseline behavior (no perturbations)
    baseline_result = solve_steady_state(network, Dict{String, Tuple{Float64, Float64}}(), ss_params)
    if baseline_result.converged
        baseline_activities = collect(values(baseline_result.node_activities)) .* 100
        baseline_mean = mean(baseline_activities)
        baseline_std = std(baseline_activities)
        
        println("  📊 Baseline behavior (no perturbations):")
        println("    Mean activity: $(round(baseline_mean, digits=1))%")
        println("    Std deviation: $(round(baseline_std, digits=1))%")
        
        # Baseline should be reasonable
        reasonable_baseline = (baseline_mean > 15.0) && (baseline_mean < 35.0)  # Around 20% ± 15%
        not_too_variable = baseline_std < 20.0  # Not extremely variable
        
        if !reasonable_baseline || !not_too_variable
            println("    ⚠️  Baseline behavior may be unrealistic")
            all_valid = false
        end
    else
        println("  ❌ Baseline simulation failed to converge")
        all_valid = false
    end
    
    if all_valid
        println("  ✅ All biological ranges and bounds are valid")
        return true
    else
        println("  ❌ Some biological range violations detected")
        return false
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_biological_realism()
end