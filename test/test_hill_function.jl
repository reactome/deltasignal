#!/usr/bin/env julia

# Test hill function behavior
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
include("../src/core/hill_functions.jl")
include("../src/core/aggregators.jl")

function test_hill_behavior()
    println("🔬 Testing Hill Function Behavior")
    println("=" ^ 50)
    
    # Current parameters from reaction_model.jl
    h = 1.5
    K = 0.3
    
    println("Current Hill parameters: h=$h, K=$K")
    println()
    
    # Test single input propagation
    println("Single Input Propagation (PI3K → AKT):")
    test_inputs = [0.01, 0.3, 0.5, 0.7, 1.0]  # 1%, 30%, 50%, 70%, 100%
    
    for input in test_inputs
        output = hill_activation(input, h, K)
        println("  Input: $(input*100)% → Output: $(round(output*100, digits=1))%")
    end
    
    println()
    println("Expected behavior:")
    println("  - Single input: Output should approximately equal input for biological relevance")
    println("  - Input 30% should give output ~30%, not $(round(hill_activation(0.3, h, K)*100, digits=1))%")
    println()
    
    # Test multiple input aggregation  
    println("Multiple Input Aggregation (A=50%, B=50%):")
    inputs = [0.5, 0.5]
    weights = [0.5, 0.5]
    
    # Current geometric mean
    geo_result = geometric_mean_aggregator(inputs, weights)
    hill_result = hill_activation(geo_result, h, K)
    
    println("  Geometric mean: $(round(geo_result*100, digits=1))%")
    println("  After hill function: $(round(hill_result*100, digits=1))%")
    println("  Expected: Should be close to 50% (both inputs equal)")
    
    # Test with better parameters
    println()
    println("🔧 Testing Better Parameters:")
    h_better = 2.0
    K_better = 0.1  # 10% threshold - balance between propagation and nonlinearity
    
    println("Better Hill parameters: h=$h_better, K=$K_better")
    println()
    
    println("Single Input Propagation (improved):")
    for input in test_inputs
        output = hill_activation(input, h_better, K_better)
        println("  Input: $(input*100)% → Output: $(round(output*100, digits=1))%")
    end
    
    println()
    println("Multiple Input Aggregation (A=50%, B=50%, improved):")
    geo_result = geometric_mean_aggregator(inputs, weights)
    hill_result = hill_activation(geo_result, h_better, K_better)
    
    println("  Geometric mean: $(round(geo_result*100, digits=1))%")
    println("  After hill function: $(round(hill_result*100, digits=1))%")
    
    println()
    println("💡 Recommendations:")
    println("  1. Use K ≈ 0.01 (1% threshold) instead of 0.3 (30%)")
    println("  2. For single inputs, consider direct propagation")
    println("  3. For multiple inputs, geometric mean is reasonable but needs low K")
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_hill_behavior()
end