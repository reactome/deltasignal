#!/usr/bin/env julia

# Test the new propagation behavior
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
include("../src/io/tsv_parser.jl")  # This includes the network types
include("../src/core/hill_functions.jl")
include("../src/core/aggregators.jl")
include("../src/core/sensitivity.jl")
include("../src/core/reaction_model.jl")

function test_single_vs_multi_input()
    println("🧪 Testing Single vs Multi-Input Propagation")
    println("=" ^ 50)
    
    # Test single input propagation (PI3K → AKT)
    println("Single Input Case (PI3K → AKT):")
    
    # Create a reaction with single activator
    single_reaction = Reaction(
        "AKT",
        ["PI3K"],  # Single activator
        String[],  # No inhibitors
        String[],  # No substrates
        String[],  # No products
        create_default_reaction_params(1, 0, 0),
        false
    )
    
    test_levels = [0.01, 0.1, 0.3, 0.5, 0.7, 1.0]
    
    for pi3k_level in test_levels
        activities = Dict("PI3K" => pi3k_level)
        akt_output = compute_reaction_output(single_reaction, activities)
        
        println("  PI3K: $(round(pi3k_level*100, digits=1))% → AKT: $(round(akt_output*100, digits=1))%")
    end
    
    println()
    println("Multi-Input Case (A + B → C):")
    
    # Create a reaction with multiple activators
    multi_reaction = Reaction(
        "C",
        ["A", "B"],  # Two activators
        String[],    # No inhibitors
        String[],    # No substrates
        String[],    # No products
        create_default_reaction_params(2, 0, 0),
        false
    )
    
    test_cases = [
        (0.3, 0.3),  # Both at 30%
        (0.5, 0.5),  # Both at 50%
        (0.8, 0.2),  # Unbalanced
        (1.0, 1.0),  # Both at 100%
    ]
    
    for (a_level, b_level) in test_cases
        activities = Dict("A" => a_level, "B" => b_level)
        c_output = compute_reaction_output(multi_reaction, activities)
        
        expected_geo_mean = sqrt(a_level * b_level)
        
        println("  A: $(round(a_level*100, digits=1))%, B: $(round(b_level*100, digits=1))% → C: $(round(c_output*100, digits=1))% (geo mean would be $(round(expected_geo_mean*100, digits=1))%)")
    end
    
    println()
    println("Expected Behavior:")
    println("  ✓ Single input: Output ≈ Input (with small nonlinearity)")
    println("  ✓ Multi-input: Geometric mean + Hill function for integration")
    
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_single_vs_multi_input()
end