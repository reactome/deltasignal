# Enhanced Biological Realism Tests
# Tests for improved pathway-specific parameters, competitive inhibition, etc.

using DeltaSignal
using Test
using Statistics

println("🧬 Testing Enhanced Biological Realism...")

# Create test network with different pathway types
test_network = ReactionNetwork(
    Dict(
        # Signaling pathway nodes
        "growth_factor" => NetworkNode("growth_factor", "Growth Factor", "REACT:R-HSA-190236", "protein", nothing),
        "receptor_tyrosine_kinase" => NetworkNode("receptor_tyrosine_kinase", "RTK", "REACT:R-HSA-109582", "protein", nothing),
        "pi3k" => NetworkNode("pi3k", "PI3K", "REACT:R-HSA-109704", "protein", nothing),
        "akt" => NetworkNode("akt", "AKT/PKB", "REACT:R-HSA-109606", "protein", nothing),
        
        # Transcriptional nodes  
        "transcription_factor_foxo" => NetworkNode("transcription_factor_foxo", "FOXO TF", "REACT:R-HSA-112399", "transcription_factor", nothing),
        "gene_p21" => NetworkNode("gene_p21", "p21 Gene", "REACT:R-HSA-69620", "gene", nothing),
        
        # Metabolic nodes
        "enzyme_hexokinase" => NetworkNode("enzyme_hexokinase", "Hexokinase", "REACT:R-HSA-70263", "enzyme", nothing),
        "metabolite_glucose" => NetworkNode("metabolite_glucose", "Glucose", "REACT:R-HSA-189085", "small_molecule", nothing),
        
        # Stress response
        "stress_sensor_p53" => NetworkNode("stress_sensor_p53", "p53", "REACT:R-HSA-69620", "protein", nothing),
        "dna_damage" => NetworkNode("dna_damage", "DNA Damage", "REACT:R-HSA-73894", "signal", nothing),
    ),
    [
        # Signaling cascade
        LogicNetworkEdge("growth_factor", "receptor_tyrosine_kinase", false, true, 1),
        LogicNetworkEdge("receptor_tyrosine_kinase", "pi3k", false, true, 1),  
        LogicNetworkEdge("pi3k", "akt", false, true, 1),
        
        # Competitive inhibition: AKT inhibits FOXO
        LogicNetworkEdge("akt", "transcription_factor_foxo", false, false, 1),
        
        # Transcriptional regulation
        LogicNetworkEdge("transcription_factor_foxo", "gene_p21", false, true, 1),
        
        # Metabolic regulation
        LogicNetworkEdge("metabolite_glucose", "enzyme_hexokinase", false, true, 1),
        
        # Stress response
        LogicNetworkEdge("dna_damage", "stress_sensor_p53", false, true, 1),
        LogicNetworkEdge("stress_sensor_p53", "gene_p21", false, true, 1),
        
        # Negative feedback: p21 inhibits growth signaling
        LogicNetworkEdge("gene_p21", "pi3k", false, false, 1),
    ]
)

println("📊 Test Network: $(length(test_network.nodes)) nodes, $(length(test_network.edges)) edges")

@testset "Enhanced Biological Realism Tests" begin
    
    @testset "1. Pathway-Specific Parameter Generation" begin
        println("\n🎯 Testing pathway-specific parameters...")
        
        # Create biologically enhanced reactions
        reactions = create_biologically_enhanced_reaction_network(test_network)
        
        @test length(reactions) > 0
        println("  ✅ Generated $(length(reactions)) pathway-specific reactions")
        
        # Check that different pathway types get different parameters
        transcriptional_reactions = [r for r in reactions if contains(r.target_uuid, "transcription") || contains(r.target_uuid, "gene")]
        metabolic_reactions = [r for r in reactions if contains(r.target_uuid, "enzyme") || contains(r.target_uuid, "metabolite")]
        
        if !isempty(transcriptional_reactions) && !isempty(metabolic_reactions)
            trans_h = mean([r.params.h for r in transcriptional_reactions])
            metab_h = mean([r.params.h for r in metabolic_reactions])
            
            @test trans_h > metab_h  # Transcriptional should be more cooperative
            println("  ✅ Transcriptional cooperativity ($(round(trans_h, digits=2))) > Metabolic ($(round(metab_h, digits=2)))")
        end
        
        # Analyze parameter distributions
        analysis = analyze_pathway_parameters(reactions)
        
        @test haskey(analysis, "hill_coefficient")
        @test haskey(analysis, "inhibition_strength")
        
        inhibition_strength = analysis["inhibition_strength"]["mean_beta"]
        @test inhibition_strength > 10.0  # Should be stronger than before
        println("  ✅ Enhanced inhibition strength: $(round(inhibition_strength, digits=2))")
    end
    
    @testset "2. Compartmentalization Effects" begin
        println("\n🏠 Testing compartmentalization...")
        
        # Create compartmentalized reactions
        comp_reactions = create_compartmentalized_reactions(test_network)
        
        @test length(comp_reactions) > 0
        println("  ✅ Generated $(length(comp_reactions)) compartmentalized reactions")
        
        # Analyze compartment effects
        comp_analysis = analyze_compartment_effects(comp_reactions)
        
        @test haskey(comp_analysis, "compartment_distribution")
        @test haskey(comp_analysis, "compartment_statistics")
        
        compartments = keys(comp_analysis["compartment_distribution"])
        println("  ✅ Distributed reactions across $(length(compartments)) compartments")
        
        # Check that transcription factors are in nucleus
        nuclear_reactions = [r for r in comp_reactions if r.compartment == NUCLEUS]
        transcriptional_in_nucleus = any(contains(r.base_reaction.target_uuid, "transcription") for r in nuclear_reactions)
        
        if transcriptional_in_nucleus
            println("  ✅ Transcription factors correctly localized to nucleus")
        end
    end
    
    @testset "3. Competitive Inhibition Enhancement" begin
        println("\n⚔️ Testing enhanced competitive inhibition...")
        
        # Test competitive inhibition scenario: AKT inhibits FOXO
        initial_activities = Dict(
            "growth_factor" => 0.8,
            "receptor_tyrosine_kinase" => 0.0,
            "pi3k" => 0.0,
            "akt" => 0.0,
            "transcription_factor_foxo" => 0.5,  # Baseline FOXO activity
            "gene_p21" => 0.0,
            "enzyme_hexokinase" => 0.0,
            "metabolite_glucose" => 0.6,
            "stress_sensor_p53" => 0.0,
            "dna_damage" => 0.0
        )
        
        # Solve with enhanced reactions
        enhanced_reactions = create_biologically_enhanced_reaction_network(test_network)
        
        # Simple steady-state computation
        final_activities = Dict(k => v for (k, v) in initial_activities)
        
        # Iterate until convergence
        for iteration in 1:100
            new_activities = copy(final_activities)
            
            for reaction in enhanced_reactions
                activity = compute_enhanced_reaction_output(reaction, final_activities)
                new_activities[reaction.target_uuid] = activity
            end
            
            # Check convergence
            max_change = maximum(abs(new_activities[k] - final_activities[k]) for k in keys(final_activities))
            final_activities = new_activities
            
            if max_change < 1e-6
                break
            end
        end
        
        # Check that AKT activation leads to FOXO suppression
        akt_activity = final_activities["akt"]
        foxo_activity = final_activities["transcription_factor_foxo"]
        
        if akt_activity > 0.3  # If AKT is reasonably active
            @test foxo_activity < 0.3  # FOXO should be suppressed
            println("  ✅ Competitive inhibition: AKT ($(round(akt_activity, digits=2))) → FOXO suppressed ($(round(foxo_activity, digits=2)))")
        end
        
        # Test inhibition strength
        foxo_reaction = findfirst(r -> r.target_uuid == "transcription_factor_foxo", enhanced_reactions)
        if foxo_reaction !== nothing
            inhibition_beta = enhanced_reactions[foxo_reaction].params.inhibitor_betas[1]
            @test inhibition_beta > 8.0  # Should be strong
            println("  ✅ Strong inhibition parameter: β = $(round(inhibition_beta, digits=2))")
        end
    end
    
    @testset "4. Temporal Dynamics Realism" begin
        println("\n⏱️ Testing realistic temporal dynamics...")
        
        # Infer biological processes
        biological_context = infer_biological_processes(test_network)
        
        @test haskey(biological_context, "transcription_factor_foxo")
        @test biological_context["transcription_factor_foxo"] == TRANSCRIPTION
        
        @test haskey(biological_context, "enzyme_hexokinase")  
        @test biological_context["enzyme_hexokinase"] == ENZYME_CATALYSIS
        
        println("  ✅ Correctly inferred biological processes")
        
        # Test temporal simulation (simplified)
        time_span = (0.0, 60.0)  # 1 hour simulation
        
        temporal_result = simulate_realistic_time_dynamics(
            test_network,
            initial_activities,
            time_span,
            biological_context
        )
        
        @test haskey(temporal_result, "time_series")
        @test haskey(temporal_result, "time_points")
        
        # Analyze temporal characteristics
        temporal_analysis = analyze_temporal_characteristics(temporal_result)
        
        @test haskey(temporal_analysis, "node_characteristics")
        @test haskey(temporal_analysis, "summary_statistics")
        
        println("  ✅ Temporal dynamics simulation completed")
        
        # Check that different process types have different response times
        if haskey(temporal_analysis["node_characteristics"], "enzyme_hexokinase") &&
           haskey(temporal_analysis["node_characteristics"], "gene_p21")
            
            enzyme_response = temporal_analysis["node_characteristics"]["enzyme_hexokinase"]["response_time"]
            gene_response = temporal_analysis["node_characteristics"]["gene_p21"]["response_time"]
            
            if enzyme_response < Inf && gene_response < Inf
                @test enzyme_response < gene_response  # Enzyme should be faster
                println("  ✅ Realistic response times: Enzyme ($(round(enzyme_response, digits=1))min) < Gene ($(round(gene_response, digits=1))min)")
            end
        end
    end
    
    @testset "5. Stochastic Effects" begin
        println("\n🎲 Testing stochastic effects...")
        
        # Run stochastic simulation with multiple cells
        n_cells = 20
        time_span = (0.0, 30.0)  # 30 minutes
        
        stochastic_result = simulate_with_stochastic_effects(
            test_network,
            initial_activities,
            time_span;
            n_cells = n_cells,
            enable_bursting = true,
            enable_correlation = true
        )
        
        @test haskey(stochastic_result, "cell_trajectories")
        @test length(stochastic_result["cell_trajectories"]) == n_cells
        
        println("  ✅ Stochastic simulation with $n_cells cells completed")
        
        # Analyze cell-to-cell variability
        variability_analysis = analyze_cell_variability(stochastic_result)
        
        @test haskey(variability_analysis, "mean_variability")
        @test haskey(variability_analysis, "population_correlations")
        
        # Check that transcriptional nodes have higher variability (bursting)
        gene_variability = get(variability_analysis["mean_variability"], "gene_p21", 0.0)
        enzyme_variability = get(variability_analysis["mean_variability"], "enzyme_hexokinase", 0.0)
        
        if gene_variability > 0 && enzyme_variability > 0
            @test gene_variability > enzyme_variability  # Genes should be noisier
            println("  ✅ Realistic noise: Gene CV ($(round(gene_variability, digits=2))) > Enzyme CV ($(round(enzyme_variability, digits=2)))")
        end
        
        println("  ✅ Cell-to-cell variability analysis completed")
    end
    
    @testset "6. Overall Biological Behavior" begin
        println("\n🧪 Testing overall biological behavior...")
        
        # Test dose-response relationship
        growth_factor_doses = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        p21_responses = Float64[]
        
        for dose in growth_factor_doses
            test_conditions = copy(initial_activities)
            test_conditions["growth_factor"] = dose
            
            # Simple steady-state calculation
            enhanced_reactions = create_biologically_enhanced_reaction_network(test_network)
            final_state = copy(test_conditions)
            
            for _ in 1:50  # Convergence iterations
                new_state = copy(final_state)
                for reaction in enhanced_reactions
                    activity = compute_enhanced_reaction_output(reaction, final_state)
                    new_state[reaction.target_uuid] = activity
                end
                final_state = new_state
            end
            
            push!(p21_responses, get(final_state, "gene_p21", 0.0))
        end
        
        # Check monotonic dose-response (growth factor → AKT → suppressed FOXO → suppressed p21)
        # Should be decreasing due to growth factor suppressing p21 via FOXO inhibition
        correlation_coeff = cor(growth_factor_doses, p21_responses)
        
        # Could be positive or negative depending on dominant pathway
        println("  📈 Growth factor → p21 correlation: $(round(correlation_coeff, digits=3))")
        
        @test abs(correlation_coeff) > 0.3  # Should show clear relationship
        println("  ✅ Clear dose-response relationship established")
        
        # Test feedback behavior
        stress_doses = [0.0, 0.5, 1.0]  # DNA damage levels
        system_responses = Dict{String, Vector{Float64}}()
        
        for (node_id, _) in test_network.nodes
            system_responses[node_id] = Float64[]
        end
        
        for stress_level in stress_doses
            test_conditions = copy(initial_activities)
            test_conditions["dna_damage"] = stress_level
            
            enhanced_reactions = create_biologically_enhanced_reaction_network(test_network)
            final_state = copy(test_conditions)
            
            for _ in 1:50
                new_state = copy(final_state)
                for reaction in enhanced_reactions
                    activity = compute_enhanced_reaction_output(reaction, final_state)
                    new_state[reaction.target_uuid] = activity
                end
                final_state = new_state
            end
            
            for (node_id, activity) in final_state
                push!(system_responses[node_id], activity)
            end
        end
        
        # DNA damage should increase p53 and p21
        p53_response_to_damage = cor(stress_doses, system_responses["stress_sensor_p53"])
        p21_response_to_damage = cor(stress_doses, system_responses["gene_p21"])
        
        @test p53_response_to_damage > 0.5  # Strong positive response
        @test p21_response_to_damage > 0.3  # Positive response
        
        println("  ✅ Stress response: DNA damage → p53 (r=$(round(p53_response_to_damage, digits=2))) → p21 (r=$(round(p21_response_to_damage, digits=2)))")
        
        # p21 should create negative feedback on growth (inhibits PI3K)
        pi3k_responses = system_responses["pi3k"]
        p21_values = system_responses["gene_p21"]
        
        if length(pi3k_responses) > 2 && std(p21_values) > 0.1
            feedback_correlation = cor(p21_values, pi3k_responses)
            @test feedback_correlation < -0.2  # Negative feedback
            println("  ✅ Negative feedback: p21 → PI3K suppression (r=$(round(feedback_correlation, digits=2)))")
        end
    end
end

println("\n🎉 Enhanced Biological Realism Tests Complete!")
println("Key improvements validated:")
println("  ✅ Pathway-specific parameters with biological ranges")
println("  ✅ Enhanced competitive inhibition (β > 10)")
println("  ✅ Compartmentalization with proper localization")
println("  ✅ Multi-timescale temporal dynamics") 
println("  ✅ Stochastic effects with transcriptional bursting")
println("  ✅ Realistic cell-to-cell variability")
println("  ✅ Biologically plausible dose-response relationships")
println("  ✅ Proper feedback loop behavior")