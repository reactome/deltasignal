#!/usr/bin/env julia

# Basic test of DeltaSignal parsing functionality
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV
using DataFrames
using JSON3

# Add the src directory to the load path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
include("../src/io/tsv_parser.jl")

function test_basic_parsing()
    println("🧬 Testing DeltaSignal TSV parsing...")
    
    try
        # Test parsing logic network
        logic_path = joinpath(@__DIR__, "..", "examples", "sample_logic_network.tsv")
        uuid_path = joinpath(@__DIR__, "..", "examples", "sample_uuid_mapping.tsv")
        set_path = joinpath(@__DIR__, "..", "examples", "sample_set_mappings.tsv")
        
        println("📁 Loading files:")
        println("  Logic network: $logic_path")
        println("  UUID mapping: $uuid_path")
        println("  Set mappings: $set_path")
        
        edges = parse_logic_network(logic_path)
        println("✅ Parsed $(length(edges)) edges")
        
        nodes = parse_uuid_mapping(uuid_path)
        println("✅ Parsed $(length(nodes)) nodes")
        
        set_mappings = parse_set_mappings(set_path)
        println("✅ Parsed $(length(set_mappings)) set mappings")
        
        network = create_reaction_network(edges, nodes, set_mappings)
        println("✅ Created reaction network")
        
        # Print some details
        println("\n📊 Network details:")
        println("  Nodes: $(length(network.nodes))")
        println("  Edges: $(length(network.edges))")
        println("  Set mappings: $(length(network.set_mappings))")
        
        # Show first few nodes
        println("\n🔍 Sample nodes:")
        for (i, (uuid, node)) in enumerate(network.nodes)
            if i <= 3
                println("  $uuid -> $(node.display_name) ($(node.entity_type))")
            end
        end
        
        # Show first few edges
        println("\n🔍 Sample edges:")
        for (i, edge) in enumerate(network.edges)
            if i <= 3
                logic_type = edge.is_and ? "AND" : "OR"
                effect = edge.is_positive ? "+" : "-"
                println("  $(edge.parent_uuid) --($effect/$logic_type)--> $(edge.child_uuid)")
            end
        end
        
        println("\n✨ Basic parsing test completed successfully!")
        return true
        
    catch e
        println("❌ Test failed: $e")
        return false
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_basic_parsing()
end