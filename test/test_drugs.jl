# specs/032: drug-derived entities held at baseline (DS_DRUG_MODE=inert).
#
# Fixture = the RAF shape traced in specs/032: a kinase X and a drug D (a root)
# assemble into the drug-bound complex C, and C inhibits the reaction by which
# X makes its target T (T = X AND W). Left to propagate, raising X raises C,
# which cancels the rise at T: 2x reads 1x. With drugs inert, C stays at
# baseline and T follows X.
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
edge(s, t, pos, et; and = true) = DS.LogicNetworkEdge(s, t, and, pos, 1.0, et)

function fixture(drug_stids)
    nodes = Dict(node.(["X", "D", "C", "W", "T"]))
    edges = [edge("X", "C", true, "assembly"), edge("D", "C", true, "assembly"),
             edge("X", "T", true, "input"), edge("W", "T", true, "input"),
             edge("C", "T", false, "regulator"; and = false)]
    DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                       Dict{String, Set{String}}(), drug_stids)
end

function solve(net; mode, x = 2.0, extra = Dict{String, Tuple{Float64, Float64}}())
    with_env("DS_DRUG_MODE" => mode, "DS_SELF_INHIBITOR_WEIGHT" => "off",
             "DS_ASSEMBLY_LIMITING" => nothing) do
        obs = merge(Dict("X" => (x, 1.0), "W" => (1.0, 1.0)), extra)
        DS.solve_steady_state(net, obs, DS.SteadyStateParams(1.0, 0.1, 500, 1e-6, "penalty"))
    end
end
fold(r, u) = r.node_activities[u] / BL

@testset "inert drugs (specs/032)" begin

@testset "default is propagate, byte-identical to a network with no drug list" begin
    with_env("DS_DRUG_MODE" => nothing) do
        @test DS.drug_mode() == "propagate"
    end
    a = solve(fixture(Set(["R-C", "R-D"])); mode = nothing)
    b = solve(fixture(nothing); mode = nothing)
    @test a.node_activities == b.node_activities
    @test a.diagnostics["drug_rule"] == "propagate"
    @test a.diagnostics["drugs_held"] == 0
    # The mechanism the rule exists for: the drug-bound complex rises with X
    # and its inhibition cancels X's rise at T.
    @test fold(a, "C") ≈ 2.0 rtol = 1e-6
    @test fold(a, "T") ≈ 1.0 rtol = 1e-6
end

@testset "inert holds drug-derived nodes at baseline and T follows X" begin
    for x in (2.0, 0.5, 0.0, 50.0)
        r = solve(fixture(Set(["R-C", "R-D"])); mode = "inert", x = x)
        @test fold(r, "C") ≈ 1.0 rtol = 1e-9
        @test fold(r, "D") ≈ 1.0 rtol = 1e-9
        @test fold(r, "T") ≈ min(x, 100.0) rtol = 1e-6
        @test r.diagnostics["drug_rule"] == "inert"
        @test r.diagnostics["drugs_held"] == 2
    end
end

@testset "a missing drug table is reported, not pretended" begin
    r = solve(fixture(nothing); mode = "inert")
    @test r.diagnostics["drug_rule"] == "inert: no drug table"
    @test r.diagnostics["drugs_held"] == 0
    @test fold(r, "T") ≈ 1.0 rtol = 1e-6          # nothing held: same as propagate
    r = solve(fixture(Set{String}()); mode = "inert")
    @test r.diagnostics["drug_rule"] == "inert"     # a table listing none is a real answer
    @test r.diagnostics["drugs_held"] == 0
end

@testset "an explicit observation of a drug wins" begin
    r = solve(fixture(Set(["R-C", "R-D"])); mode = "inert", x = 1.0,
              extra = Dict("C" => (4.0, 1.0)))
    @test fold(r, "C") ≈ 4.0 rtol = 1e-9
    @test fold(r, "T") ≈ 0.25 rtol = 1e-6
    @test r.diagnostics["drugs_held"] == 1          # D only
end

@testset "a confidence-0 observation does not unpin a drug" begin
    # The confidence gate discards it, so it must not stop the hold either.
    r = solve(fixture(Set(["R-C", "R-D"])); mode = "inert", x = 2.0,
              extra = Dict("C" => (4.0, 0.0)))
    @test fold(r, "C") ≈ 1.0 rtol = 1e-9
    @test fold(r, "T") ≈ 2.0 rtol = 1e-6
    @test r.diagnostics["drugs_held"] == 2
end

@testset "the silo-bridge rebuild keeps the drug list" begin
    # S is one entity split into a receiving copy S1 (X -> S1, no out-edge) and a
    # feeding copy S2 (S2 -> T, no in-edge), the shape silo bridges join. Drugs
    # are resolved before that rebuild, so this checks the BEHAVIOUR (still held
    # with bridges on); the rebuild's own drug_stids pass-through is defensive
    # and not read today (a mutation dropping it survives this test).
    nodes = Dict(node.(["X", "D", "C", "W", "T"]))
    nodes["S1"] = DS.NetworkNode("S1", "R-S", "protein", nothing, "S1", BL)
    nodes["S2"] = DS.NetworkNode("S2", "R-S", "protein", nothing, "S2", BL)
    edges = [edge("X", "C", true, "assembly"), edge("D", "C", true, "assembly"),
             edge("X", "T", true, "input"), edge("W", "T", true, "input"),
             edge("C", "T", false, "regulator"; and = false),
             edge("X", "S1", true, "input"), edge("S2", "T", true, "input")]
    net = DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                             Dict{String, Set{String}}(), Set(["R-C", "R-D"]))
    with_env("DS_SILO_BRIDGE_MAX_REACH" => "50") do
        @test !isempty(DS.silo_bridge_edges(net))       # the rebuild path runs
        r = solve(net; mode = "inert")
        @test r.diagnostics["drug_rule"] == "inert"
        @test r.diagnostics["drugs_held"] == 2
        @test fold(r, "C") ≈ 1.0 rtol = 1e-9
    end
end

@testset "parse_complete_network reads drugs.csv on both return paths" begin
    mktempdir() do dir
        write(joinpath(dir, "logic_network.csv"),
              "source_id,target_id,pos_neg,and_or,edge_type,stoichiometry\n" *
              "u-d,u-c,pos,and,assembly,1\nu-x,u-c,pos,and,assembly,1\n")
        write(joinpath(dir, "stid_to_uuid_mapping.csv"),
              "uuid,stable_id\nu-d,R-ALL-9\nu-x,R-HSA-1\nu-c,R-HSA-2\n")
        ln, um = joinpath(dir, "logic_network.csv"), joinpath(dir, "stid_to_uuid_mapping.csv")
        write(joinpath(dir, "drugs.csv"), "stable_id,schema_class,name,reactome_release\nR-ALL-9,ChemicalDrug,d,97\nR-HSA-2,Complex,c,97\n")
        # no cofactors.csv: the first return path
        net = DS.parse_complete_network(ln, um)
        @test net.drug_stids == Set(["R-ALL-9", "R-HSA-2"])
        @test DS.drug_uuids(net) == Set(["u-d", "u-c"])
        # a cofactors.csv declaring one present cofactor: the second return path
        write(joinpath(dir, "cofactors.csv"),
              "stable_id,molecule,chebi_id,name,in_network,reactome_release\nR-HSA-1,ATP,1,x,1,97\n")
        net2 = DS.parse_complete_network(ln, um)
        @test net2.cofactor_stids == Set(["R-HSA-1"])
        @test net2.drug_stids == Set(["R-ALL-9", "R-HSA-2"])
        rm(joinpath(dir, "drugs.csv"))
        @test DS.parse_complete_network(ln, um).drug_stids === nothing
    end
end

@testset "a typo in the mode is an error" begin
    for bad in ("Inert", "inert ", "on", "")
        with_env("DS_DRUG_MODE" => bad) do
            @test_throws ArgumentError DS.drug_mode()
        end
    end
end

@testset "drugs.csv parsing and JSON round trip" begin
    mktempdir() do dir
        ln = joinpath(dir, "logic_network.csv")
        write(ln, "")
        @test DS.parse_drug_list(ln) === nothing
        write(joinpath(dir, "drugs.csv"), "stable_id,schema_class,name,reactome_release\n")
        @test DS.parse_drug_list(ln) == Set{String}()
        write(joinpath(dir, "drugs.csv"),
              "stable_id,schema_class,name,reactome_release\nR-ALL-1,ChemicalDrug,x,97\nR-HSA-2,Complex,y,97\n")
        @test DS.parse_drug_list(ln) == Set(["R-ALL-1", "R-HSA-2"])
        write(joinpath(dir, "drugs.csv"), "id\nR-ALL-1\n")
        @test_throws ArgumentError DS.parse_drug_list(ln)
    end
    @test DS.drug_stids_json(fixture(nothing)) === nothing
    @test DS.drug_stids_json(fixture(Set(["R-D", "R-C"]))) == ["R-C", "R-D"]
    @test DS.drug_stids_from_json(nothing) === nothing
    @test DS.drug_stids_from_json(["R-D"]) == Set(["R-D"])
    @test_throws ArgumentError DS.drug_stids_from_json("R-D")
end

# specs/048: pathogen-derived nodes (pathogens.csv) under DS_PATHOGEN_MODE, the
# same hold on its own switch. Fixture: D is a viral protein, C = X:D inhibits
# the reaction making T (the N:M:PDPK1 shape).
function pfixture(pathogen_stids; drug_stids = nothing)
    f = fixture(drug_stids)
    DS.ReactionNetwork(f.nodes, f.edges, f.set_mappings, f.cofactor_stids, f.containment,
                       f.drug_stids, f.pools, pathogen_stids)
end

function psolve(net; mode, drug = nothing, x = 2.0, extra = Dict{String, Tuple{Float64, Float64}}())
    with_env("DS_PATHOGEN_MODE" => mode) do
        solve(net; mode = drug, x = x, extra = extra)
    end
end

@testset "pathogens (specs/048)" begin
    @testset "default is inert; propagate is byte-identical to a network with no pathogen list" begin
        with_env("DS_PATHOGEN_MODE" => nothing) do
            @test DS.pathogen_mode() == "inert"
        end
        d = psolve(pfixture(Set(["R-C", "R-D"])); mode = nothing)
        @test d.node_activities == psolve(pfixture(Set(["R-C", "R-D"])); mode = "inert").node_activities
        @test d.diagnostics["pathogen_rule"] == "inert"
        n = psolve(pfixture(nothing); mode = nothing)
        @test n.diagnostics["pathogen_rule"] == "inert: no pathogen table"
        @test n.diagnostics["pathogens_held"] == 0
        a = psolve(pfixture(Set(["R-C", "R-D"])); mode = "propagate")
        b = psolve(pfixture(nothing); mode = "propagate")
        @test a.node_activities == b.node_activities
        @test a.diagnostics["pathogen_rule"] == "propagate"
        @test a.diagnostics["pathogens_held"] == 0
        @test fold(a, "T") ≈ 1.0 rtol = 1e-6
        # the pre-existing 7-argument constructor leaves the list absent
        @test fixture(nothing).pathogen_stids === nothing
    end
    @testset "inert holds pathogen-derived nodes at baseline and T follows X" begin
        for x in (2.0, 0.5, 0.0, 50.0)
            r = psolve(pfixture(Set(["R-C", "R-D"])); mode = "inert", x = x)
            @test fold(r, "C") ≈ 1.0 rtol = 1e-9
            @test fold(r, "D") ≈ 1.0 rtol = 1e-9
            @test fold(r, "T") ≈ min(x, 100.0) rtol = 1e-6
            @test r.diagnostics["pathogen_rule"] == "inert"
            @test r.diagnostics["pathogens_held"] == 2
            @test r.diagnostics["drugs_held"] == 0
        end
    end
    @testset "independent of the drug rule" begin
        # a node listed only as a drug is not held by DS_PATHOGEN_MODE, and vice versa
        r = psolve(pfixture(nothing; drug_stids = Set(["R-C"])); mode = "inert")
        @test r.diagnostics["pathogen_rule"] == "inert: no pathogen table"
        @test fold(r, "C") ≈ 2.0 rtol = 1e-6
        r = psolve(pfixture(Set(["R-C"]); drug_stids = nothing); mode = "propagate", drug = "inert")
        @test r.diagnostics["drug_rule"] == "inert: no drug table"
        @test r.diagnostics["pathogens_held"] == 0
        @test fold(r, "C") ≈ 2.0 rtol = 1e-6
    end
    @testset "a missing pathogen table is reported, not pretended" begin
        r = psolve(pfixture(nothing); mode = "inert")
        @test r.diagnostics["pathogen_rule"] == "inert: no pathogen table"
        @test r.diagnostics["pathogens_held"] == 0
        r = psolve(pfixture(Set{String}()); mode = "inert")
        @test r.diagnostics["pathogen_rule"] == "inert"
        @test r.diagnostics["pathogens_held"] == 0
    end
    @testset "an explicit observation wins; a confidence-0 one does not" begin
        r = psolve(pfixture(Set(["R-C", "R-D"])); mode = "inert", x = 1.0, extra = Dict("C" => (4.0, 1.0)))
        @test fold(r, "C") ≈ 4.0 rtol = 1e-9
        @test r.diagnostics["pathogens_held"] == 1
        r = psolve(pfixture(Set(["R-C", "R-D"])); mode = "inert", x = 2.0, extra = Dict("C" => (4.0, 0.0)))
        @test fold(r, "C") ≈ 1.0 rtol = 1e-9
        @test r.diagnostics["pathogens_held"] == 2
    end
    @testset "a typo in the mode is an error" begin
        for bad in ("Inert", "inert ", "on", "")
            with_env("DS_PATHOGEN_MODE" => bad) do
                @test_throws ArgumentError DS.pathogen_mode()
            end
        end
    end
    @testset "pathogens.csv: both parse return paths, the silo-bridge rebuild, JSON" begin
        mktempdir() do dir
            write(joinpath(dir, "logic_network.csv"),
                  "source_id,target_id,pos_neg,and_or,edge_type,stoichiometry\n" *
                  "u-d,u-c,pos,and,assembly,1\nu-x,u-c,pos,and,assembly,1\n")
            write(joinpath(dir, "stid_to_uuid_mapping.csv"),
                  "uuid,stable_id\nu-d,R-COV-9\nu-x,R-HSA-1\nu-c,R-HSA-2\n")
            ln, um = joinpath(dir, "logic_network.csv"), joinpath(dir, "stid_to_uuid_mapping.csv")
            @test DS.parse_pathogen_list(ln) === nothing
            @test DS.parse_complete_network(ln, um).pathogen_stids === nothing
            write(joinpath(dir, "pathogens.csv"), "stable_id,schema_class,name,reactome_release\n")
            @test DS.parse_pathogen_list(ln) == Set{String}()
            write(joinpath(dir, "pathogens.csv"),
                  "stable_id,schema_class,name,reactome_release\nR-COV-9,EWAS,M,97\nR-HSA-2,Complex,N:M:PDPK1,97\n")
            net = DS.parse_complete_network(ln, um)
            @test net.pathogen_stids == Set(["R-COV-9", "R-HSA-2"])
            @test DS.pathogen_uuids(net) == Set(["u-d", "u-c"])
            write(joinpath(dir, "cofactors.csv"),
                  "stable_id,molecule,chebi_id,name,in_network,reactome_release\nR-HSA-1,ATP,1,x,1,97\n")
            @test DS.parse_complete_network(ln, um).pathogen_stids == Set(["R-COV-9", "R-HSA-2"])
            write(joinpath(dir, "pathogens.csv"), "id\nR-COV-9\n")
            @test_throws ArgumentError DS.parse_pathogen_list(ln)
        end
        nodes = Dict(node.(["X", "D", "C", "W", "T"]))
        nodes["S1"] = DS.NetworkNode("S1", "R-S", "protein", nothing, "S1", BL)
        nodes["S2"] = DS.NetworkNode("S2", "R-S", "protein", nothing, "S2", BL)
        f = fixture(nothing)
        net = DS.ReactionNetwork(nodes, vcat(f.edges, [edge("X", "S1", true, "input"), edge("S2", "T", true, "input")]),
                                 f.set_mappings, f.cofactor_stids, f.containment, nothing, nothing, Set(["R-C", "R-D"]))
        with_env("DS_SILO_BRIDGE_MAX_REACH" => "50") do
            @test !isempty(DS.silo_bridge_edges(net))
            r = psolve(net; mode = "inert")
            @test r.diagnostics["pathogens_held"] == 2
            @test fold(r, "C") ≈ 1.0 rtol = 1e-9
        end
        @test DS.pathogen_stids_json(pfixture(nothing)) === nothing
        @test DS.pathogen_stids_json(pfixture(Set(["R-D", "R-C"]))) == ["R-C", "R-D"]
        @test DS.pathogen_stids_from_json(nothing) === nothing
        @test DS.pathogen_stids_from_json(["R-D"]) == Set(["R-D"])
        @test_throws ArgumentError DS.pathogen_stids_from_json("R-D")
    end
end

end
