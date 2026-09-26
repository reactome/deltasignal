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

end
