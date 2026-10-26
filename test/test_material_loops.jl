# specs/107: material-closed positive loops. Metabolite pools (DS_METABOLITE_POOLS):
# small molecules converted into one another are one conserved pool; multi-step
# carriers (DS_CARRIER_MULTISTEP): a form carried round a path of several
# reactions and released unchanged reads its supply, not its release.
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
nd(u, id = "R-" * u) = u => DS.NetworkNode(u, id, "x", nothing, u, BL)
ed(s, t, ty) = DS.LogicNetworkEdge(s, t, true, true, 1.0, ty)
fold(r, u) = r.node_activities[u] / BL
params() = DS.SteadyStateParams(1.0, 0.1, 2000, 1e-9, "penalty")
function net(nodes, edges, leaves; refs = nothing, ids = Dict{String, String}())
    DS.ReactionNetwork(Dict(nd(u, get(ids, u, "R-" * u)) for u in nodes), edges,
                       Dict{String, DS.SetExpansionMapping}(), Set{String}(), Dict{String, Set{String}}(),
                       nothing, DS.CyclePool[], nothing, DS.CarrierLoopEdge[], leaves, nothing, refs)
end
solve(n, obs; env...) = with_env(() -> DS.solve_steady_state(n, Dict(k => (v, 1.0) for (k, v) in obs), params()),
                                 [string(k) => v for (k, v) in env]...)

# S supplies A; A <-> B by r1 (catalyst E1) and r2 (catalyst E2); B is consumed into T
function metabolites()
    leaves = Dict("S" => ("simple_entity", ["R-ALL-S"]), "A" => ("simple_entity", ["R-ALL-A"]),
                  "B" => ("simple_entity", ["R-ALL-B"]), "T" => ("simple_entity", ["R-ALL-T"]),
                  "E1" => ("simple_entity", ["R-HSA-E1"]), "E2" => ("simple_entity", ["R-HSA-E2"]),
                  "rS" => ("reaction", String[]), "r1" => ("reaction", String[]),
                  "r2" => ("reaction", String[]), "rT" => ("reaction", String[]))
    edges = [ed("S", "rS", "input"), ed("rS", "A", "output"),
             ed("A", "r1", "input"), ed("r1", "B", "output"), ed("E1", "r1", "catalyst"),
             ed("B", "r2", "input"), ed("r2", "A", "output"), ed("E2", "r2", "catalyst"),
             ed("B", "rT", "input"), ed("rT", "T", "output")]
    net(collect(keys(leaves)), edges, leaves)
end

# carrier F binds X (r1), X is modified in the complex (r2), F is released (r3) with
# Xp, which signals on to T; FX and FXp carry F's protein leaf
function carrier()
    leaves = Dict("F" => ("simple_entity", ["R-HSA-F"]), "X" => ("simple_entity", ["R-HSA-X"]),
                  "FX" => ("simple_complex", ["R-HSA-F", "R-HSA-X"]),
                  "FXp" => ("simple_complex", ["R-HSA-F", "R-HSA-Xp"]),
                  "Xp" => ("simple_entity", ["R-HSA-Xp"]), "T" => ("simple_entity", ["R-HSA-T"]),
                  "r1" => ("reaction", String[]), "r2" => ("reaction", String[]),
                  "r3" => ("reaction", String[]), "r4" => ("reaction", String[]))
    edges = [ed("F", "r1", "input"), ed("X", "r1", "input"), ed("r1", "FX", "output"),
             ed("FX", "r2", "input"), ed("r2", "FXp", "output"),
             ed("FXp", "r3", "input"), ed("r3", "F", "output"), ed("r3", "Xp", "output"),
             ed("Xp", "r4", "input"), ed("r4", "T", "output")]
    net(collect(keys(leaves)), edges, leaves)
end

# a small molecule P lipidates protein L (r1) and is released by de-lipidation (r2)
function lipid()
    leaves = Dict("P" => ("simple_entity", ["R-ALL-P"]), "L" => ("simple_entity", ["R-HSA-L"]),
                  "LP" => ("simple_complex", ["R-ALL-P", "R-HSA-L"]), "L2" => ("simple_entity", ["R-HSA-L"]),
                  "r1" => ("reaction", String[]), "r2" => ("reaction", String[]))
    edges = [ed("P", "r1", "input"), ed("L", "r1", "input"), ed("r1", "LP", "output"),
             ed("LP", "r2", "input"), ed("r2", "P", "output"), ed("r2", "L2", "output")]
    net(collect(keys(leaves)), edges, leaves)
end

# GRIP2-like: G at the membrane (Gm) is endocytosed to the cytosol (Gc) and trafficked back;
# the two copies are different leaves of one reference entity
function location(; refs = nothing)
    leaves = Dict("Gm" => ("simple_entity", ["R-HSA-Gm"]), "Gc" => ("simple_entity", ["R-HSA-Gc"]),
                  "Y" => ("simple_entity", ["R-HSA-Y"]),
                  "r1" => ("reaction", String[]), "r2" => ("reaction", String[]))
    edges = [ed("Gm", "r1", "input"), ed("Y", "r1", "input"), ed("r1", "Gc", "output"),
             ed("Gc", "r2", "input"), ed("r2", "Gm", "output")]
    net(collect(keys(leaves)), edges, leaves; refs = refs)
end

@testset "material-closed loops (specs/107)" begin

@testset "defaults are off, and off adds nothing" begin
    with_env("DS_METABOLITE_POOLS" => nothing, "DS_CARRIER_MULTISTEP" => nothing) do
        @test DS.metabolite_pools_mode() == "off"
        @test DS.carrier_multistep_mode() == "off"
        n = carrier()
        m, rep = DS.apply_material_loops(n)
        @test m === n
        @test rep["metabolite_pools_rule"] == "off" && rep["carrier_multistep_rule"] == "off"
    end
    @test_throws ArgumentError with_env(DS.metabolite_pools_mode, "DS_METABOLITE_POOLS" => "yes")
    @test_throws ArgumentError with_env(DS.carrier_multistep_mode, "DS_CARRIER_MULTISTEP" => "loose")
    # the solve reports the rules
    r = solve(carrier(), Dict("X" => 2.0))
    @test r.diagnostics["metabolite_pools_rule"] == "off"
    @test r.diagnostics["carrier_multistep_rule"] == "off"
end

@testset "rule 1: a reversible metabolite pair is one pool" begin
    ps = DS.metabolite_pools(metabolites())
    @test length(ps) == 1
    @test ps[1].forms == ["A", "B"] && ps[1].bases == ["A", "B"]
    @test sort([(q.from, q.to, q.enzyme) for q in ps[1].paths]) == [("A", "B", true), ("B", "A", true)]
    # a protein is a simple_entity too, but not a small molecule: no pool
    n = metabolites()
    n.node_leaves["A"] = ("simple_entity", ["R-HSA-A"])
    @test isempty(DS.metabolite_pools(n))
    # a second non-held input makes it a two-substrate reaction, not a conversion
    n = metabolites()
    n.node_leaves["Q"] = ("simple_entity", ["R-ALL-Q"]); n.nodes["Q"] = last(nd("Q"))
    push!(n.edges, ed("Q", "r1", "input"))
    @test isempty(DS.metabolite_pools(n))
    # without a node table nothing is built, and the solve says so
    n = metabolites()
    m = DS.ReactionNetwork(n.nodes, n.edges, n.set_mappings)
    @test isempty(DS.metabolite_pools(m))
    @test last(DS.apply_material_loops(m; mpools = "on", multistep = "off"))["metabolite_pools_rule"] == "on: no node table"
end

@testset "rule 1: the pool follows supply, and the enzymes set the split" begin
    on(obs) = solve(metabolites(), obs; DS_METABOLITE_POOLS = "on")
    r = on(Dict("S" => 2.0))
    @test r.diagnostics["metabolite_pools"] == 1
    @test fold(r, "A") ≈ 2.0 rtol = 1e-6
    @test fold(r, "B") ≈ 2.0 rtol = 1e-6
    # E1 knocked out: nothing converts A to B
    r = on(Dict("E1" => 0.0))
    @test fold(r, "A") ≈ 2.0 rtol = 1e-4
    @test fold(r, "B") < 1e-3
    @test fold(r, "T") < 1e-3
    # E1 doubled: B : A = 2 : 1 at a fixed total
    r = on(Dict("E1" => 2.0))
    @test fold(r, "A") ≈ 2 / 3 rtol = 1e-6
    @test fold(r, "B") ≈ 4 / 3 rtol = 1e-6
    # at rest every form is at baseline
    r = on(Dict{String, Float64}())
    @test fold(r, "A") ≈ 1.0 rtol = 1e-9
    @test fold(r, "B") ≈ 1.0 rtol = 1e-9
end

@testset "rule 2: a carrier released after several steps reads its supply" begin
    rows, _ = DS.multistep_carriers(carrier(), "strict")
    @test rows == [("F", "r3", true)]
    # the iterated loop has gain X: it rails at 2x
    off = solve(carrier(), Dict("X" => 2.0))
    @test fold(off, "F") > 50
    r = solve(carrier(), Dict("X" => 2.0); DS_CARRIER_MULTISTEP = "strict")
    @test r.diagnostics["carrier_multistep_rows"] == 1
    @test r.diagnostics["carrier_multistep_forms"] == 1
    @test fold(r, "F") ≈ 1.0 rtol = 1e-9
    @test fold(r, "T") ≈ 2.0 rtol = 1e-6
    # a form the generator already lists is not added again
    n = carrier(); push!(n.carrier_loops, ("F", "r3", true))
    @test isempty(first(DS.multistep_carriers(n, "strict")))
    # an outside producer makes the loop open
    n = carrier()
    n.node_leaves["S"] = ("simple_entity", ["R-HSA-S"]); n.nodes["S"] = last(nd("S"))
    push!(n.edges, ed("S", "F", "output"))
    @test first(DS.multistep_carriers(n, "strict")) == [("F", "r3", false)]
end

@testset "rule 3: a small-molecule carrier, and copies across compartments" begin
    @test isempty(first(DS.multistep_carriers(lipid(), "strict")))
    @test first(DS.multistep_carriers(lipid(), "wide")) == [("P", "r2", true)]
    # Gm and Gc are different leaves: only the reference table joins them
    @test isempty(first(DS.multistep_carriers(location(), "strict")))
    rows, with_refs = DS.multistep_carriers(location(), "wide")
    @test isempty(rows) && !with_refs
    @test last(DS.apply_material_loops(location(); mpools = "off", multistep = "wide"))["carrier_multistep_rule"] ==
          "wide: no reference table"
    refs = Dict{String, DS.LeafRef}("R-HSA-Gm" => ("R-HSA-G", nothing), "R-HSA-Gc" => ("R-HSA-G", nothing))
    rows, with_refs = DS.multistep_carriers(location(; refs = refs), "wide")
    @test with_refs
    @test rows == [("Gc", "r1", true), ("Gm", "r2", true)]
    # strict ignores the table
    @test isempty(first(DS.multistep_carriers(location(; refs = refs), "strict")))
end

@testset "rule 3b: compartment copies join, modified forms and free ubiquitin do not" begin
    same = Dict{String, DS.LeafRef}("R-HSA-Gm" => ("R-HSA-G", ""), "R-HSA-Gc" => ("R-HSA-G", ""))
    @test first(DS.multistep_carriers(location(; refs = same), "compartment")) == [("Gc", "r1", true), ("Gm", "r2", true)]
    # Gc is a phosphorylated form (SMAD2 -> p-SMAD2): wide joins them, compartment does not
    modified = Dict{String, DS.LeafRef}("R-HSA-Gm" => ("R-HSA-G", ""), "R-HSA-Gc" => ("R-HSA-G", "O-phospho-L-serine at 465"))
    @test first(DS.multistep_carriers(location(; refs = modified), "wide")) == [("Gc", "r1", true), ("Gm", "r2", true)]
    rows, with_refs = DS.multistep_carriers(location(; refs = modified), "compartment")
    @test isempty(rows) && with_refs
    # a table without the modifications column joins nothing, and the solve says so
    old = Dict{String, DS.LeafRef}("R-HSA-Gm" => ("R-HSA-G", nothing), "R-HSA-Gc" => ("R-HSA-G", nothing))
    rows, with_refs = DS.multistep_carriers(location(; refs = old), "compartment")
    @test isempty(rows) && !with_refs
    @test last(DS.apply_material_loops(location(; refs = old); mpools = "off", multistep = "compartment"))["carrier_multistep_rule"] ==
          "compartment: no modification column"
    @test last(DS.apply_material_loops(location(); mpools = "off", multistep = "compartment"))["carrier_multistep_rule"] ==
          "compartment: no reference table"
    # free ubiquitin (Ub [nucleoplasm]) is a carrier under wide, never under compartment
    n = carrier(); n.nodes["F"] = last(nd("F", "R-HSA-68524"))
    @test first(DS.multistep_carriers(n, "wide")) == [("F", "r3", true)]
    @test isempty(first(DS.multistep_carriers(n, "compartment")))
    # the protein carrier is unaffected
    @test first(DS.multistep_carriers(carrier(), "compartment")) == [("F", "r3", true)]
end

@testset "leaf_refs.csv: parsed when present, absent otherwise" begin
    mktempdir() do d
        f = joinpath(d, "logic_network.csv")
        @test DS.parse_leaf_refs(f) === nothing
        write(joinpath(d, "leaf_refs.csv"), "leaf_stid,reference_stid\nR-HSA-1,R-HSA-9\nR-ALL-2,R-ALL-2\n")
        @test DS.parse_leaf_refs(f) == Dict("R-HSA-1" => ("R-HSA-9", nothing), "R-ALL-2" => ("R-ALL-2", nothing))
        # amendment 1: modifications, "" for none
        write(joinpath(d, "leaf_refs.csv"),
              "leaf_stid,reference_stid,modifications\nR-HSA-1,R-HSA-9,\nR-HSA-3,R-HSA-9,\"O-phospho-L-serine at 465; O-phospho-L-serine at 467\"\n")
        @test DS.parse_leaf_refs(f) == Dict("R-HSA-1" => ("R-HSA-9", ""),
                                            "R-HSA-3" => ("R-HSA-9", "O-phospho-L-serine at 465; O-phospho-L-serine at 467"))
        write(joinpath(d, "leaf_refs.csv"), "leaf,ref\nR-HSA-1,R-HSA-9\n")
        @test_throws ArgumentError DS.parse_leaf_refs(f)
    end
end

end
