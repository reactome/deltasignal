#!/usr/bin/env julia
#
# specs/140: the pin protocol, ported from bench/benchmark_vs_mpbiopath.py's
# `pin_occurrences` (DS_PIN_SCOPE=root_cycle, DS_PIN_ROOT_IGNORE=depletion).
# These fixtures COMPUTE the quantity each rule reasons about (roots from
# edges, reachability from the graph), never pin it; the parity gate against
# the bench itself is bench/analysis/140/, run on the frozen catalog.

using Test

include(joinpath(@__DIR__, "..", "src", "DeltaSignal.jl"))
using .DeltaSignal

edge(p, c; type="input") = LogicNetworkEdge(p, c, false, true, 1.0, type)

# A PinIndex built by hand: stid → uuids in file order, uuid → diagram id.
index(pairs::Pair...; diagram=Dict{String,String}()) =
    DeltaSignal.PinIndex(Dict(String(k) => String[v...] for (k, v) in pairs), diagram)

@testset "pin protocol (specs/140)" begin

    @testset "roots are pinned; downstream copies are computed, not set" begin
        # A (root copy) → R → B (downstream copy of the same leaf).
        edges = [edge("A", "R"), edge("R", "B"; type="output")]
        idx = index("SA" => ["A", "B"]; diagram=Dict("A" => "SA", "B" => "SB"))
        r = resolve_pins(["SA"], idx, edges; root_ignore="depletion")
        @test r.pins == ["A"]
        @test r.basis == "root"
        @test r.occurrences == 2
    end

    @testset "is-or-contains: the root complex containing the gene is pinned" begin
        # The gene has no node of its own; its starting form is complex C.
        edges = [edge("C", "R"), edge("R", "X"; type="output")]
        idx = index("SA" => ["C"]; diagram=Dict("C" => "SC"))
        r = resolve_pins(["SA"], idx, edges)
        @test r.pins == ["C"]
        @test r.basis == "root"
    end

    @testset "an unknown stable id resolves to none" begin
        r = resolve_pins(["nope"], index("SA" => ["A"]), [edge("A", "R")])
        @test r.pins == String[]
        @test r.basis == "none"
        @test r.occurrences == 0
    end

    @testset "a query's stable ids share one occurrence list (dedup)" begin
        edges = [edge("A", "R")]
        idx = index("SA" => ["A"], "SC" => ["A"])
        r = resolve_pins(["SA", "SC"], idx, edges)
        @test r.pins == ["A"]
        @test r.occurrences == 1
    end

    @testset "depletion-fed roots (specs/058) and the dead-end exception" begin
        # D is fed only by depletion and feeds a reaction: a root by default.
        edges = [edge("Z", "D"; type="depletion"), edge("D", "R")]
        idx = index("SD" => ["D"]; diagram=Dict("D" => "SD"))
        r = resolve_pins(["SD"], idx, edges)
        @test r.pins == ["D"] && r.basis == "root"
        # Under DS_PIN_ROOT_IGNORE=none the depletion edge counts, so D is not
        # a root; the fallback accepts it (its one producer is derived).
        r = resolve_pins(["SD"], idx, edges; root_ignore="")
        @test r.pins == ["D"] && r.basis == "root_cycle"
        # A depletion-fed DEAD END is not a root (week review 2026-10-15) and
        # not a fallback candidate (no out-edge): nothing to pin.
        edges2 = [edge("Z", "D2"; type="depletion"), edge("Z", "W")]
        r = resolve_pins(["SD2"], index("SD2" => ["D2"]), edges2)
        @test r.pins == String[] && r.basis == "none" && r.occurrences == 1
    end

    @testset "recycle and context_bridge edges never decide a root" begin
        # E's only in-edge is a recycle edge (specs/081): still a root.
        edges = [edge("W", "E"; type="recycle"), edge("E", "R")]
        @test resolve_pins(["SE"], index("SE" => ["E"]), edges).basis == "root"
        # F's only in-edge is a context bridge (specs/056): still a root.
        edges = [edge("W", "F"; type="context_bridge"), edge("F", "R")]
        @test resolve_pins(["SF"], index("SF" => ["F"]), edges).basis == "root"
        # Ordering subtlety ported from the bench: a context_bridge OUT-edge
        # still counts as an out-edge, so a depletion-fed node with only a
        # bridge out-edge is a root, not a dead end …
        edges = [edge("Z", "G"; type="depletion"), edge("G", "O"; type="context_bridge")]
        @test resolve_pins(["SG"], index("SG" => ["G"]), edges).basis == "root"
        # … while a recycle out-edge counts for nothing, so the same node with
        # only a recycle out-edge is a dead end.
        edges = [edge("Z", "H"; type="depletion"), edge("H", "O"; type="recycle")]
        @test resolve_pins(["SH"], index("SH" => ["H"]), edges).basis == "none"
    end

    @testset "root_cycle fallback: a catalytic cycle's enzyme, own copy preferred" begin
        # S → R1 ← E; R1 → ES → R2 → {E, P}: the cycle regenerates E, so the
        # gene E has no root. Both E and ES qualify (every producer is their
        # own downstream); the gene's OWN entity is preferred over the complex.
        edges = [edge("S", "R1"), edge("E", "R1"; type="catalyst"),
                 edge("R1", "ES"; type="output"), edge("ES", "R2"),
                 edge("R2", "E"; type="output"), edge("R2", "P"; type="output")]
        idx = index("SE" => ["E", "ES"], "SS" => ["S", "ES"];
                    diagram=Dict("E" => "SE", "ES" => "SES", "S" => "SS"))
        r = resolve_pins(["SE"], idx, edges)
        @test r.pins == ["E"]
        @test r.basis == "root_cycle"
        @test r.occurrences == 2
        # S itself is a plain root; the fallback never fires for it.
        @test resolve_pins(["SS"], idx, edges).pins == ["S"]
    end

    @testset "root_cycle fallback: entry filter and mutual regeneration" begin
        # C1 regenerates itself (C1 → R1 → C1) and hands off to C2's own loop
        # by a derived dissociation edge; C2 qualifies too, but C1 reaches C2
        # and not the reverse, so C1 is the entry.
        edges = [edge("C1", "R1"), edge("R1", "C1"; type="output"),
                 edge("R1", "C2"; type="dissociation"),
                 edge("C2", "R2"), edge("R2", "C2"; type="output")]
        idx = index("SX" => ["C1", "C2"]; diagram=Dict("C1" => "SC1", "C2" => "SC2"))
        r = resolve_pins(["SX"], idx, edges)
        @test r.pins == ["C1"]
        @test r.basis == "root_cycle"
        # Mutual regeneration: each is downstream of the other, no entry —
        # the whole cycle's pool is pinned.
        edges = [edge("C1", "R1"), edge("R1", "C2"; type="output"),
                 edge("C2", "R2"), edge("R2", "C1"; type="output")]
        r = resolve_pins(["SX"], idx, edges)
        @test r.pins == ["C1", "C2"]
    end

    @testset "build_pin_index reads nodes.csv the way the bench does" begin
        mktempdir() do d
            path = joinpath(d, "nodes.csv")
            write(path, join([
                "uuid,node_kind,member_leaves,diagram_entity_id",
                "n1,entity,SA|SB,SC",
                "n2,set_pool,SA,SP",
                "n3,entity,SA|SA,SA",
                "n4,reaction,,",
            ], "\n") * "\n")
            idx = build_pin_index(path)
            @test idx !== nothing
            # indexed by diagram id AND every member leaf, file order kept
            @test idx.stid_to_uuids["SA"] == ["n1", "n3"]
            @test idx.stid_to_uuids["SB"] == ["n1"]
            @test idx.stid_to_uuids["SC"] == ["n1"]
            # a set POOL node is never a pin target (specs/033) …
            @test all(!("n2" in v) for v in values(idx.stid_to_uuids))
            # … but its diagram id is still known (the bench's ident covers
            # every row).
            @test idx.diagram_entity["n2"] == "SP"
            @test !haskey(idx.diagram_entity, "n4")
            # a bundle without the table refuses to approximate
            @test build_pin_index(joinpath(d, "missing.csv")) === nothing
            # a file that is not a node table is an error, not an empty index
            bad = joinpath(d, "bad.csv")
            write(bad, "uuid,foo\nx,1\n")
            @test_throws ArgumentError build_pin_index(bad)
        end
    end

    @testset "DS_PIN_ROOT_IGNORE is read as the bench reads it" begin
        withenv("DS_PIN_ROOT_IGNORE" => nothing) do
            @test pin_root_ignore() == "depletion"
        end
        withenv("DS_PIN_ROOT_IGNORE" => "none") do
            @test pin_root_ignore() == ""
        end
        withenv("DS_PIN_ROOT_IGNORE" => "dissociation") do
            @test pin_root_ignore() == "dissociation"
        end
    end
end
