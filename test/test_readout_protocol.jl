#!/usr/bin/env julia
#
# specs/141: the readout protocol, ported from the bench — the resolution
# ladder and entity_mean of bench/benchmark_vs_mpbiopath.py (specs/105) and
# the korkut total-abundance rule of bench/korkut_bench.py (specs/134,
# abundance_a3). These fixtures COMPUTE the quantity each rule reasons about
# (groups from tables and edges, the carried set from edge walks), never pin
# it; the parity gate against the bench itself is bench/analysis/141/, run on
# the frozen catalog.

using Test

include(joinpath(@__DIR__, "..", "src", "DeltaSignal.jl"))
using .DeltaSignal

edge(p, c; type="input") = LogicNetworkEdge(p, c, false, true, 1.0, type)

pin_index(pairs::Pair...; diagram=Dict{String,String}()) =
    DeltaSignal.PinIndex(Dict(String(k) => String[v...] for (k, v) in pairs), diagram)

readout_index(; proxies=Dict{String,Vector{String}}(), cofactors=Set{String}(),
              has_proxy=true, has_cof=true) =
    DeltaSignal.ReadoutIndex(proxies, cofactors, has_proxy, has_cof)

@testset "readout protocol (specs/141)" begin

    @testset "direct: one group of every variant copy, first candidate wins" begin
        idx = pin_index("SE" => ["E1", "E2"], "SF" => ["F1"])
        r = resolve_readout(["SE", "SF"], String[], idx, readout_index(), LogicNetworkEdge[])
        @test r.basis == "direct"
        @test r.groups == [["E1", "E2"]]
        # The historical second candidate is reached only when the first resolves nothing.
        r2 = resolve_readout(["nope", "SF"], String[], idx, readout_index(), LogicNetworkEdge[])
        @test r2.groups == [["F1"]]
        # entity-mean over one group is the group mean; a missing uuid reads baseline.
        @test readout_entity_mean(r.groups, Dict("E1" => 0.02, "E2" => 0.04)) ≈ 0.03
        @test readout_entity_mean(r.groups, Dict("E1" => 0.03)) ≈ (0.03 + 0.01) / 2
    end

    @testset "proxy: grouped by non-cofactor output entities, not flat" begin
        # SP has no node of its own; three proxy reactions, R1 and R3 making
        # the same entity SX, R2 making SY. Entity-mean weighs SX and SY once
        # each — the flat mean over proxies is a DIFFERENT number, so this
        # fixture fails under the retired rule.
        edges = [edge("R1", "O1"; type="output"), edge("R2", "O2"; type="output"),
                 edge("R3", "O1b"; type="output")]
        idx = pin_index(; diagram=Dict("O1" => "SX", "O1b" => "SX", "O2" => "SY"))
        ridx = readout_index(proxies=Dict("SP" => ["R1", "R2", "R3"]))
        r = resolve_readout(["SP"], String[], idx, ridx, edges)
        @test r.basis == "proxy"
        @test r.groups == [["R1", "R3"], ["R2"]]
        acts = Dict("R1" => 0.01, "R3" => 0.03, "R2" => 0.06)
        @test readout_entity_mean(r.groups, acts) ≈ (0.02 + 0.06) / 2
        @test readout_entity_mean(r.groups, acts) != (0.01 + 0.03 + 0.06) / 3
    end

    @testset "proxy: cofactor outputs do not group; an outputless proxy is its own group" begin
        edges = [edge("R1", "ATP"; type="output"), edge("R2", "O2"; type="output")]
        idx = pin_index(; diagram=Dict("ATP" => "SATP", "O2" => "SY"))
        ridx = readout_index(proxies=Dict("SP" => ["R1", "R2", "R3"]),
                             cofactors=Set(["SATP"]))
        r = resolve_readout(["SP"], String[], idx, ridx, edges)
        # R1's only output is a cofactor and R3 has none: each keys by its own uuid.
        @test r.groups == [["R1"], ["R2"], ["R3"]]
    end

    @testset "set members: one group per member, unequal copy counts weigh once" begin
        idx = pin_index("M1" => ["a"], "M2" => ["b", "c"])
        r = resolve_readout(["SSET"], ["M1", "M2", "M2", "MISSING"], idx,
                            readout_index(), LogicNetworkEdge[])
        @test r.basis == "set_members"
        @test r.groups == [["a"], ["b", "c"]]   # dedup; unresolved member dropped
        acts = Dict("a" => 0.01, "b" => 0.05, "c" => 0.05)
        @test readout_entity_mean(r.groups, acts) ≈ (0.01 + 0.05) / 2
        @test readout_entity_mean(r.groups, acts) != (0.01 + 0.05 + 0.05) / 3
    end

    @testset "nothing resolves: basis none, no value invented" begin
        r = resolve_readout(["nope"], ["also-nope"], pin_index("SE" => ["E1"]),
                            readout_index(), LogicNetworkEdge[])
        @test r.basis == "none"
        @test isempty(r.groups)
        @test readout_entity_mean(r.groups, Dict{String,Float64}()) === nothing
    end

    @testset "the direct step outranks proxies; proxies outrank members" begin
        idx = pin_index("SE" => ["E1"], "M1" => ["a"])
        ridx = readout_index(proxies=Dict("SE" => ["R1"], "SQ" => ["R2"]))
        @test resolve_readout(["SE"], ["M1"], idx, ridx, LogicNetworkEdge[]).basis == "direct"
        @test resolve_readout(["SQ"], ["M1"], idx, ridx, LogicNetworkEdge[]).basis == "proxy"
    end

    @testset "abundance: made-from copies drop, recycle does not make" begin
        # Protein P: starting copy S; complex C assembled from S; released
        # copy Rel dissociated back out of C, recycled to S. Under forms all
        # three are the readout; under abundance only S is (assembly makes C,
        # dissociation makes Rel, and the recycle edge into S is not "made
        # from" — specs/134 amendment 3). The fixture fails under a1/a2
        # semantics, where recycle counted and S dropped too.
        idx = pin_index("SP" => ["S", "C", "Rel"])
        edges = [edge("S", "C"; type="assembly"),
                 edge("C", "Rel"; type="dissociation"),
                 edge("Rel", "S"; type="recycle")]
        occ = resolve_abundance_occurrences(["SP"], idx)
        @test occ == ["S", "C", "Rel"]
        @test abundance_copies(occ, edges) == ["S"]
    end

    @testset "abundance: a carried set makes its reaction's product a non-copy" begin
        # S feeds a set pool SET by set_member (SET carries the protein,
        # amendment 2); reaction R takes SET as input and outputs C2 — C2 is
        # made from the protein through the carried set, so it drops; a
        # reaction product fed only by OTHER inputs stays.
        idx = pin_index("SP" => ["S", "C2", "D"])
        edges = [edge("S", "SET"; type="set_member"),
                 edge("SET", "R"; type="input"), edge("R", "C2"; type="output"),
                 edge("X", "R2"; type="input"), edge("R2", "D"; type="output")]
        @test abundance_copies(["S", "C2", "D"], edges) == ["S", "D"]
    end

    @testset "abundance: the union is over ALL stable ids, deduplicated, order kept" begin
        idx = pin_index("S1" => ["a", "b"], "S2" => ["b", "c"])
        @test resolve_abundance_occurrences(["S1", "S2"], idx) == ["a", "b", "c"]
    end

    @testset "abundance statistic: mean over unpinned; all pinned leaves scoring" begin
        acts = Dict("S" => 0.5, "T" => 0.3)
        @test abundance_value(["S", "T"], Set{String}(), acts) ==
              (status = "scored", value = 0.4)
        @test abundance_value(["S", "T"], Set(["S"]), acts) ==
              (status = "scored", value = 0.3)
        @test abundance_value(["S"], Set(["S"]), acts).status == "trivial"
        @test abundance_value(String[], Set{String}(), acts).status == "no_copy"
        # A kept copy missing from the activities reads baseline, as the bench reads it.
        @test abundance_value(["S", "U"], Set{String}(), acts).value ≈ (0.5 + 0.01) / 2
    end

    @testset "the call boundaries are the bench's: strict <0.85, inclusive >=1.15" begin
        @test readout_call(0.84) == "down"
        @test readout_call(0.85) == "normal"
        @test readout_call(1.1499) == "normal"
        @test readout_call(1.15) == "up"
        @test readout_call(0.0) == "down"
        @test readout_call(100.0) == "up"
    end

    @testset "build_readout_index reads the bundle tables; absence is recorded" begin
        dir = mktempdir()
        write(joinpath(dir, "entity_reaction_proxy_mapping.csv"),
              "entity_stable_id,proxy_uuid,proxy_role\nSP,R1,producing\nSP,R2,consuming\nSQ,R3,producing\n")
        write(joinpath(dir, "cofactors.csv"), "stable_id,name\nSATP,ATP\n")
        ridx = build_readout_index(dir)
        @test ridx.proxies == Dict("SP" => ["R1", "R2"], "SQ" => ["R3"])
        @test ridx.cofactors == Set(["SATP"])
        @test ridx.has_proxy_table && ridx.has_cofactor_table

        empty_dir = mktempdir()
        bare = build_readout_index(empty_dir)
        @test isempty(bare.proxies) && isempty(bare.cofactors)
        @test !bare.has_proxy_table && !bare.has_cofactor_table

        bad = mktempdir()
        write(joinpath(bad, "entity_reaction_proxy_mapping.csv"), "foo,bar\n1,2\n")
        @test_throws ArgumentError build_readout_index(bad)
    end

end
