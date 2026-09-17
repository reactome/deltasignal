#!/usr/bin/env julia
#
# Observation-file reading, which the CLI does before every solve.
#
# `read_observations` was inline in `execute_solve_command` and therefore
# untestable without spawning a process, so none of its behaviour was pinned:
# not the range checks, not the column aliases, and not what happens when a
# node is listed twice.
#
# It was listed twice by silent last-write-wins. The documented input format
# says every extra column -- including `condition` -- is ignored, so a file
# holding several conditions is applied as ONE simultaneous perturbation set.
# That invites exactly this file:
#
#     P53,treatment_A,0.3,0.9
#     P53,treatment_B,2.5,0.9
#
# which gave P53 = 2.5 and discarded the knockdown, with no warning. Two
# contradictory measurements are a conflict, not a set.

using Test

include(joinpath(@__DIR__, "..", "cli", "deltasignal.jl"))

"""Write `content` to a temp CSV, run `f` on the path, always clean up."""
function with_csv(f, content::String)
    path = tempname() * ".csv"
    write(path, content)
    try
        return f(path)
    finally
        rm(path; force=true)
    end
end

@testset "CLI observation reading" begin

    @testset "a well-formed file loads" begin
        obs = with_csv("node,condition,value,confidence\nEGFR,A,1.5,0.9\nMYC,A,2.1,0.8\n") do p
            read_observations(p)
        end
        @test length(obs) == 2
        @test obs["EGFR"] == (1.5, 0.9)
        @test obs["MYC"] == (2.1, 0.8)
    end

    @testset "column aliases are accepted" begin
        # node_uuid/node and activity/value are both documented spellings.
        a = with_csv("node_uuid,activity,confidence\nX,2.0,1.0\n") do p; read_observations(p); end
        b = with_csv("node,value,confidence\nX,2.0,1.0\n") do p; read_observations(p); end
        @test a == b
        @test a["X"] == (2.0, 1.0)
    end

    @testset "a node listed twice with different values is a conflict" begin
        # The defect. Previously P53 came back as 2.5 with the 0.3 discarded.
        conflicting = "node,condition,value,confidence\nP53,A,0.3,0.9\nP53,B,2.5,0.9\n"
        @test_throws ErrorException with_csv(p -> read_observations(p), conflicting)

        err = try
            with_csv(p -> read_observations(p), conflicting)
            ""
        catch e
            sprint(showerror, e)
        end
        @test occursin("P53", err)
        @test occursin("more than once", err)
        # Both values must appear, so the conflict can be found in the file.
        @test occursin("0.3", err)
        @test occursin("2.5", err)

        # A differing CONFIDENCE is a conflict too, not just a differing value.
        @test_throws ErrorException with_csv(p -> read_observations(p),
            "node,value,confidence\nP53,0.3,0.9\nP53,0.3,0.1\n")
    end

    @testset "an exactly repeated row is not a conflict" begin
        # Harmless duplication must not break a file that already worked.
        obs = with_csv("node,condition,value,confidence\nP53,A,0.3,0.9\nP53,B,0.3,0.9\n") do p
            read_observations(p)
        end
        @test obs["P53"] == (0.3, 0.9)
        @test length(obs) == 1
    end

    @testset "out-of-range values fail loudly" begin
        # Activity is 0-100 and is range-checked, not clamped: -100 and 0 both
        # report 0 downstream but produce completely different networks.
        @test_throws ErrorException with_csv(p -> read_observations(p),
            "node,value,confidence\nX,-1,1.0\n")
        @test_throws ErrorException with_csv(p -> read_observations(p),
            "node,value,confidence\nX,101,1.0\n")
        @test_throws ErrorException with_csv(p -> read_observations(p),
            "node,value,confidence\nX,1.0,1.5\n")
        @test_throws ErrorException with_csv(p -> read_observations(p),
            "node,value,confidence\nX,1.0,-0.1\n")

        # The boundaries themselves are legal.
        @test with_csv(p -> read_observations(p),
            "node,value,confidence\nX,0,0\n")["X"] == (0.0, 0.0)
        @test with_csv(p -> read_observations(p),
            "node,value,confidence\nX,100,1\n")["X"] == (100.0, 1.0)
    end

    @testset "extra columns are ignored" begin
        # `condition` in particular: the file is one perturbation set.
        obs = with_csv("node,condition,replicate,value,confidence\nX,treatment,3,2.0,1.0\n") do p
            read_observations(p)
        end
        @test obs == Dict("X" => (2.0, 1.0))
    end
end
