#!/usr/bin/env julia
#
# specs/143: the served scope statement. Two kinds of guard:
#
# 1. The shipped file itself — every number must carry its provenance (build,
#    solver, date, spec), every limitation its spec reference, and no CPPA
#    content may appear (the data agreement keeps its numbers private). A
#    statement failing these cannot ship, the way check_doc_counts.py stops a
#    wrong test table.
# 2. The handler — serves the file verbatim with the served-catalog
#    comparison, and refuses (allow-listed message) on a missing or malformed
#    file rather than fabricating a statement.

using Test

include(joinpath(@__DIR__, "..", "src", "DeltaSignal.jl"))
using .DeltaSignal

const SRV = DeltaSignal
const SCOPE_FILE = joinpath(@__DIR__, "..", "src", "api", "scope.json")

@testset "scope statement (specs/143)" begin

    @testset "the shipped statement carries its provenance" begin
        raw = read(SCOPE_FILE, String)
        st = SRV.JSON3.read(raw)
        @test !isempty(String(st[:version]))
        @test !isempty(String(st[:provenance_build]))
        @test length(st[:summary]) >= 1
        @test length(st[:validated]) >= 3
        for v in st[:validated]
            for k in (:axis, :statement, :value, :catalog_build, :solver_commit, :date, :spec)
                @test haskey(v, k) && !isempty(String(v[k]))
            end
        end
        @test length(st[:limitations]) >= 5
        for l in st[:limitations]
            for k in (:name, :statement, :spec)
                @test haskey(l, k) && !isempty(String(l[k]))
            end
        end
        @test length(st[:protocols]) == 3
        for pr in st[:protocols]
            for k in (:name, :endpoint, :statement, :spec)
                @test haskey(pr, k) && !isempty(String(pr[k]))
            end
        end
        # FR-007: no CPPA content at all — its numbers are private under the
        # data agreement, and the simplest enforceable rule is absence.
        @test !occursin("CPPA", raw)
        # The recorded limitation classes are all covered (FR-003).
        names = Set(String(l[:name]) for l in st[:limitations])
        for required in ("knowledge_ceiling", "drug_response", "cell_context",
                         "magnitude_rails", "switch_pathways", "synthesis_stock",
                         "convergence")
            @test required in names
        end
    end

    @testset "the handler serves the statement with the catalog comparison" begin
        req = SRV.HTTP.Request("GET", "/api/scope")
        mktempdir() do catalog
            # A catalog whose BUILD.json matches the statement's build.
            st = SRV.JSON3.read(read(SCOPE_FILE, String))
            write(joinpath(catalog, "BUILD.json"),
                  """{"build_id": "$(st[:provenance_build])", "status": "complete"}""")
            resp = SRV.scope_handler(req; catalog_dir = catalog)
            @test resp.status == 200
            body = SRV.JSON3.read(String(resp.body))
            @test body[:measured_on_served_catalog] == true
            @test String(body[:served_catalog]) == String(st[:provenance_build])
            # served verbatim: the statement in the response is the file's
            @test body[:statement][:provenance_build] == st[:provenance_build]
            @test length(body[:statement][:validated]) == length(st[:validated])
        end
        mktempdir() do catalog
            write(joinpath(catalog, "BUILD.json"),
                  """{"build_id": "20991231-0000_notthis", "status": "complete"}""")
            body = SRV.JSON3.read(String(SRV.scope_handler(req; catalog_dir = catalog).body))
            @test body[:measured_on_served_catalog] == false
            @test String(body[:served_catalog]) == "20991231-0000_notthis"
        end
        # No usable catalog information: served null, match false — never a guess.
        body = SRV.JSON3.read(String(SRV.scope_handler(req;
            catalog_dir = joinpath(tempdir(), "no-such-catalog-143")).body))
        @test body[:measured_on_served_catalog] == false
        @test body[:served_catalog] === nothing
    end

    @testset "a missing or malformed statement refuses, never fabricates" begin
        req = SRV.HTTP.Request("GET", "/api/scope")
        resp = SRV.scope_handler(req; path = joinpath(tempdir(), "absent-scope-143.json"))
        @test resp.status == 400
        @test String(SRV.JSON3.read(String(resp.body))[:message]) == SRV.SCOPE_STATEMENT_ERROR
        mktempdir() do dir
            bad = joinpath(dir, "scope.json")
            write(bad, "not json {")
            resp = SRV.scope_handler(req; path = bad)
            @test resp.status == 400
            @test String(SRV.JSON3.read(String(resp.body))[:message]) == SRV.SCOPE_STATEMENT_ERROR
            # Valid JSON but no provenance build: equally unusable.
            write(bad, """{"version": "x", "summary": []}""")
            resp = SRV.scope_handler(req; path = bad)
            @test resp.status == 400
        end
    end
end
