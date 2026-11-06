#!/usr/bin/env julia
#
# The API's error sanitizer, and the maintenance hazard inside it.
#
# `user_facing_error` deliberately does NOT echo arbitrary client input back
# into a response body, so it allow-lists the messages it is willing to
# forward and replaces everything else with a generic "Invalid input.".
#
# The hazard: adding a new `*_ERROR` constant and throwing it is not enough.
# If it is missing from the allow-list, the specific message is thrown,
# discarded, and silently downgraded. OBS_UNKNOWN_ERROR was stranded that way
# from the moment it was written. It names the catalog-mismatch failure mode
# -- observations whose uuids come from a different catalog build than the
# network -- which is the hardest failure in this project to diagnose from
# outside, and the client only ever saw "Invalid input.".
#
# The first test below is the guard: it enumerates the constants rather than
# listing them by hand, so a newly added one fails here until it is forwarded.

using Test

include(joinpath(@__DIR__, "..", "src", "DeltaSignal.jl"))
using .DeltaSignal

const SRV = DeltaSignal

@testset "API error sanitization" begin

    @testset "every error constant reaches the client" begin
        # Enumerated from the module, NOT hand-listed: a constant added later
        # is covered automatically.
        err_names = [n for n in Base.names(SRV; all=true)
                 if endswith(String(n), "_ERROR") && isdefined(SRV, n)]
        @test !isempty(err_names)

        for n in err_names
            msg = getfield(SRV, n)
            msg isa String || continue
            status, body = SRV.user_facing_error(ArgumentError(msg))
            @test status == 400
            @test body == msg
        end
    end

    @testset "the stranded one specifically" begin
        # Regression: this returned "Invalid input." and the catalog-mismatch
        # hint never reached anyone.
        _, body = SRV.user_facing_error(ArgumentError(SRV.OBS_UNKNOWN_ERROR))
        @test body == SRV.OBS_UNKNOWN_ERROR
        @test occursin("same catalog build", body)
    end

    @testset "arbitrary messages are still not echoed" begin
        # The reason the allow-list exists. Client-controlled text must not be
        # reflected into the response body.
        status, body = SRV.user_facing_error(ArgumentError("<script>alert(1)</script>"))
        @test status == 400
        @test body == "Invalid input."
        @test !occursin("script", body)

        _, body2 = SRV.user_facing_error(ArgumentError("pathway R-HSA-999999 is bad"))
        @test !occursin("R-HSA-999999", body2)
    end

    @testset "non-ArgumentError exceptions map to their own messages" begin
        @test SRV.user_facing_error(KeyError(:missing))[2] == "Required field missing from input."
        @test SRV.user_facing_error(SystemError("open"))[2] == "Could not read input file."
        for e in (KeyError(:x), SystemError("open"), ArgumentError("x"))
            @test SRV.user_facing_error(e)[1] == 400
        end
    end

    # The dev API served a two-month-old catalog while every benchmark ran
    # against newer ones, and nothing said so. /api/health now reports which
    # catalog is being served, and must not look healthy while serving nothing.
    @testset "health reports the catalog it is serving" begin
        mktempdir() do root
            # a recorded, complete build
            built = joinpath(root, "built")
            for p in ("R-HSA-1", "R-HSA-2")
                mkpath(joinpath(built, p))
                write(joinpath(built, p, "logic_network.csv"), "source_id,target_id\n")
            end
            write(joinpath(built, "BUILD.json"),
                  """{"build_id":"20260925-1200_abc1234","status":"complete",
                      "generator_commit":"abc1234","pathways_built":2}""")
            info = SRV.catalog_provenance(built)
            @test info["build_id"] == "20260925-1200_abc1234"
            @test info["generator_commit"] == "abc1234"
            @test info["status"] == "complete"
            @test info["pathways"] == 2
            @test !haskey(info, "warning")

            # a catalog with no manifest is reported as such, never guessed at
            legacy = joinpath(root, "legacy")
            mkpath(joinpath(legacy, "R-HSA-9"))
            write(joinpath(legacy, "R-HSA-9", "logic_network.csv"), "source_id,target_id\n")
            linfo = SRV.catalog_provenance(legacy)
            @test linfo["build_id"] == "unrecorded"
            @test linfo["pathways"] == 1

            # docker creates an EMPTY directory when a mount path is missing:
            # that must carry a warning, not read as a healthy catalog
            empty = joinpath(root, "empty"); mkpath(empty)
            einfo = SRV.catalog_provenance(empty)
            @test einfo["pathways"] == 0
            @test haskey(einfo, "warning")

            # a pathway directory without a network does not count as served
            partial = joinpath(root, "partial"); mkpath(joinpath(partial, "R-HSA-5"))
            @test SRV.catalog_provenance(partial)["pathways"] == 0

            # a path that does not exist at all
            @test SRV.catalog_provenance(joinpath(root, "nope"))["pathways"] == 0

            # a truncated BUILD.json is reported, not raised (kills the mutant
            # that removes the JSON try/catch)
            trunc = joinpath(root, "trunc")
            mkpath(joinpath(trunc, "R-HSA-3"))
            write(joinpath(trunc, "R-HSA-3", "logic_network.csv"), "source_id,target_id\n")
            write(joinpath(trunc, "BUILD.json"), """{"build_id": "20260925-12""")
            @test SRV.catalog_provenance(trunc)["build_id"] == "unreadable BUILD.json"

            # a self-referencing symlink makes isfile throw ELOOP; health must
            # still answer rather than return 500
            loopy = joinpath(root, "loopy"); mkpath(loopy)
            symlink("R-HSA-7", joinpath(loopy, "R-HSA-7"))
            linfo2 = SRV.catalog_provenance(loopy)
            @test linfo2 isa Dict
            @test haskey(linfo2, "warning")
        end
    end

    # The whole point is that the RUNNING service reports its catalog, so test
    # the handler itself, not only the helper (kills the mutant that drops the
    # "catalog" key from the response).
    # Duplicate keys in the observations object: JSON allows them, and the
    # handler used to keep the last. Conflict is a 400 naming the problem; an
    # exact repeat is fine (it reaches the membership check instead).
    @testset "/api/solve rejects conflicting observations for one node" begin
        req(b) = SRV.HTTP.Request("POST", "/api/solve", ["Content-Type" => "application/json"], b)
        resp = SRV.solve_handler(req("""{"observations": {"X": [80.0, 1.0], "X": [0.0, 1.0]}}"""))
        @test resp.status == 400
        @test SRV.JSON3.read(String(resp.body))[:message] == SRV.OBS_CONFLICT_ERROR
        resp = SRV.solve_handler(req("""{"observations": {"X": [80.0, 1.0], "X": [80.0, 1.0]}}"""))
        @test resp.status == 400
        @test SRV.JSON3.read(String(resp.body))[:message] == SRV.OBS_UNKNOWN_ERROR
    end

    # specs/052: the solve response carries the zero-supply report, so an arm
    # can show which model ran (the first arm recorded "unknown": the field
    # was computed but never returned).
    @testset "/api/solve reports zero_supply_rule" begin
        post(path, b) = SRV.HTTP.Request("POST", path, ["Content-Type" => "application/json"], b)
        parsed = SRV.JSON3.read(String(SRV.parse_handler(post("/api/parse", "")).body))
        u = String(parsed[:nodes][1][:uuid])
        body = """{"network_id": "$(parsed[:network_id])", "observations": {"$u": [0.0, 1.0]}}"""
        for (mode, want) in ((nothing, "off"), ("hold", "hold: no node table"))
            old = get(ENV, "DS_ZERO_SUPPLY_MODE", nothing)
            mode === nothing ? delete!(ENV, "DS_ZERO_SUPPLY_MODE") : (ENV["DS_ZERO_SUPPLY_MODE"] = mode)
            try
                r = SRV.JSON3.read(String(SRV.solve_handler(post("/api/solve", body)).body))
                @test r[:zero_supply_rule] == want
                @test r[:zero_supply_held] == 0 && r[:zero_supply_leaves] == 0
            finally
                old === nothing ? delete!(ENV, "DS_ZERO_SUPPLY_MODE") : (ENV["DS_ZERO_SUPPLY_MODE"] = old)
            end
        end
    end

    # specs/122: edge_redirects on a cached network solve exactly as the same network POSTed with each cut edge
    # taking its input from a block node pinned at 0; an edge the network does not have is refused.
    @testset "/api/solve edge_redirects" begin
        post(path, b) = SRV.HTTP.Request("POST", path, ["Content-Type" => "application/json"], b)
        parsed = SRV.JSON3.read(String(SRV.parse_handler(post("/api/parse", "")).body))
        e = parsed[:edges][1]
        p, c = String(e[:parent_uuid]), String(e[:child_uuid])
        u = String(parsed[:nodes][end][:uuid])
        body = """{"network_id": "$(parsed[:network_id])", "observations": {"$u": [50.0, 1.0]},
                   "edge_redirects": [["$p", "$c"]]}"""
        r = SRV.JSON3.read(String(SRV.solve_handler(post("/api/solve", body)).body))
        @test r[:status] == "success"
        @test r[:edge_redirects] == 1
        b = "drugblock:$p:$c"
        src = first(n for n in parsed[:nodes] if String(n[:uuid]) == p)
        nodes = [Dict(String(k) => v for (k, v) in pairs(n)) for n in parsed[:nodes]]
        push!(nodes, Dict("uuid" => b, "name" => "drug block", "reactome_id" => nothing,
                          "entity_type" => src[:entity_type], "baseline" => src[:baseline], "set_id" => nothing))
        edges = [Dict(String(k) => v for (k, v) in pairs(x)) for x in parsed[:edges]]
        edges[1]["parent_uuid"] = b
        inline = SRV.JSON3.write(Dict("network" => Dict("nodes" => nodes, "edges" => edges,
                                                          "pathways" => parsed[:pathways]),
                                      "observations" => Dict(u => [50.0, 1.0], b => [0.0, 1.0])))
        q = SRV.JSON3.read(String(SRV.solve_handler(post("/api/solve", inline)).body))
        @test q[:status] == "success"
        @test Dict(r[:node_activities]) == Dict(q[:node_activities])
        bad = """{"network_id": "$(parsed[:network_id])", "observations": {"$u": [50.0, 1.0]},
                  "edge_redirects": [["no-such-node", "$c"]]}"""
        resp = SRV.solve_handler(post("/api/solve", bad))
        @test resp.status == 400
        @test occursin("edge_redirects", String(resp.body))
    end

    # specs/044: an uploaded network is saved under its part name, so the
    # sibling rule cannot find boundary_edges.csv; it is its own part.
    @testset "/api/parse reads an uploaded boundary_edges part" begin
        hdr = "source_id,target_id,pos_neg,and_or,edge_type,stoichiometry\n"
        part(name, text) = SRV.HTTP.Multipart("$(name).csv", IOBuffer(text), "text/csv")
        # fresh parts per request: a Multipart's buffer is consumed when sent
        fields() = Dict(
            "logic_network" => part("logic_network", hdr * "u-x,u-r,pos,and,input,1\nu-r,u-y,pos,,output,1\n"),
            "uuid_mapping" => part("map", "uuid,stable_id\nu-a,R-HSA-1\nu-x,R-HSA-2\nu-r,R-HSA-3\nu-y,R-HSA-4\n"))
        function upload(fields)
            form = SRV.HTTP.Form(fields)
            SRV.HTTP.Request("POST", "/api/parse",
                ["Content-Type" => "multipart/form-data; boundary=$(form.boundary)"], read(form))
        end
        nedges(resp) = length(SRV.JSON3.read(String(resp.body))[:edges])
        r = SRV.parse_handler(upload(fields()))
        @test r.status == 200
        @test nedges(r) == 2
        f = fields()
        f["boundary_edges"] = part("boundary_edges", hdr * "u-a,u-x,pos,and,assembly,1\n")
        r = SRV.parse_handler(upload(f))
        @test r.status == 200
        @test nedges(r) == 3
        # week review 2026-10-15: bundle tables come as their own parts and are
        # read as the catalog's are; without them a rule reports no table
        drugs(resp) = SRV.JSON3.read(String(resp.body))[:drug_stids]
        @test drugs(SRV.parse_handler(upload(fields()))) === nothing
        f = fields()
        f["drugs"] = part("drugs", "stable_id,schema_class,name,reactome_release\nR-HSA-2,ChemicalDrug,x,97\n")
        r = SRV.parse_handler(upload(f))
        @test r.status == 200
        @test collect(drugs(r)) == ["R-HSA-2"]
    end

    # specs/140: /api/pins resolves a perturbation the way the benchmarks did.
    @testset "/api/pins resolves the protocol of record" begin
        hdr = "source_id,target_id,pos_neg,and_or,edge_type,stoichiometry\n"
        part(name, text) = SRV.HTTP.Multipart("$(name).csv", IOBuffer(text), "text/csv")
        # u-x (root, leaf R-HSA-2) → u-r → u-y (downstream copy of the SAME leaf)
        fields() = Dict(
            "logic_network" => part("logic_network", hdr * "u-x,u-r,pos,and,input,1\nu-r,u-y,pos,,output,1\n"),
            "uuid_mapping" => part("map", "uuid,stable_id\nu-x,R-HSA-2\nu-r,R-HSA-3\nu-y,R-HSA-4\n"))
        nodes_part() = part("nodes", "uuid,node_kind,member_leaves,diagram_entity_id\n" *
            "u-x,entity,R-HSA-2,R-HSA-2\nu-r,reaction,,\nu-y,entity,R-HSA-2,R-HSA-4\n")
        function upload(fields)
            form = SRV.HTTP.Form(fields)
            SRV.HTTP.Request("POST", "/api/parse",
                ["Content-Type" => "multipart/form-data; boundary=$(form.boundary)"], read(form))
        end
        post(b) = SRV.HTTP.Request("POST", "/api/pins", ["Content-Type" => "application/json"], b)

        f = fields(); f["nodes"] = nodes_part()
        parsed = SRV.JSON3.read(String(SRV.parse_handler(upload(f)).body))
        nid = String(parsed[:network_id])

        resp = SRV.pins_handler(post("""{"network_id": "$nid",
            "queries": [{"id": "q1", "stable_ids": ["R-HSA-2"]},
                        {"id": "q2", "stable_ids": ["R-HSA-999"]}]}"""))
        @test resp.status == 200
        body = SRV.JSON3.read(String(resp.body))
        # the protocol is on the record, not inferred
        @test body[:protocol][:pin_scope] == "root_cycle"
        @test body[:protocol][:root_ignore] == "depletion"
        r1 = body[:results][1]
        # only the ROOT copy of the leaf is pinned; the downstream copy of the
        # same leaf is an occurrence that must be computed, not set
        @test [String(p[:uuid]) for p in r1[:pins]] == ["u-x"]
        @test r1[:basis] == "root"
        @test r1[:occurrences] == 2
        r2 = body[:results][2]
        @test isempty(r2[:pins]) && r2[:basis] == "none" && r2[:occurrences] == 0

        # a network parsed WITHOUT a node table refuses: resolving by exact-id
        # matching would be a second, unvalidated protocol
        parsed2 = SRV.JSON3.read(String(SRV.parse_handler(upload(fields())).body))
        resp = SRV.pins_handler(post("""{"network_id": "$(parsed2[:network_id])",
            "queries": [{"stable_ids": ["R-HSA-2"]}]}"""))
        @test resp.status == 400
        @test String(SRV.JSON3.read(String(resp.body))[:message]) == SRV.PINS_NO_TABLE_ERROR

        # an unknown network_id is the existing re-parse 400, not a fall-through
        resp = SRV.pins_handler(post("""{"network_id": "pw:nope",
            "queries": [{"stable_ids": ["R-HSA-2"]}]}"""))
        @test resp.status == 400
        @test occursin("network_id", String(SRV.JSON3.read(String(resp.body))[:message]))

        # malformed requests name the shape they need
        for bad in ("""{"network_id": "$nid"}""",
                    """{"network_id": "$nid", "queries": []}""",
                    """{"network_id": "$nid", "queries": [{"stable_ids": "R-HSA-2"}]}""",
                    """{"network_id": "$nid", "queries": [{"stable_ids": [2]}]}""")
            resp = SRV.pins_handler(post(bad))
            @test resp.status == 400
            @test String(SRV.JSON3.read(String(resp.body))[:message]) == SRV.PINS_QUERY_ERROR
        end
    end

    # specs/141: /api/readouts reads an entity the way the benchmarks did.
    @testset "/api/readouts resolves and aggregates the protocol of record" begin
        hdr = "source_id,target_id,pos_neg,and_or,edge_type,stoichiometry\n"
        part(name, text) = SRV.HTTP.Multipart("$(name).csv", IOBuffer(text), "text/csv")
        # u-x (leaf R-HSA-2) → u-r → u-y (downstream copy of the same leaf,
        # diagram entity R-HSA-4)
        fields() = Dict(
            "logic_network" => part("logic_network", hdr * "u-x,u-r,pos,and,input,1\nu-r,u-y,pos,,output,1\n"),
            "uuid_mapping" => part("map", "uuid,stable_id\nu-x,R-HSA-2\nu-r,R-HSA-3\nu-y,R-HSA-4\n"))
        nodes_part() = part("nodes", "uuid,node_kind,member_leaves,diagram_entity_id\n" *
            "u-x,entity,R-HSA-2,R-HSA-2\nu-r,reaction,,\nu-y,entity,R-HSA-2,R-HSA-4\n")
        proxies_part() = part("proxies",
            "entity_stable_id,proxy_uuid,proxy_role\nR-HSA-9,u-r,producing\n")
        cofactors_part() = part("cofactors", "stable_id,name\nR-ALL-ATP,ATP\n")
        function upload(fields)
            form = SRV.HTTP.Form(fields)
            SRV.HTTP.Request("POST", "/api/parse",
                ["Content-Type" => "multipart/form-data; boundary=$(form.boundary)"], read(form))
        end
        post(b) = SRV.HTTP.Request("POST", "/api/readouts", ["Content-Type" => "application/json"], b)

        f = fields(); f["nodes"] = nodes_part(); f["proxies"] = proxies_part()
        f["cofactors"] = cofactors_part()
        parsed = SRV.JSON3.read(String(SRV.parse_handler(upload(f)).body))
        nid = String(parsed[:network_id])

        resp = SRV.readouts_handler(post("""{"network_id": "$nid",
            "activities": {"u-x": 0.02, "u-y": 0.04},
            "queries": [{"id": "direct", "stable_ids": ["R-HSA-2"]},
                        {"id": "proxy", "stable_ids": ["R-HSA-9"]},
                        {"id": "set", "stable_ids": ["R-HSA-888"], "members": ["R-HSA-4"]},
                        {"id": "nope", "stable_ids": ["R-HSA-999"]}]}"""))
        @test resp.status == 200
        body = SRV.JSON3.read(String(resp.body))
        # the protocol is on the record, not inferred
        @test body[:protocol][:ko_agg] == "mean_entity"
        @test body[:protocol][:total_readout] == "abundance_a3"
        @test isempty(body[:missing_tables])
        r1 = body[:results][1]
        # one curated-entity group of both variant copies; entity-mean value
        @test [[String(u) for u in g] for g in r1[:groups]] == [["u-x", "u-y"]]
        @test r1[:basis] == "direct"
        @test r1[:value] ≈ 0.03
        @test r1[:fold] ≈ 3.0
        @test r1[:call] == "up"
        r2 = body[:results][2]
        # no node of its own: the producing reaction's flux is the proxy
        @test r2[:basis] == "proxy"
        @test [[String(u) for u in g] for g in r2[:groups]] == [["u-r"]]
        r3 = body[:results][3]
        @test r3[:basis] == "set_members"
        @test [[String(u) for u in g] for g in r3[:groups]] == [["u-y"]]
        @test r3[:value] ≈ 0.04
        r4 = body[:results][4]
        @test r4[:basis] == "none" && isempty(r4[:groups]) && r4[:value] === nothing && r4[:call] === nothing

        # abundance: the downstream copy is made from the pinned starting
        # copy's reaction, so only the starting copy is total protein; pinning
        # it leaves the readout trivial, never scored back
        resp = SRV.readouts_handler(post("""{"network_id": "$nid",
            "activities": {"u-x": 0.5, "u-y": 0.9}, "pinned": ["u-x"],
            "queries": [{"id": "total", "stable_ids": ["R-HSA-2"], "mode": "abundance"}]}"""))
        @test resp.status == 200
        ra = SRV.JSON3.read(String(resp.body))[:results][1]
        @test ra[:occurrences] == 2
        @test [[String(u) for u in g] for g in ra[:groups]] == [["u-x"]]
        @test [String(u) for u in ra[:pinned]] == ["u-x"]
        @test ra[:status] == "trivial" && ra[:value] === nothing && ra[:call] === nothing
        resp = SRV.readouts_handler(post("""{"network_id": "$nid",
            "activities": {"u-x": 0.5},
            "queries": [{"id": "total", "stable_ids": ["R-HSA-2"], "mode": "abundance"}]}"""))
        ra = SRV.JSON3.read(String(resp.body))[:results][1]
        @test ra[:status] == "scored" && ra[:value] ≈ 0.5 && ra[:call] == "up"

        # without activities: resolution only, no value invented
        resp = SRV.readouts_handler(post("""{"network_id": "$nid",
            "queries": [{"stable_ids": ["R-HSA-2"]}]}"""))
        r = SRV.JSON3.read(String(resp.body))[:results][1]
        @test r[:basis] == "direct" && !haskey(r, :value)

        # absent optional tables are NAMED, not approximated over
        f2 = fields(); f2["nodes"] = nodes_part()
        parsed2 = SRV.JSON3.read(String(SRV.parse_handler(upload(f2)).body))
        resp = SRV.readouts_handler(post("""{"network_id": "$(parsed2[:network_id])",
            "queries": [{"stable_ids": ["R-HSA-9"]}]}"""))
        body2 = SRV.JSON3.read(String(resp.body))
        @test sort([String(t) for t in body2[:missing_tables]]) ==
              ["cofactors.csv", "entity_reaction_proxy_mapping.csv"]
        # the proxy step cannot run, so the query falls through honestly
        @test body2[:results][1][:basis] == "none"

        # a network parsed WITHOUT a node table refuses: flat reactome_id
        # matching would be a second, unvalidated protocol
        parsed3 = SRV.JSON3.read(String(SRV.parse_handler(upload(fields())).body))
        resp = SRV.readouts_handler(post("""{"network_id": "$(parsed3[:network_id])",
            "queries": [{"stable_ids": ["R-HSA-2"]}]}"""))
        @test resp.status == 400
        @test String(SRV.JSON3.read(String(resp.body))[:message]) == SRV.READOUTS_NO_TABLE_ERROR

        # an unknown network_id is the existing re-parse 400, not a fall-through
        resp = SRV.readouts_handler(post("""{"network_id": "pw:nope",
            "queries": [{"stable_ids": ["R-HSA-2"]}]}"""))
        @test resp.status == 400
        @test occursin("network_id", String(SRV.JSON3.read(String(resp.body))[:message]))

        # malformed requests name the shape they need
        for bad in ("""{"network_id": "$nid"}""",
                    """{"network_id": "$nid", "queries": []}""",
                    """{"network_id": "$nid", "queries": [{"stable_ids": "R-HSA-2"}]}""",
                    """{"network_id": "$nid", "queries": [{"stable_ids": []}]}""",
                    """{"network_id": "$nid", "queries": [{"stable_ids": ["R-HSA-2"], "mode": "entity"}]}""",
                    """{"network_id": "$nid", "activities": [1], "queries": [{"stable_ids": ["R-HSA-2"]}]}""")
            resp = SRV.readouts_handler(post(bad))
            @test resp.status == 400
            @test String(SRV.JSON3.read(String(resp.body))[:message]) == SRV.READOUTS_QUERY_ERROR
        end
    end

    @testset "a bad structural switch stops the server at start (review 2026-11-05)" begin
        for (k, bad) in ("DS_CONTEXT_BRIDGE_ORDER" => "late", "DS_SELF_ACTIVATOR_LEAVES" => "yes")
            old = get(ENV, k, nothing)
            ENV[k] = bad
            try
                @test_throws ArgumentError SRV.resolve_server_env()
            finally
                old === nothing ? delete!(ENV, k) : (ENV[k] = old)
            end
        end
        @test SRV.resolve_server_env() === nothing
    end

    @testset "/api/health answers 200 and carries the catalog" begin
        resp = SRV.health_handler(nothing)
        @test resp.status == 200
        body = SRV.JSON3.read(String(resp.body))
        @test haskey(body, :catalog)
        @test haskey(body[:catalog], :build_id)
        @test body[:status] == "ok"
    end
end
