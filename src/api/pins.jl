# specs/140: pin resolution — the protocol of record, in the engine.
#
# Every benchmark number this project cites was measured by perturbing ONLY
# root inputs that are, or contain, the gene (DS_PIN_SCOPE=root_cycle with
# DS_PIN_ROOT_IGNORE=depletion; specs/023, 024, 058, 061), with every
# downstream copy COMPUTED, never set. That resolution lived only in
# bench/benchmark_vs_mpbiopath.py (`pin_occurrences`), which no API client
# runs: the Pathway Browser pinned every copy of the clicked entity, the
# `all`-family protocol specs/023 retired as not comparable. This file is the
# port. The parity gate (specs/140 FR-003, bench/analysis/140/) holds the two
# implementations equal on the frozen catalog, case for case; a drift here is
# a second protocol, which is the disease being treated.

# In-edges that the root_cycle fallback treats as derived rather than as a
# producer (bench DERIVED_IN_EDGES).
const PIN_DERIVED_IN_EDGES = ("dissociation", "depletion")

"""Index of pin targets, from the generator's `nodes.csv`: a node is indexed
by its `diagram_entity_id` and every stable id in its `member_leaves`
(is-OR-contains); `set_pool` nodes are never pin targets (specs/033: the pool
is the set as one catalyst, not a pinnable occurrence). `diagram_entity`
backs the fallback's own-entity preference (specs/024)."""
struct PinIndex
    stid_to_uuids::Dict{String, Vector{String}}
    diagram_entity::Dict{String, String}
end

"""Build a `PinIndex` from the `nodes.csv` beside a bundle's logic network,
or `nothing` when the bundle has none (a pre-specs/046 build, or the sample
network): resolution is then refused, not approximated."""
function build_pin_index(nodes_csv_path::AbstractString)::Union{Nothing, PinIndex}
    isfile(nodes_csv_path) || return nothing
    df = CSV.read(nodes_csv_path, DataFrame; types = String)
    for c in ("uuid", "node_kind", "member_leaves")
        c in names(df) || throw(ArgumentError("$nodes_csv_path has no `$c` column; it is not a generator node table"))
    end
    has_diagram = "diagram_entity_id" in names(df)
    stid_to_uuids = Dict{String, Vector{String}}()
    diagram_entity = Dict{String, String}()
    seen = Set{Tuple{String, String}}()
    for r in eachrow(df)
        ismissing(r.uuid) && throw(ArgumentError("$nodes_csv_path: a row has no uuid"))
        uuid = String(r.uuid)
        kind = ismissing(r.node_kind) ? "" : String(r.node_kind)
        de = has_diagram && !ismissing(r.diagram_entity_id) ? strip(String(r.diagram_entity_id)) : ""
        isempty(de) || (diagram_entity[uuid] = String(de))
        kind == "set_pool" && continue
        keys = Set{String}()
        isempty(de) || push!(keys, String(de))
        if !ismissing(r.member_leaves)
            for m in split(String(r.member_leaves), '|')
                m = strip(m)
                isempty(m) || push!(keys, String(m))
            end
        end
        for k in keys
            (k, uuid) in seen && continue
            push!(seen, (k, uuid))
            push!(get!(Vector{String}, stid_to_uuids, k), uuid)
        end
    end
    return PinIndex(stid_to_uuids, diagram_entity)
end

"""The in-degrees that decide which nodes are roots, the bench's
`root_indegree` ported edge for edge: `recycle` edges count for nothing
(specs/081), a `context_bridge` edge's source has an out-edge but its target
gains no in-degree (specs/056: roots are read without context), and an edge of
type `root_ignore` (default `depletion`, specs/058) is ignored — except that a
node fed ONLY by ignored edges and with no out-edge is not a root (week review
2026-10-15: a depletion-fed dead end pins nothing that moves)."""
function pin_root_indegree(edges::Vector{LogicNetworkEdge};
                           root_ignore::AbstractString = "depletion")::Dict{String, Int}
    indeg = Dict{String, Int}()
    ignored_in = Set{String}()
    has_out = Set{String}()
    for e in edges
        e.edge_type == RECYCLE_EDGE && continue
        push!(has_out, e.parent_uuid)
        e.edge_type == CONTEXT_BRIDGE_EDGE && continue
        if !isempty(root_ignore) && e.edge_type == root_ignore
            push!(ignored_in, e.child_uuid)
            continue
        end
        indeg[e.child_uuid] = get(indeg, e.child_uuid, 0) + 1
    end
    for u in setdiff(ignored_in, has_out)
        get(indeg, u, 0) == 0 && (indeg[u] = 1)
    end
    return indeg
end

"""`DS_PIN_ROOT_IGNORE` as the bench reads it: `depletion` by default
(specs/061 adoption); `none` restores the specs/024 rule."""
function pin_root_ignore()::String
    v = get(ENV, "DS_PIN_ROOT_IGNORE", "depletion")
    return v == "none" ? "" : v
end

# The occurrences of a query: every node indexed by any of its stable ids,
# file order kept, duplicates dropped (bench: gene_to_uuids built per stid,
# then deduplicated by root_occurrences / recycled_root_occurrences).
function _pin_occurrence_uuids(stids::Vector{String}, index::PinIndex)::Vector{String}
    uuids = String[]
    seen = Set{String}()
    for sid in stids, u in get(index.stid_to_uuids, sid, String[])
        u in seen && continue
        push!(seen, u)
        push!(uuids, u)
    end
    return uuids
end

# specs/024 fallback, the bench's `recycled_root_occurrences`: for a query
# with occurrences but no root, the nodes with at least one out-edge whose
# every producer is downstream of them or reaches them only by a derived
# dissociation/depletion edge; the query's own entities preferred, else the
# qualifying nodes no other qualifying node reaches (or, if they all
# regenerate each other, all of them). `inc`/`fwd` are built without recycle
# and context_bridge edges, as the bench builds them.
function _recycled_root_occurrences(uuids::Vector{String}, own_stids::Set{String},
                                    index::PinIndex, edges::Vector{LogicNetworkEdge})::Vector{String}
    inc = Dict{String, Vector{Tuple{String, String}}}()
    fwd = Dict{String, Vector{String}}()
    for e in edges
        (e.edge_type == CONTEXT_BRIDGE_EDGE || e.edge_type == RECYCLE_EDGE) && continue
        push!(get!(Vector{Tuple{String, String}}, inc, e.child_uuid), (e.parent_uuid, e.edge_type))
        push!(get!(Vector{String}, fwd, e.parent_uuid), e.child_uuid)
    end

    cache = Dict{String, Set{String}}()
    function reach(u::String)::Set{String}
        get!(cache, u) do
            seen = Set{String}()
            stack = copy(get(fwd, u, String[]))
            while !isempty(stack)
                n = pop!(stack)
                n in seen && continue
                push!(seen, n)
                append!(stack, get(fwd, n, String[]))
            end
            seen
        end
    end

    function fed_only_by_itself(u::String)::Bool
        # A node with no out-edges (a dissociation sink) cannot carry a
        # perturbation anywhere, so it is never a root (review of PR #73).
        isempty(get(fwd, u, String[])) && return false
        return all(et in PIN_DERIVED_IN_EDGES || s in reach(u)
                   for (s, et) in get(inc, u, Tuple{String, String}[]))
    end

    own = [u for u in uuids if get(index.diagram_entity, u, "") in own_stids && fed_only_by_itself(u)]
    isempty(own) || return own
    cand = [u for u in uuids if fed_only_by_itself(u)]
    entry = [u for u in cand if !any(u in reach(v) for v in cand if v != u)]
    return isempty(entry) ? cand : entry
end

"""Resolve a perturbation query (a list of stable ids: the clicked entity, or
every form of a gene) to the nodes the protocol of record pins. Returns
`(pins, basis, occurrences)`: the pinned uuids; `"root"`, `"root_cycle"` (the
specs/024 fallback fired) or `"none"` (no pinnable form — the honest answer a
client must surface instead of pinning mid-pathway copies); and how many
occurrences the query had before the root filter."""
function resolve_pins(stids::Vector{String}, index::PinIndex, edges::Vector{LogicNetworkEdge};
                      root_ignore::AbstractString = pin_root_ignore())
    uuids = _pin_occurrence_uuids(stids, index)
    isempty(uuids) && return (pins = String[], basis = "none", occurrences = 0)
    indeg = pin_root_indegree(edges; root_ignore = root_ignore)
    roots = [u for u in uuids if get(indeg, u, 0) == 0]
    isempty(roots) || return (pins = roots, basis = "root", occurrences = length(uuids))
    pins = _recycled_root_occurrences(uuids, Set(stids), index, edges)
    return (pins = pins, basis = isempty(pins) ? "none" : "root_cycle",
            occurrences = length(uuids))
end
