# specs/141: readout resolution and aggregation — the protocol of record, in
# the engine.
#
# Every citable figure reads a readout by a specific protocol that lived only
# in bench Python: resolve the curated entity (direct stable-id match, else
# proxy-reaction flux grouped by output entity, else set members), weigh each
# curated entity once (entity-mean, specs/105), and — for a total-protein
# question — keep only copies nothing makes from the protein itself
# (abundance_a3, specs/134). The Pathway Browser instead flat-means every
# node copy whose reactome_id matches, and shows nothing for expanded species
# or sets. This file is the port; the parity gate (specs/141 FR-003/FR-004,
# bench/analysis/141/) holds engine and bench equal on every benchmark
# readout of the frozen catalog. A drift here is a second protocol, which is
# the disease being treated.

const READOUT_PROXY_TABLE = "entity_reaction_proxy_mapping.csv"
const READOUT_COFACTOR_TABLE = "cofactors.csv"

# specs/134 amendment 2: identity edges — a set, pool or identity copy of the
# protein's copies IS the protein (bench CARRY_EDGES).
const ABUNDANCE_CARRY_EDGES = ("released_twin", "set_member", "specialization",
                               "variant_pool", "variant_split")
# Made-from edges under abundance_a3 (bench ABUNDANCE_RULES: MADE_FROM_EDGES
# minus recycle — amendment 3: a recycle edge returns a released copy to the
# starting copy, it does not make the starting copy from the protein).
const ABUNDANCE_MADE_FROM_EDGES = ("assembly", "dissociation", "released_twin",
                                   "set_member", "specialization",
                                   "variant_pool", "variant_split")

# The validated direction cutoffs (bench DS_DOWN_CUTOFF / DS_UP_CUTOFF
# defaults): fold < 0.85 is down, fold >= 1.15 is up.
const READOUT_DOWN_CUTOFF = 0.85
const READOUT_UP_CUTOFF = 1.15

"""The per-bundle tables the proxy step reads: entity stable id → proxy
reaction uuids (`entity_reaction_proxy_mapping.csv`, every role, row order and
duplicates kept, as the bench's `load_entity_reaction_proxies` keeps them) and
the cofactor stable ids (`cofactors.csv`). A missing table resolves as the
bench resolves it — no proxies / no cofactors — but is recorded so the API can
SAY it was missing instead of approximating silently."""
struct ReadoutIndex
    proxies::Dict{String, Vector{String}}
    cofactors::Set{String}
    has_proxy_table::Bool
    has_cofactor_table::Bool
end

"""Build a `ReadoutIndex` from the tables beside a bundle's logic network.
Both tables are optional (the direct and set-member steps need neither); a
present table with the wrong columns is an error, not an empty index."""
function build_readout_index(bundle_dir::AbstractString)::ReadoutIndex
    proxies = Dict{String, Vector{String}}()
    proxy_path = joinpath(bundle_dir, READOUT_PROXY_TABLE)
    has_proxy = isfile(proxy_path)
    if has_proxy
        df = CSV.read(proxy_path, DataFrame; types = String)
        for c in ("entity_stable_id", "proxy_uuid")
            c in names(df) || throw(ArgumentError("$proxy_path has no `$c` column; it is not a proxy table"))
        end
        for r in eachrow(df)
            (ismissing(r.entity_stable_id) || ismissing(r.proxy_uuid)) && continue
            push!(get!(Vector{String}, proxies, String(r.entity_stable_id)), String(r.proxy_uuid))
        end
    end
    cofactors = Set{String}()
    cof_path = joinpath(bundle_dir, READOUT_COFACTOR_TABLE)
    has_cof = isfile(cof_path)
    if has_cof
        df = CSV.read(cof_path, DataFrame; types = String)
        "stable_id" in names(df) || throw(ArgumentError("$cof_path has no `stable_id` column; it is not a cofactor table"))
        for r in eachrow(df)
            ismissing(r.stable_id) || push!(cofactors, String(r.stable_id))
        end
    end
    return ReadoutIndex(proxies, cofactors, has_proxy, has_cof)
end

# Reaction uuid → the non-cofactor entities its output edges make (bench
# `load_output_entities`: the entity map is EVERY nodes.csv row, set_pool
# included, which is exactly `PinIndex.diagram_entity`).
function _readout_output_entities(edges::Vector{LogicNetworkEdge}, pin_index::PinIndex,
                                  cofactors::Set{String})::Dict{String, Set{String}}
    out = Dict{String, Set{String}}()
    for e in edges
        e.edge_type == "output" || continue
        ent = get(pin_index.diagram_entity, e.child_uuid, "")
        (isempty(ent) || ent in cofactors) && continue
        push!(get!(Set{String}, out, e.parent_uuid), ent)
    end
    return out
end

"""Resolve a readout by the bench's ladder (`run_pathway`'s key-output block):
the first candidate stable id with a direct match (one group); else the first
candidate with proxy reactions, grouped by the sorted non-cofactor output
entities each proxy makes (a proxy making none groups under its own uuid);
else the client-supplied set members, direct-only, one group per member that
resolves, in order, deduplicated. Returns `(basis, groups)` with basis one of
`"direct"`, `"proxy"`, `"set_members"`, `"none"`."""
function resolve_readout(stids::Vector{String}, members::Vector{String},
                         pin_index::PinIndex, ridx::ReadoutIndex,
                         edges::Vector{LogicNetworkEdge})
    for sid in stids
        us = get(pin_index.stid_to_uuids, sid, String[])
        isempty(us) || return (basis = "direct", groups = [copy(us)])
    end
    for sid in stids
        us = get(ridx.proxies, sid, String[])
        isempty(us) && continue
        outent = _readout_output_entities(edges, pin_index, ridx.cofactors)
        key_order = String[]
        by_out = Dict{String, Vector{String}}()
        for u in us
            ents = get(outent, u, Set{String}())
            k = isempty(ents) ? u : join(sort(collect(ents)), "|")
            haskey(by_out, k) || push!(key_order, k)
            push!(get!(Vector{String}, by_out, k), u)
        end
        return (basis = "proxy", groups = [by_out[k] for k in key_order])
    end
    groups = Vector{Vector{String}}()
    seen = Set{String}()
    for m in members
        m in seen && continue
        push!(seen, m)
        g = get(pin_index.stid_to_uuids, m, String[])
        isempty(g) && continue   # the bench appends the empty group, then filters it
        push!(groups, copy(g))
    end
    isempty(groups) && return (basis = "none", groups = Vector{Vector{String}}())
    return (basis = "set_members", groups = groups)
end

"""specs/105 `entity_mean`: the mean over non-empty groups of each group's
mean activity, a missing uuid reading `baseline` (the bench's
`activities.get(u, BASELINE)`). `nothing` when every group is empty."""
function readout_entity_mean(groups::Vector{Vector{String}}, activities;
                             baseline::Float64 = 0.01)::Union{Nothing, Float64}
    means = Float64[]
    for g in groups
        isempty(g) && continue
        push!(means, sum(Float64(get(activities, u, baseline)) for u in g) / length(g))
    end
    isempty(means) && return nothing
    return sum(means) / length(means)
end

"""The validated direction call on a fold over rest (bench `classify`):
strictly below 0.85 is `"down"`, at or above 1.15 is `"up"`."""
readout_call(fold::Float64)::String =
    fold < READOUT_DOWN_CUTOFF ? "down" : (fold >= READOUT_UP_CUTOFF ? "up" : "normal")

"""The korkut total-readout resolution (`readout_nodes`): the union of direct
matches over ALL the readout's stable ids, order kept, deduplicated — one flat
group, no proxy or set-member steps."""
function resolve_abundance_occurrences(stids::Vector{String}, pin_index::PinIndex)::Vector{String}
    nodes = String[]
    seen = Set{String}()
    for sid in stids, u in get(pin_index.stid_to_uuids, sid, String[])
        u in seen && continue
        push!(seen, u)
        push!(nodes, u)
    end
    return nodes
end

"""specs/134 `abundance_copies` under `abundance_a3`: the members of `nodes`
that nothing makes from `nodes`. Close the set over carry edges (a set or
pool of the protein's copies is the protein, amendment 2); drop a node that a
made-from edge feeds from the carried set (recycle excluded, amendment 3) or
that a reaction produces while taking an input in the carried set. Catalyst,
regulator and depletion edges are not sources. Order kept."""
function abundance_copies(nodes::Vector{String}, edges::Vector{LogicNetworkEdge})::Vector{String}
    inputs = Dict{String, Set{String}}()
    producers = Dict{String, Set{String}}()
    ident = Dict{String, Set{String}}()
    for e in edges
        if e.edge_type == "input"
            push!(get!(Set{String}, inputs, e.child_uuid), e.parent_uuid)
        elseif e.edge_type == "output"
            push!(get!(Set{String}, producers, e.child_uuid), e.parent_uuid)
        end
        if e.edge_type in ABUNDANCE_CARRY_EDGES
            push!(get!(Set{String}, ident, e.parent_uuid), e.child_uuid)
        end
    end
    carry = Set{String}(nodes)
    todo = collect(carry)
    while !isempty(todo)
        v = pop!(todo)
        for w in get(ident, v, Set{String}())
            w in carry && continue
            push!(carry, w)
            push!(todo, w)
        end
    end
    made = Set{String}()
    for e in edges
        if e.edge_type in ABUNDANCE_MADE_FROM_EDGES && e.parent_uuid in carry
            push!(made, e.child_uuid)
        end
    end
    empty_set = Set{String}()
    return [u for u in nodes if !(u in made) &&
            !any(!isempty(intersect(get(inputs, r, empty_set), carry))
                 for r in get(producers, u, empty_set))]
end

"""The korkut statistic for an abundance readout: the mean over UNPINNED kept
copies, as a fold over rest. Returns `(status, value)` where status is
`"scored"`, `"trivial"` (every kept copy pinned — the readout leaves scoring)
or `"no_copy"` (nothing kept); value is `nothing` unless scored."""
function abundance_value(kept::Vector{String}, pinned::Set{String}, activities;
                         baseline::Float64 = 0.01)
    isempty(kept) && return (status = "no_copy", value = nothing)
    free = [u for u in kept if !(u in pinned)]
    isempty(free) && return (status = "trivial", value = nothing)
    v = sum(Float64(get(activities, u, baseline)) for u in free) / length(free)
    return (status = "scored", value = v)
end
