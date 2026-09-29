"""
Self-fed inputs (specs/040, rule A) and leaf-level self-contained inhibitors
(specs/040, rule B). One principle: a reaction must not read its own signal
twice.

**Rule A.** In a cyclic component that the solver iterates, an input that is
fed by nothing except what it feeds carries no information its consumers do
not already have. Multiplying it in multiplies the product into itself: the
RAF/MAP kinase trace has three p-MEK dimers in the "RAF activating kinases"
set pool and 22 scaffold leaves at the next step, all produced only by the
dissociation downstream, so the loop gain is fold^24 and the component finds
the all-zero root under every perturbation.

Per component, when its iteration starts (upstream is solved):
- *signal-carrying entries*: nodes pinned off baseline, and unpinned targets
  with an activator outside the component that is off baseline. A pin AT
  baseline (an inert cofactor, a held drug) is a constant: not an entry, and
  the search below does not pass through it. Pool-owned nodes (specs/039
  states, intermediates, step copies) are never candidates.
- Out(u): the in-component targets u feeds. A set-pool node is transparent,
  so its members' Out is the pool's Out.
- u is **self-fed** iff every activator path from an entry to u inside the
  component passes through Out(u): fed only by what it feeds.
- Every edge from a self-fed u into Out(u) reads the component-entry value
  (the specs/018 `activator_break` / `supply` plumbing). No entry reaches u
  except through its consumers, so that value is exactly baseline. A pin
  between consumer and u (MEK knocked out) makes u entry-fed, so it reads live.

The set is a dominator set: a property of the graph and the entry state,
not of node labels or visit order. It is computed once per component per
solve; it does not change between sweeps.

  DS_SELF_FED_MODE=off    (default) previous behaviour.
  DS_SELF_FED_MODE=entry  the rule above (specs/040 rule A, measured and
                          refuted: it held TP53's and PIP3's signal routes).
  DS_SELF_FED_MODE=multi  rule A2 (specs/040 amendment 1). An edge from u into
                          a step is held only if BOTH hold:
                          (i) u is dominated, not merely unreached: some entry
                              reaches u by activator paths with Out(u)
                              passable, and none with Out(u) blocked. Under
                              `entry` "every path passes through Out(u)" was
                              vacuously true for a loop whose signal enters
                              through a negative edge (TP53, PIP3), and every
                              node of it was held.
                          (ii) the recycled fold is multiplied into the step:
                              the step's target is a set-pool node, or the
                              step reads >= 2 dominated AND inputs (RAF's 21
                              scaffold leaves; its 3 p-MEK dimers in a pool).
                              A width-1 recycling (PIP2 <-> PIP3) stays live.
                          The solve reports `self_fed_rule`, `self_fed_nodes`
                          (nodes with a held edge) and `self_fed_edges_held`.

**Rule B.** specs/022 flags an inhibitor that CONTAINS its reaction's input
(whole-node stId containment) and is computed from it. The PEBP1 inhibitor of
"MAP2Ks and MAPKs bind to the activated RAF complex" is built from the
RAF1-only sibling variant of the step's input, so neither test holds, yet it
tracks the input exactly and cancels it. Under DS_SELF_INHIBITOR_LEAVES=1 a
second clause is added: the inhibitor and the input share a non-cofactor LEAF
(a contained stId that contains nothing else), and a carrier of that leaf
reaches both — they are built from the same species. Same formula, same
weight, still weaken-only; specs/022's pairs are a subset. The solve reports
`self_inhibitor_leaves` and `self_inhibitor_leaf_pairs` (slots flagged by the
new clause alone).

Both rules are exact at baseline (every held edge reads fold 1 there) and
byte-identical off: they add flags and pairs, never nodes or edges.
"""

const DS_VALID_SELF_FED_MODES = Set(["off", "entry", "multi"])

function self_fed_mode()::String
    value = get(ENV, "DS_SELF_FED_MODE", "off")
    value in DS_VALID_SELF_FED_MODES || throw(ArgumentError(
        "DS_SELF_FED_MODE must be one of $(join(sort(collect(DS_VALID_SELF_FED_MODES)), ", ")); got $(repr(value))"))
    return value
end

const DS_VALID_SELF_INHIBITOR_LEAVES = Set(["0", "1"])

function self_inhibitor_leaves()::Bool
    value = get(ENV, "DS_SELF_INHIBITOR_LEAVES", "1")   # default on since 2026-09-28 (specs/040 result)
    value in DS_VALID_SELF_INHIBITOR_LEAVES || throw(ArgumentError(
        "DS_SELF_INHIBITOR_LEAVES must be 0 or 1; got $(repr(value))"))
    return value == "1"
end

# "Off baseline" for entry detection. Upstream values at rest are products of
# exact 1.0 folds, so anything above float noise is a signal.
const SELF_FED_TOL = 1e-12

"""
The self-fed nodes of one component (specs/040 rule A).

- `nodes_c`: the component's node indices;
- `act_out[u]`: the IN-COMPONENT targets of u's activator edges;
- `is_pool`: set-pool nodes (transparent: their targets stand in for them);
- `skip`: nodes that are never candidates (set pools, pool-owned nodes);
- `entry`: the signal-carrying entries;
- `blocked`: nodes pinned at baseline: constants, never candidates, and the
  search may not pass through them.

u is self-fed iff no entry reaches u along in-component activator edges without
passing through Out(u). With `dominated_only` (rule A2, condition i), u must
ALSO be reached by some entry once Out(u) is passable: a node no entry reaches
at all is unreached, not dominated, and stays live. Deterministic: the result
is a set, and the search visits nodes in sorted order.
"""
function self_fed_nodes(nodes_c::Vector{Int}, act_out::Dict{Int, Vector{Int}},
                        is_pool::Set{Int}, skip::Set{Int}, entry::Set{Int},
                        blocked::Set{Int}; dominated_only::Bool = false)::Set{Int}
    out = Set{Int}()
    isempty(entry) && return out
    seen = Set{Int}()
    outs = Set{Int}()
    for u in sort(nodes_c)
        (u in entry || u in skip || u in blocked) && continue   # a pinned node is a constant
        empty!(outs)
        for t in get(act_out, u, Int[])
            if t in is_pool
                for t2 in get(act_out, t, Int[]); push!(outs, t2); end
            else
                push!(outs, t)
            end
        end
        isempty(outs) && continue
        empty!(seen); union!(seen, entry)
        stack = sort(collect(entry))
        reached = false
        while !isempty(stack)
            v = pop!(stack)
            if v == u; reached = true; break; end
            (v in outs || v in blocked) && continue
            for w in get(act_out, v, Int[])
                w in seen && continue
                push!(seen, w); push!(stack, w)
            end
        end
        reached && continue
        if dominated_only
            # Condition (i): with Out(u) passable, does any entry reach u at all?
            empty!(seen); union!(seen, entry)
            stack = sort(collect(entry))
            any_reach = false
            while !isempty(stack)
                v = pop!(stack)
                if v == u; any_reach = true; break; end
                v in blocked && continue
                for w in get(act_out, v, Int[])
                    w in seen && continue
                    push!(seen, w); push!(stack, w)
                end
            end
            any_reach || continue   # unreached: not dominated, stays live
        end
        push!(out, u)
    end
    return out
end

"""
Whether an indexed reaction is a set POOL node (every input a `set_member`
edge, specs/033): the structural test `compute_reaction_output_vec` applies.
"""
is_set_pool_reaction(r::IndexedReaction)::Bool =
    !isempty(r.activator_group) && all(==(-1), r.activator_group) &&
    isempty(r.inhibitor_indices) && isempty(r.depletion_indices)
