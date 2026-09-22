# Research: Loop as a conserved pool

All five unknowns resolved from the code in `src/solvers/steady_state.jl`
(`solve_scc_ordered!`, lines 569–820 on this branch) and
`src/core/reaction_model.jl`.

## R1. Evaluating a member reaction with in-component inputs at baseline

**Decision**: build `xpool = copy(x)`; for every member node `m` of component
`c` (`comp_id[m] == c`) that is not in `obs_set`, set `xpool[m] =
baseline_vec[m]`; evaluate `compute_reaction_output_vec(xpool, r; config)` for
each member reaction `r` and take `e_r = out / baseline_vec[r.target_idx]`.

**Rationale**: `compute_reaction_output_vec` reads activators
(`x[rxn.activator_indices[k]]`), inhibitors (`x[rxn.inhibitor_indices[k]]`,
line 1714) and depleters (`x[rxn.depletion_indices[k]]`, line 1735) from the
same vector, so one overlaid copy covers every input class. External inputs
keep their solved upstream values (the component is visited in topological
order, so they are final); internal ones read baseline, i.e. fold 1, which is
exactly "no internal gain". Pinned members keep their pinned value and
therefore act as entries.

**Alternatives**: the `supply` overlay used by `DS_SCC_BREAK_CATALYST` freezes
only catalyst back-edges at the component's *entry state* (which is baseline
for interior nodes anyway — specs/016 §5); it is a per-edge mechanism and
cannot express "all internal inputs". Rejected.

## R2. Detecting internal negative edges per component

**Decision**: for member reaction `r` (target `t`), an internal negative edge
exists iff any `i in r.inhibitor_indices` or `i in r.depletion_indices` has
`comp_id[i] == comp_id[t]`. Computed once per component from
`comp_rxns[c]` (built at line 690).

**Rationale**: these are the only two negative input classes the propagator
has (`pos_neg == "neg"` edges land in `inhibitors` or, for `edge_type ==
"depletion"`, in `depletions` — `reaction_model.jl:205–210`). Activators are
positive regardless of edge_type. `comp_neg_frac[c]` (line 697) already counts
negative in-component edges for the `DS_SCC_NEG_MODE` layer and can be reused
as the gate: `comp_neg_frac[c] > 0`.

## R3. Parity with multiple entries and multiple paths

**Decision**: build the component's internal signed node graph once:
`adj[src] = [(tgt, sign)]` for every member reaction `r` and every input index
`i` with `comp_id[i] == c`, `sign = +1` for activators, `-1` for inhibitors and
depleters. For each entry node `e` (target of a reaction with `e_r != 1`), BFS
from `e` assigning `s[e] = +1` and `s[v] = s[u] * sign(u→v)`; a node reached
with two different signs from the same entry marks the component
*inconsistent* → fall back to the fixed point. Member value:
`baseline[m] * prod_e e_fold^(s_e[m])`.

**Rationale**: in an SCC every node is reachable from every entry, so every
member gets a sign; inconsistency is exactly an odd negative cycle reachable
from that entry, which is genuine negative feedback and is the case the spec
leaves to iteration. Multiple entries combine as a product per the stated rule,
each with its own parity.

**Determinism**: fold multiplication is not associative in floating point, so
the product is taken over entry folds sorted by value (ties by exponent sign),
never by node index — bit-identical under relabelling (FR-008).

## R4. Surfacing pooled / iterated counts

**Decision**: `solve_scc_ordered!` returns a third value, a `NamedTuple`
`(pooled, iterated, fallback_negative, fallback_inconsistent, pooled_nodes)`;
the call site (line 925) stores them in `SolverResult.diagnostics` under
`"scc_method"`, `"scc_pooled"`, `"scc_iterated"`, `"scc_fallback_negative"`,
`"scc_fallback_inconsistent"`, `"scc_pooled_nodes"`. The API solve response
gains an additive `"scc"` object carrying the same keys (documented in
`docs/API.md`; nothing existing changes). The CLI provenance already records
`scc_method`.

**Rationale**: `diagnostics::Dict{String,Any}` exists for this; the API does
not currently echo diagnostics at all (`server.jl:762–769`), so an additive
object is the smallest contract change. A benchmark arm verifies the mode from
the server env *and* from a probe solve's `scc` counts, so an arm cannot run
the old solver by accident (constitution V).

## R5. Reusable fixtures

- `test/test_loop_elasticity.jl::posloop()` — U → A, A → r1 → B → r2 → A: the
  positive two-node loop with one external drive (US1 scenarios 1–2).
- `test/test_solver_determinism.jl::relabel`, `RELABELLINGS`,
  `solve_by_original` — the relabelling harness (FR-008, SC-002); its
  `cyclic_net()` has a single root and one component.
- New fixtures needed: four-entry AND loop (scenario 3), OR-entry loop
  (scenario 4), the X → M, M ⊣ T parity fixture (US3), an odd-negative-cycle
  fixture (US3 scenario 3), a three-component network for the count report
  (US2 scenario 2). All go in a new `test/test_loop_pool.jl` (the repo's rule:
  assertions in a file that has them).

## Constitution check

- I (processing is the solver's job): no generator change. ✓
- II (measure, then claim): both variants measured on the wide set, both axes,
  held-out, p; predictions committed first. ✓
- III (negative results are results): the decision rule in SC-004 makes a
  negative outcome a recorded result. ✓
- IV (default-OFF): `DS_SCC_METHOD` default `fixed_point` unchanged;
  differential check against the pre-feature tree in the tasks. ✓
- V (honest reporting): fallback counts in diagnostics and API. ✓
- VI (API boundary): additive field only. ✓

## Results

### Finding from the fixtures, before any arm

A cycle with exactly one negative edge is an odd negative cycle by definition,
and MDM2 ⊣ TP53 → MDM2 is exactly that. So `pool_parity` cannot pool the
component the spec named as its motivating case; it pools only components in
which every cycle has an even number of negatives (double inhibition =
positive feedback), and `pool` pools only components with no internal negative
at all. Every giant component in the catalog (TP53, WNT, DSB) carries
depletion or inhibitor edges, so under both spec'd variants they fall back to
the fixed point and the arm measures the old solver there. `pool_all` — the
literal reading of the proposal, internal negatives contribute nothing — is
the only variant that pools them, and it moves TP53 *with* MDM2, the wrong
direction for the AKT cases. The AKT → TP53 case therefore needs a
hierarchical treatment (pool the positive recycling sub-cycles, iterate the
negative feedback over the pooled units), which is out of scope here and
recorded as the follow-up.

### Pre-registration (committed before the arms ran)

Arms on `cat_os`, production defaults, baseline `ab_onesided.tsv` (23,908
cases), `DS_SCC_METHOD ∈ {pool, pool_parity, pool_all}`; each arm's probe solve
of TP53 must report `scc.pooled > 0` for `pool_all` and the counts for the
others, or the arm is void. Curator held-out / tuning with concentration
columns, McNemar p, false change, per-pathway net; experimental axis for each
variant; the 102 TP53 AKT1/AKT2-KO cases scored directly; relabel churn on
one cyclic pathway.

- **P1 (coverage of the rule).** `pool` and `pool_parity` pool few of the
  cyclic components in loop-heavy pathways (< 20% of TP53's, DSB's and WNT's
  components; the small positive recycling loops elsewhere are pooled); their
  held-out nets are within ±15 of zero and their false-change counts within
  ±10 of the baseline's 1,430 — mostly a null measurement of the old solver.
- **P2 (pool_all).** `pool_all` pools every cyclic component. It removes the
  coin-flips: on the relabelled catalog (`cat_perm`) it moves 0 predictions
  where `fixed_point` moved 14. False change in the loop-heavy pathways falls
  by ≥ 100 (rails and collapses cannot form). Held-out net vs baseline is
  **positive**, ≥ +20, distributed over ≥ 5 pathways and ≥ 10 genes.
- **P3 (the price).** Under `pool_all` the 102 AKT1/AKT2-KO TP53 cases stay
  wrong (≤ 10 correct; today 2), and TP53's tuning net is negative or flat:
  members downstream of an internal inhibitor move with the entry.
- **P4 (experimental).** `pool_all` is not negative on the experimental axis
  (conditioned net ≥ −5).
- **P5 (decision rule).** Adopt a variant only if P2 and P4 both hold; if P2's
  sign is set by one pathway (> 50% of the net), report it as such and do not
  adopt. If `pool_all` loses held-out, the pool rule is refuted as a blanket
  treatment and the hierarchical variant is the only path left.

### Review, 2026-09-19 (two independent reviewers, code and data only), and one trace

**Pre-registration error, corrected before scoring.** P3 and SC-003 said the
102 AKT1/AKT2-KO TP53 cases are "currently 2 correct". In `ab_onesided.tsv`
they are **100 of 102** (94 UP→UP, 6 DOWN→DOWN, 2 UP→NORMAL); the 2 was the
shared catalog's figure (specs/016). P3 is re-stated: **no variant may lose
more than 10 of the 100**; its original form would have scored a −90 TP53
regression as "held". The baseline solve for that case is itself
non-converged (500 iterations), so the 100 are a rail landing on the right
side — the knife-edge premise stands, the starting point was misdescribed.

**pool ≡ pool_parity on this catalog.** Census of `cat_os` (reviewer's
`scc_census.py`, Tarjan over the edge lists): 258 cyclic components in 77
pathways; 60 carry an internal negative edge; **0 of the 60 are
sign-balanced**. Both variants pool the same 198 small positive loops (5,578
of 12,749 cyclic nodes) and fall back on every giant; their headlines are
identical (20,048/24,100). Three arms are two.

**P1 failed for `pool`, and not as a null result.** Held-out **−53** (15
fixed / 68 broke, p < 0.0001), tuning −16, false change on held-out keys
984 → 1,046. The broken cases sit at exactly **0.0 (27), 100.0 (24) or 80.0
(10)** where the baseline read exactly 1.0: the fixed point held these small
positive loops at baseline; the pool is what rails them. So the spec's "Why"
premise — any leak rails a positive cycle — is contradicted by the data for
the very loops these variants pool.

**Traced (benchmark-faithful gene resolution): HRAS KO → R-HSA-354126 in
Signaling by MET.** HRAS resolves to 4 nodes; one sits in a 19-node component
with no internal negative edge — the RAS GTP/GDP recycling cycle
(R-HSA-8851827 / 8851877 / 8851899) that also carries the other isoforms'
routes (R-HSA-8875568, 5674631, 354074, 8875591). Fixed point: the HRAS
reactions read 0, **those other members stay at 1.000**, the readout stays
NORMAL — the curators' call (isoform redundancy). Pool: the pinned HRAS node is
an entry with fold 0, **0 absorbs, all 19 members read 0**, the readout reads
DOWN. Product-across-entries is AND semantics applied to alternatives: a
component that bundles redundant routes is treated as if every route were
co-required. This is the mechanism behind MET −20, SCF-KIT −12, DAP12 −12,
and it is a property of the rule, not a bug.

**Solver review (reviewer A), verified and acted on:**

| # | finding | action |
|---|---|---|
| 1 | A pinned member does not sever the pool; entries on both sides of a pin multiply into one fold (the fixed point treats a pin as a boundary). `fourentry` with C pinned at baseline + E1 = 3x: fixed point D = 1, pool D = 3. | Recorded as a consequence of the rule (the granularity loss that was accepted); documented in `pool_component!`; not changed. |
| 2 | Exact `f != 1.0` entry test: hill_sat's smooth cap returns bl·(1 + 3e-15) at pure baseline, so ~97% of member reactions were "entries" (TP53 808 of 837). | **Fixed**: relative tolerance 1e-9; test asserts a resting loop pools with a real product of 1. `pool`/`pool_parity` arms ran with the exact test (factors of 1 + 3e-15, numerically inert); later arms carry the fix. |
| 3 | `converged` is false for every solve that pools a real entry (a pooled state is not a fixed point of F). | Documented; asserted in the test. The benchmark's "both arms converged" statistic is void under pooling. |
| 4 | The rule overrides well-posed contracting loops (OR-producer loop: unique stable fixed point 2.0, pool 5/3), double-counts one signal entering by two member reactions (4 vs the fixed point's rail), and lets two entries of opposite parity cancel exactly. | Consequences of the product rule, now stated as such in the docstring rather than pinned as truths. |
| 5–6 | Heterogeneous baselines create phantom entries (unreachable from the parser); `DS_SCC_BREAK_CATALYST` is a no-op inside a pooled component. | Documented. |
| 7 | Test gaps: dead code in FR-008, `ref` re-solve, pinned-member test weak, `converged` never asserted, file not in CI. | Fixed the first four; CI list and CLAUDE.md table updated in this commit. |
| 8 | `scc.method` reads `"flat"` under `DS_SCC_SOLVE=0`, not a `DS_SCC_METHOD` value; `pooled_nodes` counts reaction nodes. | Documented in docs/API.md. |

Default bit-identity vs the base branch was verified by the reviewer on the
TP53 catalog network with 25 pins (0 activity diffs) — the differential check
T015 is thereby done by an independent hand.

**Methodology review (reviewer B), verified:** the `--max-edges` cap differs
between baseline (20,000) and arms (40,000), adding CD28 (192 cases, 97.9%) to
the arms' headline — paired comparisons are on the 23,908 shared keys and are
unaffected; a `fixed_point` control on the current tree and cap is queued so
the −53 is attributed to the rule and not to the tree. The spec cited specs
013–016, which are on a sibling branch (#55), not an ancestor of this one. P5
pre-committed to the hierarchical follow-up whichever way the result fell;
that sentence is withdrawn — the follow-up is argued from the trace above, not
from the decision rule.

**Where this leaves the proposal before the remaining arms land.** The
product-of-entries pool is refuted as a blanket rule by a traced case: it
destroys the redundancy the fixed point preserves. What survives of the idea
is the part the trace does not touch — that a *pure recycling* cycle (one
species cycling through states, one entry) should not amplify — and the
remaining arms (`pool_all`, experimental, relabel churn, control) are run for
the record, not for a decision.

### Scoring, all arms landed (2026-09-19 23:15)

Control: `fixed_point` on the current tree with `--max-edges 40000`
(`ab_fixed_point_ctrl.tsv`) is **bit-identical** to `ab_onesided.tsv` on all
23,908 shared cases (0 predictions, 0 `pred_ui` differ), so the tree carries
no effect and CD28's 192 extra cases only move headlines. The `pool_all`
re-run on the committed tree (`fadc0b2`, clean) is bit-identical to the
original `pool_all` arm, so the mid-run entry-test edit was value-neutral and
the churn number below is attributable.

| variant | held-out (fixed/broke, p) | tuning | false change (23,908 keys) | experimental (849, conditioned) | AKT-KO TP53 (of 100) | relabel churn |
|---|---|---|---|---|---|---|
| `pool` | **−53** (15/68, <1e-4) | −16 | 1,437 → 1,518 | 0 (612 → 612) | 100 | — |
| `pool_parity` | identical to `pool` (0 values differ) | | | 0 | 100 | — |
| `pool_all` | **−326** (168/494, <1e-4), 20 pw, 89 genes | −399 | 1,437 → 1,771 | **−17** (24/41, p 0.046) | **2** | **0 of 24,100** (fixed point: 14) |

- **P1** — fails: not a null measurement; −53 by railing small positive loops the fixed point held at baseline.
- **P2** — fails: held-out −326, false change +334 (predicted −100 in loop-heavy pathways; observed +214 there). The one clause that held: 0 relabel churn.
- **P3** (corrected) — fails: `pool_all` loses 98 of the 100 AKT-KO cases. (The original wording, built on the wrong baseline, would have scored this as held.)
- **P4** — fails for `pool_all` (−17, p 0.046); `pool` / `pool_parity` neutral.
- **P5** — nothing adopted; the product-of-entries pool is refuted as a blanket treatment on both axes.

Mechanism, from the trace: the components the rule pools bundle alternative routes (RAS isoforms in one GTPase cycle) and, in the giants, genuine negative feedback; a product over entries is AND semantics applied to alternatives, and 0 absorbs. The follow-up is not a variant of the pool: specs/018 shows the giant components are welded by our own derived edges (assembly, depletion), which is a different problem with a different rule.

Process note for the record: `src/` was edited twice while arms mounting the live tree were running (the entry-test fix during the pool arms; the specs/018 implementation during the control and re-run). Both were shown harmless after the fact — the servers had loaded their code at start (`code: fadc0b2 dirty_src=0` in the arm log), and the affected dumps are bit-identical to their comparators — but the rule (arms mount a pinned worktree; no `src/` edits until the queue is empty) was written after the first and broken by the second. It stands.
