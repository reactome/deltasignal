# The solver's answer depends on node names, not just on the graph

**Created**: 2026-09-18
**Status**: IN PROGRESS — defect measured and reproduced; fix implemented behind a flag, not yet decided
**Flag**: `DS_SCC_SWEEP` (`gauss_seidel` default = current behaviour, `jacobi` = order-free)

## How it was found

Adversarially reviewing LNG #89 (emit one-sided reactions). The PR's headline
gain was 82.95% → 83.36%, and `Transcriptional_Regulation_by_TP53` contributed
+84 of it. But TP53 gained **no** one-sided reactions — only 11 of 93 pathways
changed structurally and TP53 is not among them. Its network is canonically
identical across the two catalogs: mapping every UUID to its stable id gives the
same 5,418-edge multiset, same entities, types, signs, stoichiometries and
reaction ids.

A pathway that did not change produced 112 different predictions.

## The defect

**The solve is not a function of the graph.** Built a catalog that is `cat_os`
with every UUID renamed, verified an exact isomorphism (identical stable-id
canonical form for **92 of 92** pathways, 0 changed). Same solver, same config,
same observations.

Five labellings of that one network:

| labelling | accuracy | vs original | net | held-out net |
|---|---|---|---|---|
| original | 83.36% | — | — | — |
| perm | 83.42% | 14 cases | +14 | +0 |
| p101 | 83.41% | 53 cases | +13 | +5 |
| **p202** | **83.01%** | **112 cases** | **−84** | **+0** |
| p303 | 83.47% | 28 cases | +28 | +0 |

**Spread 0.47pp, stdev 0.187pp, from node names alone.** Not systematic: p202
is −84, the same magnitude as the +84 that LNG #89 appeared to gain on TP53.

Resolution is identical in every churned case (`n_gene_uuids`/`n_ko_uuids`
unchanged), so this is the solver, not the benchmark's readout matching.

## Mechanism

`steady_state.jl` builds node indices with `collect(keys(network.nodes))` and
`reaction_model.jl` iterates `target_groups`; both are Julia `Dict`s, so index
assignment and reaction order follow UUID hashing. That order reaches the sweep
inside a cyclic component, which is **Gauss-Seidel** — `x[t] = nv` writes back
mid-sweep, so each reaction reads values its neighbours already updated this
pass. The visit order therefore selects which fixed point a multi-root
component relaxes into. This is the all-zero-root problem of `specs/004`
showing up as a measurement artifact rather than as a modelling error.

**Confirmed at the pathway level: cycles are necessary.** Zero acyclic pathways
churn. All 10 churning pathways contain a cycle (median max-SCC 218 against 26
for the 62 quiet ones). Not sufficient — `DNA_Double-Strand_Break_Repair` has a
1127-node SCC and zero churn, which is what a unique fixed point looks like.

**No case-level discriminator, stated so nobody re-derives it.** Within TP53,
neither "readout is inside an SCC" (churned 3.6% vs stable 3.9%) nor "readout
is downstream of an SCC" (92.9% vs 92.1%) separates churned from stable cases.
That is the same wall recorded for saturation, path length and readout
in-degree.

## What this invalidates, and what survives

- **The headline accuracy figure carries a ±0.5pp label-noise band.** LNG #89's
  +0.41pp headline sits inside it and should not have been quoted as an effect.
- **The tuning-half +84 attributed to #89 is this defect**, not the fix.
- **Held-out numbers survive.** Worst case across four relabellings is +5 of
  18,808 (0.03pp); three of the four are exactly +0. The held-out +14 for
  #89 (16 fixed / 2 broke, p=0.0013) lands entirely in structurally-changed
  pathways and stands unaltered. The held-out protocol was already absorbing
  this, which is an argument for keeping every reported number on that split.

## The fix under test

`DS_SCC_SWEEP=jacobi` evaluates every reaction in a sweep against the state at
the sweep's start and commits the new values together, so the result is
invariant to the order of `rs`. Default stays `gauss_seidel`; a typo throws
rather than silently selecting a scheme.

This removes node naming as the thing that chooses between roots. It does
**not** explain why a component has more than one root — that is still
`specs/004`.

Assertions: `test/test_solver_determinism.jl` (80) and the `DS_SCC_SWEEP`
guard rails in `test/test_config_validation.jl` (175 → 187).

**Fixture caveat, recorded so a green run is not over-read.** The seven-node
cyclic fixture has a single fixed point. On it Gauss-Seidel deviates across
relabellings by only 2e-9 to 4e-9 — at the 1e-8 convergence tolerance, not at a
level that could flip a prediction — while Jacobi is exactly 0.0. The fixture
proves exact order-freedom; it does not reproduce the catalog's prediction
flips, which need a component that does not converge inside the iteration
budget.

## Open, measurements running

1. Is the churn a **convergence artifact**? `DS_MAX_ITERS=10000` on both
   catalogs. If a 20x budget removes it, the defect is "we stop before the
   root" rather than "there are several roots".
2. What does the order-free sweep **cost**? Jacobi on both catalogs, scored on
   the wide curator set with the held-out split.
3. Does Jacobi actually **converge** on the large SCCs, or does it trade
   order-dependence for non-convergence?

A decision on the default waits on 1–3.

## Amendment 1 (2026-10-03): the default sweep, pre-registered before the arm

**Why now.** A real-catalog relabelling happened by accident. Two builds have
networks identical by stable id, with curated rows in the same order; only
the minted uuids differ:
- canonical `20260928-1110_06ccb63`;
- `20261003-1058_1491276_f7ctrl`, the same generator content after the specs/044 split.

Under the default Gauss-Seidel sweep they differ by **158 TP53 predictions
(−145 cases)**, plus RUNX2 (45) and Fanconi (2); held-out is +2. TP53 is a
tuning pathway, so every tuning-half number is a coin flip across rebuilds.
specs/044 also showed the cost: a +96 TP53 "gain" in the F7 arm that vanished
under `jacobi`.

**Arm:** canonical split build `20260928-1110_06ccb63_split044`, whose
predictions are proven identical to canonical, scored with
`DS_SCC_SWEEP=jacobi`. It is compared with the existing `f7ctrl_jacobi` and
with both Gauss-Seidel scorings.

**Predictions:**
1. **Label invariance on the real catalog.** Jacobi on canonical against
   Jacobi on `f7ctrl` changes **0** predictions on both axes (Gauss-Seidel:
   205). Any non-zero count means order-dependence remains somewhere else
   (iteration over a Dict keyed by uuid outside the sweep, for example) and
   must be traced.
2. **Cost**, Jacobi against Gauss-Seidel on the same build:
   - held-out within ±15, the noise floor;
   - experimental within ±15.
   - On `f7ctrl` this was already seen once: curator all 85.22% vs 85.20%,
     experimental −2. That is not a pre-registration, so it is quoted here only
     as the reason for running.
3. **Convergence.** Report `Converged: N of M solves` for both sweeps. Jacobi
   may converge less often on the large components. A drop of more than 1% of
   solves is a finding to trace, not to accept.

**Decision rule:** make `jacobi` the default if (1) holds and (2) is within the
floor on both axes. Otherwise record it and keep Gauss-Seidel.

### Result (2026-10-03, solver ab0f44d / a952354)

Arms:
- `canon_jacobi`: build `20260928-1110_06ccb63_split044`, `DS_SCC_SWEEP=jacobi`;
- `f7ctrl_jacobi`;
- the two Gauss-Seidel scorings `split044` and `f7ctrl`.

1. **Label invariance holds exactly on the real catalog.** Jacobi on canonical
   against Jacobi on `f7ctrl` changes 0 predictions and 0 numeric values on
   both axes. Gauss-Seidel changes 205 + 14 predictions and 2,193 + 189
   values. So uuid-driven sweep order is the whole cause; nothing else in the
   solver depends on labels.
2. **Cost: fails.** Jacobi against Gauss-Seidel on the canonical build:
   - held-out **−94** (19 fixed / 113 broken, p = 1.8e-17);
   - tuning −42;
   - experimental −6 (0/6, p = 0.031).

   The held-out loss is DNA Double-Strand Break Repair −86 (KPNA2, MDC1,
   MRE11, NBN, KAT5) and DSB Response −14. Tuning is TP53 −42, nearly all MDM4
   (−49).
3. **Convergence: fails.** Converged solves:
   - curator: 1,669 → 1,609 of 1,725 (−3.5%);
   - experimental: 218 → 185 of 244 (−13.5%).

   Jacobi trades order dependence for non-convergence in the large
   components. That answers open question 3 above.

The near-zero difference seen first on `f7ctrl` (85.22% vs 85.20%) was not a
cost estimate: that build's Gauss-Seidel had landed in TP53's worse basin
(−145 against canonical).

**Decision: Gauss-Seidel stays the default.** Jacobi as it stands is not
adoptable. The label dependence remains open. The direct fix it points to is
to keep Gauss-Seidel and take its order from a **label-free key** (the
reaction's stable id and a structural signature, with ties broken by
structure) instead of from the uuid. That keeps the convergence properties and
removes the dependence on minted uuids. It is untried, and needs its own
pre-registration. Until then, every tuning-pathway number (TP53, RUNX2) is
quoted with the ±150-case rebuild band measured here, and the held-out split
remains the number of record.
