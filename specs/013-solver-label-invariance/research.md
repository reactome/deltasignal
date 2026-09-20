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
