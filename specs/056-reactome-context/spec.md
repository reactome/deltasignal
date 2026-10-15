# specs/056: Reactome context — each pathway solved with the pathways that feed it

**Status:** **claim not supported** (2026-10-14; research.md R1): both arms
lose on both axes. Registered before either arm ran and before any case was
read for it. Requested by Adam, 2026-10-14. This is a **second protocol**,
reported next to the publication protocol and never in place of it.

## Why

Adam (2026-10-14): MP-BioPath's protocol solves each pathway alone, so some
"wrong" predictions may be biology that Reactome curates in another pathway.
If solving a pathway with the pathways that feed it fixes such cases,
including through pathways outside the test sets, that shows the tool is not
wrong, only bounded by the pathway boundary. This must not become the
headline, because the publication did not use it.

History, on a much older solver (memory, 2026-07-16):
- 82% of failing experimental cases had a curated route somewhere in the 92
  catalog pathways;
- merging all 92 by node identity lost (experimental macro-F1 0.66 → 0.48),
  with hubs (active ERK, RAS, PIP3) flooded and readouts that should not change
  predicted to change (NORMAL-F1 0.39 → 0.18); a hub-excluded stitch at 5
  pathways also lost (0.52);
- the old joint networks (2026-05) helped the cell-cycle pathways and broke
  PIP3; they were chosen per pathway after seeing results, so are not cited.

What has changed since: welds removed (specs/018, 030), variant nodes (046),
cycle balance (039), carrier loops (049), cofactors, pathogens and drugs held
(007, 048, 055). The July result is a reason for caution, not a measurement
of this solver.

## The rule (`bench/context_build.py`, default off: a separate catalog build)

From a built catalog, each pathway P becomes a bundle of P and its
neighbours:
1. **Identity.** An entity node is identified by its full `nodes.csv` row
   minus the uuid. Reactions, set pools, cofactors, drugs and pathogens are
   never bridged.
2. **Hubs.** An entity whose diagram entity is in more than **5** of the 92
   pathways is never bridged (361 of 11,688 entities). Fixed; not swept.
3. **Reduction.** A neighbour Q first loses every reaction P already has (its
   reaction nodes, their edges, and edges attributed to it), so it brings
   only what P lacks. Without this, half of the first build's bridges closed
   cycles, mostly a parent pathway re-linked to its own sub-pathway's copies
   of the same reactions.
4. **Bridges.** `pos,or,context_bridge` from a copy produced in one pathway
   (an in-edge there, after reduction) to a copy of the same identity that
   is a root in the other. Several producers into one root read as their
   mean. No uuid is merged; nothing inside a pathway is bridged; only P and
   one neighbour at a time (a star around P).
5. **Neighbours** are the pathways with at least one bridge INTO P.
6. **Pruning.** A neighbour keeps only the ancestors, inside it, of its
   bridges into P, closed over its pools and carrier loops. Checked on three
   pathways against the unpruned bundle, pinning 12 own roots each to 0 and
   80 (`bench/analysis/056/prune_exact.jl`): two identical, IL-2 within 5e-7
   (internal scale; 5e-5 of baseline), converged flags equal.

Two variants:
- **`ctxin` (feed-in).** No bridge leaves P. Neighbours only feed P's roots,
  so the context cannot close a cycle through P.
- **`ctxfb` (feedback).** Bridges also run from P's products to neighbour
  roots, so a signal can leave P and come back (the shape of RAF's ERK
  negative feedback, though that pathway, R-HSA-5674499, is not in the
  catalog).

## Benchmark protocol under context

- **Pins.** The gene's root forms, read with bridges ignored, in P and in
  every neighbour part of the bundle (a knockout is absent everywhere). A
  root of a reduced neighbour includes a node whose only producer was one of
  P's reactions; it is also fed by P's copy through a bridge (`ctxfb`).
- **Readouts** are P's own nodes only (`context_own.csv`), so what is measured
  is unchanged.
- Solver at code defaults. The bench log's `Protocol:` line names the context
  (`context=<base>+hub5[,feed-in]`).

## Measured scope (structural, canonical build `20261012-1401_01206ca`)

| | `ctxin` (`…ctx5in`) | `ctxfb` (`…ctx5`) |
|---|---|---|
| builder commit | aa3a691 | aa3a691 |
| pathways with a neighbour | 79 of 92 | 79 |
| mean neighbours | 4.9 | 4.9 |
| bridges | 1,040 | 1,417 |
| bridges closing a cycle | 0 | 383, in 26 pathways (DNA repair family, TP53) |
| nodes, all bundles | 213,171 → 479,929 | same |

Experimental pathways: RAF 16 neighbours, Cell Cycle Checkpoints 12, S phase
13, Mitotic G1 6, TP53 16, PIP3 19.

## Arms

All three by `scripts/run_arm.sh`, curator and experimental axes, solver at
the registration commit.

| arm | catalog | compared against |
|---|---|---|
| `ctx0` | canonical `20261012-1401_01206ca` | — (control; must reproduce `results/64d23cb` on every prediction) |
| `ctxin` | `--catalog 20261012-1401_01206ca.ctx5in` | `ctx0` |
| `ctxfb` | `--catalog 20261012-1401_01206ca.ctx5` | `ctx0` |

## Predictions

- **C-P1 (structural).** The 13 pathways with no neighbour predict exactly as
  `ctx0`, in both arms.
- **C-P2.** Cases invalid in `ctx0` whose gene has a root only in a neighbour
  become valid. They are counted and scored separately, never in the paired
  comparison.
- **C-P3.** The gain, if any, is on the experimental axis, in the cell-cycle
  pathways (Mitotic G1, S phase, checkpoints), whose 2026-05 failures were
  cross-pathway (CCND1:CDK4 phosphorylates RB1 in G1 phase, not in Mitotic
  G1). No direction is predicted for RAF, TP53 or PIP3.
- **C-P4 (named risk).** Readouts that should not change are predicted to
  change: the July mechanism. Its measure is the change in correct NORMAL
  predictions on both axes. For `ctxfb`, the DNA repair family and TP53 carry
  the cycle-closing bridges and are where a loss would come from.
- **C-P5.** Curators predicted each pathway from inside it, so curator
  held-out may not gain and may lose.

## Decision rule

Nothing here changes the default model or the headline protocol. Both arms
are reported whatever they show, next to the within-pathway numbers.

An arm **supports the claim** ("solving with the pathways that feed it
explains cases the pathway alone gets wrong") against `ctx0`, on cases valid
in both, when:
- experimental net > 0 at McNemar p < 0.05, read with
  `--set-aside RAF_MAP_kinase_cascade` (specs/054 D1; both figures reported);
- curator held-out net loss, if any, not significant (p ≥ 0.05);
- no pathway loses more than 10 cases on either axis;
- C-P1 holds and `ctx0` reproduces `results/64d23cb`;
- every experimental fix is traced to a named route that crosses a bridge
  (neighbour pathway and reaction named). A fix whose route crosses no bridge
  is not counted for the claim.

The fixes are also counted by whether the bridge's neighbour is a curator or
experimental test pathway or not.

If an arm supports the claim, a separate registration scores it **once** on
korkut2015 (TEST, specs/051, values unread) against the within-pathway
protocol. Until then the claim is a SELECT-set result, not a test result.
No hub cutoff, neighbour choice or per-pathway selection is tuned against
these results; a changed rule is a new registration.
