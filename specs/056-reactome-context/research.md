# specs/056 research

Registered 2026-10-14 (spec.md, f9c82b3). All three arms at solver f9c82b3,
`scripts/run_arm.sh`:
- `ctx0`: `builds/20261012-1401_01206ca/results/f9c82b3/ctx0`;
- `ctxin`: `builds/20261012-1401_01206ca.ctx5in/results/f9c82b3/ctxin`;
- `ctxfb`: `builds/20261012-1401_01206ca.ctx5/results/f9c82b3/ctxfb`.

## R1. Neither arm supports the claim

`ctx0` reproduces `results/64d23cb` on every curator and experimental
prediction.

**How the arms are paired.** `holdout_report.py` pairs a case only when both
arms pinned the same number of nodes. A context arm pins the gene in its
neighbours too, by design, so that rule dropped 6,724 curator cases, the very
ones context can move. The spec registered pairing on cases valid in both,
so the report gained `--pair-on-valid` (tested), and the figures below use
it: `--cases ctx0 --compare <arm> --set-aside RAF_MAP_kinase_cascade
--pair-on-valid`.

| | `ctxin` (feed-in) | `ctxfb` (feedback) |
|---|---|---|
| curator held-out (18,481) | **−164** (35 / 199, p < 1e-4), 22 of 70 pathways | **−577** (39 / 616), 24 of 70 |
| curator tuning (5,030) | −37 (67 / 104, p 0.006) | −423 (204 / 627) |
| experimental, all (845) | −14 (15 / 29, p 0.049) | −17 (41 / 58, p 0.11) |
| experimental, RAF set aside (796) | **−14** (15 / 29, p 0.049) | **−30** (28 / 58, p 0.0016) |
| valid only under context | curator 8 (5 correct), experimental 0 | curator 8 (3 correct), experimental 0 |
| experimental solves converged | 194 / 244 (`ctx0` 215) | 182 / 244 |

Gate: both fail the first condition (experimental net > 0) and the second
(curator held-out loss is significant). C-P1 was not separately checked
because the decision does not depend on it. Nothing is adopted, and the
korkut2015 scoring is not run.

## R2. Where the cases moved (`ctxin`)

**Experimental fixes (15)** are where C-P3 predicted: S phase (CCND1, CDK4),
Mitotic G1 (AKT1, PPP2R1A), checkpoints (ATM), and PIP3 (TP53). Not traced
individually.

**Experimental losses (29)**: 25 are three knockouts, all reading a
DOWN-expected readout UP (railed at 100×): EGFR in PIP3 (8) and ERBB2 (5),
PDPK1 in PIP3 (7), ERBB3 in ERBB2 (5). Two traced
(`bench/analysis/056/trace_bridges.jl`):
- **EGFR knockout, PIP3.** In Signaling by EGFR (R-HSA-177929) the knockout
  leaves the free ligands EGF (R-HSA-179863), EREG(63-108) (R-HSA-1233233)
  and four others unbound, and they read **10×**, the depletion
  de-repression ceiling (specs/011). The bridges carry them into PIP3's
  ligand roots, where they drive the other ERBB receptors and AKT rails.
  Alone, PIP3 holds those roots at baseline.
- **PDPK1 knockout, PIP3.** In Signaling by Rho GTPases (R-HSA-194315) the
  knockout drains free KDM1A (R-HSA-996769) to **0**; that copy is released
  only (its producers all release it from a complex). The bridge carries 0
  into PIP3's KDM1A root, and AKT rails.

Both import a neighbour's free form as if it were supply, and in both the
value is one the model reads at a bound, not a curated amount.

**Curator losses**: of 303 broken (held-out and tuning), 227 were expected
NO CHANGE and now change (C-P4, the July mechanism). Three traced:
- DSB repair, BRIP1 knockout: in Transcriptional Regulation by TP53 the
  knockout lowers free ATM, ATR, MDC1, CHEK1 and RAD51D to 0.25; bridged into
  DSB repair's roots, the readout goes to 0. The free-form drain class again.
- SCF-KIT, JAK2 knockout: in IL-4/IL-13 signaling (R-HSA-6785807) MMP9, a
  JAK–STAT transcription target, falls to ~0; bridged into KIT's MMP9 root,
  the readout reads 0.25. Plausibly real biology the curator, predicting
  from inside KIT signaling, would not mark.
- RUNX2, CCND1 overexpression: CCND1 is also pinned in IL-4/IL-13, which
  raises TWIST1 5×; bridged into RUNX2, the readout reads 1.17, just over
  the 1.15 cutoff.

**Is there one rule?** No. A release-only filter (the specs/053 "closed"
copy, 238 of 1,040 bridges) would remove the KDM1A bridge but not EGF, whose
copy is fed. Any filter chosen now would be fitted to these cases; none is
registered.

## What this means

- A third measurement of cross-pathway context (2026-05 joint networks,
  2026-07 stitch, this) loses on both axes. This time the hubs were excluded,
  the overlapping reactions removed and no cycles allowed (`ctxin`).
- The cell-cycle fixes are real, and they are the cases Adam had in mind.
  But on the experimental axis they are outnumbered two to one. The losses
  come from importing a neighbour's free forms, whose levels sit at the
  model's bounds (the 10× de-repression ceiling, drains to 0) rather than
  carrying curated supply.
- The feedback variant is worse everywhere except RAF, where 12 RAS
  overexpression cases flip right. RAF is set aside (specs/054).
