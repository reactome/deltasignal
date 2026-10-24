# specs/094 research: arm lel (part B)

Canonical `20261022-2011_233665c`, solver `231eee9` (main's code) with
`--server DS_LOOP_ELASTICITY=0.5 DS_LOOP_ELASTICITY_WIDTH=0.2 DS_LOOP_ELASTICITY_HI=0.95`,
against `results/8162f5f/dflt091`. Results `results/231eee9/lel`. Read with
`holdout_report.py --set-aside RAF_MAP_kinase_cascade --pair-on-valid`.

| split | net (fixed / broken) | p |
|---|---|---|
| curator held-out | **−40 (8 / 48)** | < 1e-4 |
| curator tuning | **−126 (4 / 130)** | < 1e-4 |
| experimental (RAF set aside) | **−42 (4 / 46)** | < 1e-4 |

By pathway (curator): **PIP3 −112**, MHC class I −19, TP53 −12, RUNX2 −10, intrinsic apoptosis −6,
VEGF −5, caspase activation −4; DSB Repair +4, Fanconi +3. Experimental: PIP3 −41, TP53 −1.
Every solve converges (1,725 / 1,725, 244 / 244). **The gate fails** on every count.

## Predictions

1. PTEN overexpression in PIP3, at least 6 of 9 experimental cases fixed: **refuted, 0**. It
   still reads UP (66.4×, one readout 100×) instead of 83.5×.
2. TP53 in PIP3, at least 2 of 4 fixed: **refuted, 0**.
3. Against specs/014 (held-out +29, tuning −61 on 2026-09-19 networks): held-out now −40.

## What broke

PIP3's experimental losses are receptor knockouts that were right at about 0.17× (DOWN) and now
read NORMAL: ERBB2 11, EGFR 8, KIT 7, ERBB3 7, PDGFRB 7, FGFR1 3. Their signal runs through
PIP3's large cyclic component, where every in-loop activator edge reads `fold^ε`; over many hops
the exponents compound, and a genuine knockout is compressed into the NORMAL band. The same
damping does not stop the PTEN loop railing, so the knife-edge is not what decides PTEN's
direction on today's network.

## Reading

Refuted on today's networks, more strongly than on September's. A constant loop elasticity
compresses real signals through large components and leaves the PTEN inversion in place.
specs/014's conclusion stands: if loop gain is the lever, it has to be per loop (learned), not a
constant. PTEN needs its own trace: the inversion survives a damped loop, so a second mechanism
(specs/083 noted a two-copy PTEN mRNA) is likely involved. Not traced here.
