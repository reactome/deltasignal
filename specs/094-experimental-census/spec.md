# specs/094: the remaining experimental errors, and the PIP3 knife-edge re-measured

**Status:** part A analysis, 2026-10-23. Part B pre-registered 2026-10-23, before the arm; arm `lel` **refuted**: held-out −40, tuning −126, experimental −42 (PIP3 −112 / −41: receptor knockouts compressed to NORMAL), and PTEN overexpression still reads UP. Not adopted.

Adam, 2026-10-23: "okay continue" to a census of the remaining experimental errors.

## A. Census (current canonical defaults, `results/8162f5f/dflt091`)

`specs/079-error-census/census.py` and `agreed_errors.py` (→ `agreed_errors.txt`):

- experimental: 845 valid, **220 wrong** (no_path 98, propagator_missed 83, false change 39);
- the curators agree with the experiment on 693 of 845. Where they agree we are right on 615
  (88.7%); where they disagree (152) we give the curators' answer 135 times, the experiment's 10;
- so **142 of the 220 errors are curator–experiment disagreements**, and **78** are errors on
  which every reference agrees: no_path 26, propagator_missed 46, false change 6.

The 78 by perturbation (largest first): **PTEN overexpression in PIP3 9** (all read 83.5× UP,
expected DOWN), NF1 knockout in RAF 5, RB1 knockout 4, AKT1 knockout in Mitotic G1 3, SRC
overexpression in Mitotic G1 3, RAF1 knockout 3, MDM2 knockout in Cell Cycle Checkpoints 3,
AKT1 knockout in TP53 3, MDM4 knockout in TP53 3, and 2 or 1 each for 33 more. With TP53 in
PIP3 (overexpression reads 83.5–100× UP, knockout reads DOWN: 4 more), **the PIP3 PTEN axis is
13 of the 78**, the largest coherent class.

specs/082 traced it: the curated loop AKT → MKRN1 ⊣ PTEN ⊣ PIP3 → AKT has gain exactly 1 at
baseline (every fold-multiplying hop has elasticity 1), so any perturbation rails it, and which
rail does not follow PTEN. That is specs/014's knife-edge.

## B. Pre-registered: specs/014's loop elasticity on today's networks

specs/014's best reproducible configuration (`DS_LOOP_ELASTICITY=0.5`,
`DS_LOOP_ELASTICITY_WIDTH=0.2`, `DS_LOOP_ELASTICITY_HI=0.95`: in-loop activator edges read
`fold^ε(f)`, ε rising from 0.5 near baseline to 0.95 above the band, so no loop has gain 1)
gave held-out +29, tuning −61 on 2026-09-19 networks; the default was left open. Since then
the networks and the solver changed a great deal (cycle-balance pools, carrier loops,
specialization bridges, the non-converged fallback), so the trade-off is re-measured, unchanged.

| arm | build | solver | against |
|---|---|---|---|
| `lel` | canonical `20261022-2011_233665c` | main `e1d2fc2`, `--server DS_LOOP_ELASTICITY=0.5 DS_LOOP_ELASTICITY_WIDTH=0.2 DS_LOOP_ELASTICITY_HI=0.95` | `results/8162f5f/dflt091` (same code defaults) |

**Gate (as specs/077):** curator held-out net ≥ 0, experimental net ≥ −5, no pathway below −10.

### Predictions

1. PTEN overexpression in PIP3: at least 6 of the 9 experimental cases are fixed (the loop no
   longer rails).
2. TP53 in PIP3: at least 2 of its 4 experimental cases are fixed.
3. No prediction on held-out or tuning nets. specs/014 measured held-out +29, tuning −61; the
   direction is reported against that.

## C. The gap to MP-BioPath, split by reference agreement (2026-10-23)

Adam, 2026-10-23: "I was hoping we would get the experimental results to at least equal to
where we were with the publication." At `results/b61b365/dflt095` DeltaSignal scores 634 of 849
against MP-BioPath's published 643 (`bench/analysis/experimental_gap.py` → `gap_b61b365.tsv`;
`gap_split.py` → `gap_split.txt`):

| | MP-BioPath right, DeltaSignal wrong | DeltaSignal right, MP-BioPath wrong | net |
|---|---|---|---|
| curators agree with the experiment | 35 | **46** | **DeltaSignal +11** |
| curators disagree with the experiment | **26** | 6 | MP-BioPath +20 |
| total | 61 | 52 | MP-BioPath +9 |

On the cases where every reference agrees, DeltaSignal is ahead. The whole gap is on
curator–experiment disagreements, where MP-BioPath (whose networks carry ~1,800 hand edits made
against experimental data, specs/034) matches the experiment and DeltaSignal mostly gives the
curators' answer. Under specs/096 those are not targeted, except where curated biology supports
the experiment (a justified disagreement with the curators): CDKN2A in TP53 (ARF sequesters MDM2,
not curated in that pathway) is the first candidate. The agreed classes MP-BioPath still wins:
NF1 knockout in RAF 5, RB1 knockout 4, MDM4 knockout 3, TP53 in PIP3 4, MDM2 knockout in Cell
Cycle Checkpoints 2, AKT1 knockout in TP53 2.
