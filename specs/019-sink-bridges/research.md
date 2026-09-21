# Research: acyclic sink bridges

## Findings before the rule was built (fixed catalog `cat_fix2`, 2026-09-20)

Error anatomy: 24,100 cases, 3,686 wrong (15.3%). By class: **missed change /
`no_path` 1,526**; false change 1,150; propagator missed 364 + wrong
direction 309; gene or readout not in network 337. Worst pathways: Mitotic G1
520/884 wrong (430 `no_path`), Interferon α/β 354/448 (256 `no_path`).

Oracle on `cat_fix2`: Mitotic G1 957 stable-id-connected root→terminal pairs,
530 uuid-connected, **430 severed** (= its `no_path` count); break kinds
dissociation_sink 264, simple_entity 145, reaction 17. Interferon α/β 642 →
293, 392 severed; break kinds reaction 266, dissociation_sink 78 — the
"reaction" breaks sit downstream of the sink break at R-HSA-909683 →
R-HSA-9710965.

Simulation (bridges added to copies of the bundles, oracle re-run):

| | sinks | acyclic bridges | cycle-closing skipped | fan-out median / p90 / max | uuid-connected | cyclic comps |
|---|---|---|---|---|---|---|
| Interferon α/β | 30 | 13 | 23 | 3 / – / 3 | **293 → 631 of 685** | 2 → 2 |
| Mitotic G1 | 287 | 234 | 152 | 1 / – / 19 | **530 → 705 of 960** | 4 → 4 |
| **catalog (92)** | 35,230 | **10,547** | **56,259** | 1 / 8 / 160 | severed **47.5% → 37.9%** (18,303 → 14,615 of 38,550) | **210 → 210**, largest unchanged |

Five of six candidate bridges would close a cycle; the guard is the rule. The
remaining severance is `simple_entity` and reaction-level breaks, not sinks.

## Pre-registration (committed before the regeneration)

Regenerate `cat_fix2_sb` (LNG `2387c88`, `LNG_SINK_BRIDGES=1`) and
`cat_fix2_sb8` (fan-out cap 8). Arms from a pinned DS worktree (`main`
`72015ce`), production solver, same-catalog control via
`DS_SKIP_EDGE_TYPES=sink_bridge` (zero churn), curator and experimental.

- **P1 (structure)**: cyclic components 210 in every catalog; bridges within
  10% of the simulation (10,547; cap-8 fewer); oracle severance ≈ 37.9%.
- **P2 (the two pathways)**: Interferon α/β ≥ +150 (it has 256 `no_path`
  cases and its connected routes go 293 → 631); Mitotic G1 ≥ +100 (430
  `no_path`, routes 530 → 705).
- **P3 (held-out)**: ≥ +100 vs control, p < 0.05; Interferon α/β will
  dominate it and is stated as such; Mitotic G1 is tuning.
- **P4 (the risk)**: new routes carry perturbations to readouts the curators
  call NORMAL — false change rises, but by less than half the gross gain; the
  cap-8 arm has a smaller false-change rise than the uncapped one and a
  held-out net within 30 of it.
- **P5 (experimental)**: conditioned net ≥ −10 (Mitotic G1 is experimental-
  heavy: this is where the earlier bridges lost).
- **Decision**: adopt the better of the two arms as the generator default
  only if P3, P4 and P5 hold; a P3 gain carried entirely by Interferon α/β
  with Mitotic G1 not moving is a one-pathway fix and is reported as such.
  If held-out < 0, the sink design stands as measured (specs/015) and the
  acyclic guard was not the missing piece.

## Side question (Adam, 2026-09-20): are the depletion edges necessary?

Pre-registered before the arm: `cat_fix2`, production solver, client-side
`DS_SKIP_EDGE_TYPES=depletion` vs the same catalog with them (zero churn), both
axes. Prediction: removing all 4,930 depletion edges is **negative** — held-out
≤ −40 and tuning ≤ −80, concentrated in TP53 (the AKT-KO cases run through
MDM2:TP53 ⊣ TP53; expect ≥ 60 of the 100 lost) and EGFR/GRB2 (specs/011's
case); false change *falls* (fewer routes) while missed change rises more.
If instead held-out ≥ 0, the edges are not load-bearing and the class is a
candidate for removal on fidelity grounds (Reactome never asserted them).

### Depletion ablation scored (2026-09-21; `cat_fix2`, same catalog, zero churn)

| | with depletion (control) | **without (all 4,930 edges skipped)** |
|---|---|---|
| all pathways | 84.71% / mF1 0.8137 | **83.52% / 0.7960** |
| held-out | | **−90** (38 fixed / 128 broke, p < 1e-4), 13 pathways, 43 readouts, 33 genes |
| tuning | | **−196** (TP53 −241) |
| experimental, conditioned | | **−24** (21/45, p 0.0043), 5 of 10 pathways |
| false change | 1,150 | **974 (−176)** |
| AKT-KO → TP53 (of 102) | 100 | **2** |
| predictions changed | | 815; 489 of them collapse to NORMAL (246 true DOWN, 243 true UP lost) |

Per pathway: TP53 −241, MET −60, IFN-γ −34, EGFR −25; only ERBB2 +26 and WNT +9 positive.

**Every clause of the pre-registration held.** Depletion edges are load-bearing
on both axes: removing them does exactly what was predicted — false change
falls (fewer routes) and missed change rises more than twice as far. They are
the only edge class that carries "consumption": an abundant complex drawing
down its free subunit, a phosphatase consuming its substrate. Without them a
knockout upstream of a complex cannot raise the free partner, which is how
93 of the 100 AKT-KO → TP53 cases are answered (AKT ⊣ MDM2 phosphorylation →
less MDM2:TP53 → more free TP53).

They remain **our inference, not curation** (1.5% of edges; Reactome asserts
no such relation), and the two bounded-ness defects found this month were both
in this class (the floor at zero, specs/011, +28 held-out; the own-product
double count, specs/016, inert). The honest statement for a paper: *a
consumption term is required for knockout propagation through complexes;
Reactome does not record it; we derive it from reaction stoichiometry and
bound it symmetrically.*
