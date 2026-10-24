# specs/101: traces of the remaining agreed experimental errors (analysis)

Canonical build `20261024-0801_7a754fc`, solver `b0d6737` (`results/b0d6737/dflt100`):
experimental 639/849, MP-BioPath 643. Classes where curators agree with the experiment and
MP-BioPath is right (specs/094): NF1 KO in RAF (5; RAF set aside, specs/054), RB1 KO in
Mitotic G1 (4), MDM4 KO in TP53 (3), AKT1 KO in TP53 (2).

## T1. RB1 knockout in Mitotic G1: not curated as a route

- Our prediction is exactly baseline for every G1 readout (pred_ui 1.000). Curator: 30 of
  50 RB1 cases wrong (KO should read UP, OE DOWN; tuning pathway). Experimental: 4 wrong
  (CCNA1, TK1, CCNE1 UP; FBXO5 DOWN).
- Network: RB1 [cytosol] (root, pinned) → translocation → RB1 [nucleoplasm] → (a) CDK4/6
  phosphorylation (R-HSA-69227) and (b) "RB1 binds and inhibits E2F1/2/3:DP1/2 complexes"
  (R-HSA-9018017), whose product RB1:E2F (R-HSA-68644) is a dead end.
- Neo4j (Release 97): R-HSA-68644 is no input, catalyst or regulator anywhere. The E2F set it
  binds (R-HSA-1227905) is used only by 9018017 and four promoter-binding reactions (p14ARF,
  EED, EZH2, SUZ12). The G1 target genes (CDC6, PCNA, POLA1, TK1, ...) are activated by
  E2F1:TFDP1,TFDP2:gene complexes, another entity. No GO process marks 9018017. The
  repression lives only in the reaction's name.
- MP-BioPath's 2019 network gets it through hand-added edges: RB1 (R-HSA-68642) inhibits 15
  transcription reactions (R-HSA-8961619 ...). Our propagator on that network reads 10× UP on
  every E2F target: 3 experimental right, but 3 curator-NORMAL readouts (PCNA, 68487, RRM2)
  and FBXO5 DOWN wrong.
- Not a rule candidate: it is the sequestration class (specs/063, refuted in three forms; no
  structural signature), and here even full sequestration would not reach the target genes
  (different E2F entity). **Curation-review item** (Adam's call): RB1:E2F as a negative
  regulator of E2F1-stimulated transcription in Mitotic G1.

## T2. MDM4 knockout in TP53: unbalanced copies of "MDM2 forms homo- or heterodimers"

See spec.md (pre-registration of `LNG_BALANCE_CHOICES`). Trace (`trace_mdm4.jl`, MDM4 KO,
folds against baseline): TP53 [nucleoplasm] 1.064 ← depletion by pool R-HSA-6804885 (0.855) ←
variants MDM2:MDM4 and MDM2 dimer both 0.855 ← MDM2:MDM4 (0.803) fed by two copies of
R-HSA-6804741: MDM4 + p-MDM2 (0) and p-MDM2 alone (1.607, makes MDM4 from nothing).

## T2 result 1: the `bal` build (structure)

Builds from LNG `balance-choices-101` 64a95a0: `20261024-1307_64a95a0_bal` (flag on) and
`20261024-1346_64a95a0_balctl` (flag off, control).

- R-HSA-6804741 in TP53: 4 copies → 2 (the two unbalanced pairings gone), as designed.
- `compare_build.py` canonical vs bal: 38 pathway tables differ. **Prediction 3 fails as
  written**: 3 moved pathways are not in the census (Signaling by FGFR2 R-HSA-5654738, ERBB4
  R-HSA-1236394, AUF1 R-HSA-450531) and 2 census pathways do not move (R-HSA-198933,
  R-HSA-5617472).
  - The extras are the intended class, which the census missed because it counted only
    output edges that gain a protein: set-to-set bindings emitted as a cross product. "FGFR2b
    binds to FGF" (R-HSA-190260) 240 copies → 120, FGFR2c 48 → 24, "ERBB4 binds NRGs"
    (R-HSA-1236398) and "Homodimerization of ERBB4" (R-HSA-1250220) 63 → 21 each, AUF1
    ubiquitination (R-HSA-450580) 4 → 2.
  - The two census pathways that do not move have unbalanced copies with no better-balanced
    partner on both sides, so the rule keeps them.

The control build `balctl` (flag off) is identical to canonical: 0 of 92 tables differ. The
canonical defaults arm `20261024-0801_7a754fc/results/b0d6737/dflt100` is therefore the control.

## T2 result 2: arm `bal` (`20261024-1307_64a95a0_bal/results/319dd0e/bal`, solver = main ddf0060)

Read with `holdout_report.py --cases <dflt100> --compare <bal> --pair-on-valid --set-aside
RAF_MAP_kinase_cascade` (`--cases` is the baseline; a first reading had the arms swapped).

| axis | net | fixed / broke | p | RAF set aside |
|---|---|---|---|---|
| curator held-out | **+30** | 34 / 4 | < 1e-4 | +30 (no RAF case moves) |
| curator tuning | −12 | 11 / 23 | 0.058 | −12 |
| experimental | −4 | 2 / 6 | 0.29 | −4 |

Curator all 88.40% → 88.48%; held-out in release 90.25% → 90.41%; experimental 639 → 635/849.
**Registered gate passes**: held-out +30 ≥ 0, experimental −4 ≥ −5, worst pathway
Signaling by ERBB2 −7 curator / −3 experimental (≥ −10).

Predictions:
1. **Holds.** MDM2:MDM4 (R-HSA-6804932) reads 0 under MDM4 knockout.
2. **Fails.** 0 of the 43 wrong MDM4 knockout cases fixed; 8 more broken (curator), experimental
   1 fixed / 2 broken. Trace (`trace_bal.jl` = `trace_mdm4.jl` on the bal build): with the
   heterodimer gone, p-MDM2 is no longer drained into it and reads 8.3×; the MDM2 homodimer
   (6804933) reads 5.5–100×; the ligase pool R-HSA-6804885::pool reads 2.73× (was 0.855), and
   TP53 [nucleoplasm] falls to 0.655. MDM4 knockout now reads TP53 targets DOWN (0.12–0.43).
   The network says MDM2 alone is as good a TP53 ligase as MDM2:MDM4 (one output set), so
   freeing MDM2 from MDM4 raises ligase. The curators' answer (MDM4 knockout de-represses
   TP53) needs the heterodimer to be the stronger ligase, which Reactome does not encode.
   The 3 MDM4 cases fixed are overexpression readouts that now read DOWN (0.81).
3. **Fails as written** (result 1).

What moved, by mechanism (curator):
- Fixed, Interleukin-4/13 (24, all STAT6/STAT3 knockout/overexpression at NORMAL readouts):
  "p-Y-STATs dimerize" (R-HSA-6788622, 9 → 3 copies) and R-HSA-6786097 (4 → 2) no longer
  let p-STAT6 make the STAT3 or STAT1 dimer; the readouts are other STATs' targets and
  now read exactly baseline, as the curators say. Same shape for STAT5B in IL-3/5/GM-CSF (7),
  ROBO2 (2), PDGFRA (1).
- Broken, Signaling by ERBB2 (15, of which 8 fixed): "ERBB2 forms heterodimers with
  ligand-activated EGFR/ERBB3/ERBB4" (R-HSA-1963589, 227 → 25 copies) and the two
  trans-autophosphorylations. EGFR and ERBB3 knockouts used to feed every heterodimer copy;
  now each lowers only its own heterodimer and the readouts move to ~0.94–0.96 (DOWN expected
  by curators and experiments). Right before only through the cross product.
- Broken, TP53 (8, MDM4 above); IFN α/β PTPN11 (2); TGF-beta SMAD2/3 (2, 0.60 → 0.94).

Status: structural fix that passes its registered gate on held-out gains that look right
case by case; its motivating case (MDM4) inverts instead of improving. Adoption is Adam's call.
