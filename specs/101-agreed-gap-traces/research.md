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
