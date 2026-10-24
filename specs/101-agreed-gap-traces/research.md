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
