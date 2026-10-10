# 047: Where the experimental gap to MP-BioPath is, case by case

**Status:** analysis in progress (2026-10-09). No behaviour change yet.

## Why

On the experimental axis we score 594/849 (69.96%), and MP-BioPath scores
643/849 (75.74%). Curator-axis work has gained 0.3–0.7pp per feature while
experimental drifted from 627 to 594 across generator changes. Experimental
is the axis measured against biology, so it should be the target, not a
guardrail. This spec decomposes the gap before any lever is designed.

## Data

MP-BioPath's per-case calls are published in
`mp-biopath-pathways/analysis_results/MP_BioPathReactomePathwayAccuracy-Tests-updated.tsv`
(`mp_biopath_state` against `experimental_value`).
`bench/analysis/experimental_gap.py` joins them with an arm's
`experimental_cases.tsv`:
- all 849 cases join;
- the experimental value equals our expected label on all 849;
- 643 MP-BioPath-correct reproduces.

A missing case or a label disagreement is an error, not a skip.

## Decomposition (canonical: build `20261004-1959_6990015`, solver 5979e48)

| | Cases |
|---|---|
| Both right | 540 |
| **Only MP-BioPath right** | **103** |
| Only DeltaSignal right | 54 |
| Neither | 152 |

The two are complementary: either one right covers 697/849 (82.1%).

**Per pathway** (only-MP-BioPath | only-DeltaSignal | net):

| Pathway | Only MP-BioPath | Only DeltaSignal | Net |
|---|---|---|---|
| TP53 | 32 | 18 | +14 |
| Mitotic G1 | 20 | 11 | +9 |
| PIP3 | 18 | 9 | +9 |
| Cell Cycle Checkpoints | 8 | 1 | +7 |
| RAF | 9 | 3 | +6 |
| HDR | 4 | 0 | +4 |
| Mitotic Prophase | 2 | 0 | +2 |
| S Phase | 1 | 1 | 0 |
| ERBB2 | 3 | 4 | −1 |
| WNT | 6 | 7 | −1 |

The gap is spread across pathways, not one.

**The 103 only-MP-BioPath cases:**
- **By our failure category:** propagator_missed 47, no_path 31,
  false_positive_change 24, readout not in network 1.
- **By (truth, our call):**

  | Truth → our call | Cases |
  |---|---|
  | DOWN → UP | 30 |
  | UP → NORMAL | 22 |
  | DOWN → NORMAL | 15 |
  | NORMAL → UP | 14 |
  | UP → DOWN | 11 |
  | NORMAL → DOWN | 10 |

  **41 are sign inversions**: the signal reaches the readout backwards.
- **By perturbation (top):**

  | Perturbation | Cases | Note |
  |---|---|---|
  | PTEN UP (PIP3) | 9 | all inverted |
  | RBL2 KO (Mitotic G1) | 7 | |
  | PDPK1 KO (PIP3) | 7 | inverted |
  | DAXX KO + UP (TP53) | 11 | |
  | CDKN2A UP + KO (TP53) | 8 | |
  | RB1 KO | 4 | |
  | NF1 KO | 4 | |
  | BRCA1 KO | 4 | |

## Trace 1: PTEN overexpression inverts (9 cases)

The readouts read 100× (UP) after PTEN overexpression; truth is DOWN. They
also read 100× after PTEN knockout, so they do not depend on PTEN's direction.

1. The bench pins all three PTEN roots at 80×: PTEN gene (R-HSA-5632949),
   K48polyUb-PTEN (R-HSA-8948777) and PTEN mRNA (R-HSA-2318746).
2. PTEN translation, R-HSA-8944497 ("PTEN mRNA translation is negatively
   regulated by microRNAs"), has **12 negative regulators**. Each is a
   `miR-x RISC:PTEN mRNA` complex, and each reads 64.5×.
3. Under `divide` they multiply: 80 / 64.5¹² ≈ 1.5e-20. **Overexpressing PTEN
   shuts off PTEN protein.** The PTEN copy that depletes PI(4,5)P2 reads 0,
   PIP2 is de-repressed to 6.7×, PIP3 rises, and a positive loop rails AKT to
   100×.
4. The specs/022 self-inhibitor rule is on (227 slots flagged in the
   pathway), and the complexes contain PTEN mRNA. The rule still does not
   fire, because it requires the shared input node to reach the inhibitor.
   PTEN mRNA exists as **two silo copies**:
   - root bb6f1223 feeds only translation;
   - produced copy 337e82aa (from transcription of the pinned gene) feeds the
     12 RISC-binding reactions.

   Translation reads one copy, its inhibitors are built from the other.
5. With only the mRNA pinned, the readout reads 0.00065 (DOWN, correct).

Biologically, RISC binding sequesters a fraction of the mRNA. More mRNA
means proportionally more RISC:mRNA, not less translation. The failure is
the specs/012/022 double count, unmasked by copy siloing.

**Candidate levers** (not designed yet):
- (a) the self-inhibitor rule matches the shared input by entity across copies
  rather than by node reachability;
- (b) the generator joins the transcription product to the translation input
  (the copy silo).

Either is a pre-registered arm, measured with experimental as the primary
axis.

## Trace 2: PDPK1 knockout inverts through a SARS-CoV-2 inhibitor (7 cases)

PDPK1 knockout should lower phospho-AKT; truth is DOWN, and we read 100×.

1. "PDPK1 phosphorylates AKT at T308" (R-HSA-198270) has a negative
   regulator, `N:M:PDPK1` (R-HSA-9755778): the SARS-CoV-2 N dimer and M
   protein bound to PDPK1. Reactome added it after 2019, so MP-BioPath's
   networks lack it.
2. Knocking out PDPK1 zeroes that complex, the reaction is de-repressed to the
   ceiling, and the AKT positive loop
   (R-HSA-198270 ↔ R-HSA-2317313 ↔ R-HSA-2317314) rails to 100×.
3. In a human cell without the virus the inhibitor does not exist. Holding it
   at baseline makes it a constant factor, which cancels in fold-change. That
   is the specs/032 drug rule's mechanism applied to pathogen-derived nodes.

Across the ten experimental pathways, the containment tables name only three
non-human entities, all SARS-CoV-2: N dimer, a modified N, and M. A pathogen
rule would be narrow.

Separately, the catalyst `p-S-AKT:PDPK1:PIP3` (R-HSA-2317313) reads 100× under
a PDPK1 knockout. Its PDPK1 comes in through the loop, not from the pinned
copy. That is the same copy-silo pattern as trace 1, not yet traced to a
line.

**Candidate lever:** `DS_PATHOGEN_MODE=inert`, or the generator listing
non-human-derived nodes the way `drugs.csv` lists drugs. Pre-registered and
measured with experimental as the primary axis.

## Controls running

DeltaSignal's propagator on MP-BioPath's own 2019 networks (build
`20261009-2051_mpbnet2019`, via the adapter), against canonical. Both arms
run with `DS_CYCLE_MODE=off` and `DS_SELF_INHIBITOR_WEIGHT=off`, because
MP-BioPath's networks carry no pool or containment tables. That splits the
49-case gap into a network share and a propagator share.
