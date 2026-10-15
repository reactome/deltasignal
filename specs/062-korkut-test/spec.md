# specs/062: korkut2015 as the independent check of specs/061's candidate

**Status:** registered 2026-10-15, before any korkut2015 value is read
(`expr.csv`, `pert.csv`, and supplementary tables 2–5 stay unread). Adam
(2026-10-15) asked for both independent checks; the other one, the next
Reactome release, waits: Release 98 is not published (ContentService reports
97, and no `graphdb:Release98` image exists).

## Why

specs/061 `rtwin58g2` (specs/058 pins + specs/059 twin edge + the amended
own-product rule) passes the registered gate on the curator and experimental
axes, but it was refined four times against those sets. korkut2015 (specs/051
R2) is the registered TEST set: one melanoma cell line (SkMel-133), 12 drugs
alone and in pairs, 82 RPPA protein readouts, values never read.

**Role.** TEST is "scored once per frozen release; never used to choose". The
candidate was chosen on SELECT and is frozen (`9a0b839`); this scores it and
the control once, with the verdict below fixed in advance. That confirms or
refutes one choice already made; it selects nothing among alternatives. After
this scoring the set is no longer TEST: it is re-registered FIT, as decryptM 1
was (specs/051).

## Design (from names and supplementary table 1 only)

- **Drugs** (`bench/analysis/062/drugs.tsv`): each code mapped to its drug by
  dose (every code's two doses in `expr_index.txt` equal supplementary table
  1's), then to target genes. Nutlin blocks MDM2's binding to TP53; it is
  modelled as MDM2 down (an approximation, stated here).
- **Readouts** (`bench/analysis/062/readouts.tsv`): 82 antibodies; 29 phospho
  (gene and coordinate, legacy numbering translated, checked against Release
  97 with `korkut_sites.cypher`), 44 total protein, 9 excluded (6 sites not
  curated in Release 97, the loading control, two ambiguous names).

## Protocol (bench `bench/korkut_bench.py`, mirroring specs/021's drug bench)

- A **condition** is one row of `expr_index.txt`: one drug, or a pair, at given
  doses. A **drug set** is its codes without doses.
- **Solve:** per drug set, per catalog pathway where any target gene has a
  root form: every target gene of every drug in the set pinned to
  `DS_PERTURB_UI_DOWN` at the forms the curator bench pins
  (`pin_occurrences`, so `DS_PIN_ROOT_IGNORE` applies as in the arm).
- **Prediction:** the mean over the readout's nodes (phospho: nodes that are or
  contain an EWAS of the gene phosphorylated at the coordinate; total: any EWAS
  of the gene), pinned nodes left out; classified with the curator cutoffs.
  A case is (drug set, pathway, readout); one with no readout node, or only
  pinned ones, is not scored.
- **Measured:** `expr.csv` is taken to be log2(treated / untreated), as CellBox
  documents; if the file shows otherwise, stop and amend before scoring. Per
  condition: UP if ≥ log2(1.5), DOWN if ≤ −log2(1.5), else NOT. Per drug set,
  over its conditions (doses): all regulated calls agree → that call; regulated
  calls disagree → `conflict` (not scored); none regulated → NOT.
  Secondary threshold, reported only: log2(1.15).
- Both arms run through `scripts/run_arm.sh --axis korkut` on the canonical
  build: control (defaults) and `rtwin58g2` (the specs/061 overrides), same
  code commit.

## Step 1, values-blind: is the test informative?

Before `expr.csv` is opened, both arms run with `--axis korkut_predict`
(predictions only). If fewer than **10** scored cases change predicted class
between the arms, the test cannot inform the decision: it is recorded as
uninformative, and **the values stay unread** (the set stays TEST for a later
question). The cases that change are listed (drug set, pathway, readout).

## Step 2, the verdict (only if step 1 passes)

On cases scored in both arms where the predictions differ, exact two-sided
McNemar on (control right, arm wrong) vs (control wrong, arm right):
- **confirms** if the arm's net > 0 at p < 0.05;
- **contradicts** if net < 0 at p < 0.05;
- **inconclusive** otherwise.
Also reported, not part of the verdict: both arms' accuracy and macro-F1 on all
scored cases, the always-DOWN and reachability baselines (as specs/021), the
result at the secondary threshold, and the verdict clustered by drug set (one
cell line: readouts within a drug set are not independent).

The verdict is evidence for Adam's adoption decision, not an automatic one.

## Predictions (before step 1)

- Few cases will differ. The candidate moves TP53 (AKT, MDM4, MDM2 knockouts),
  TGF-beta SMAD2/3, MET, CBL, PTCH. korkut2015 perturbs AKT (AKTi) and MDM2
  (Nutlin) and reads p53, p21, AKT and SMAD3 forms; the overlap is mostly the
  AKTi and Nutlin drug sets in TP53 pathways. I expect step 1 to be close to
  its threshold, and I would not be surprised if it fails.
