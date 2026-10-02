# specs/041 research

## Step 1: route classes of the experimental axis (2026-09-29, summarised 2026-10-02)

**Script:** `~/deltasignal-catalogs/analysis/041/routes041/routes041.py`. The
route definition is in its docstring and was fixed before any accuracy was
read:
- **Graph:** 16,344 human reactions (Reactome 97). Inputs, catalysts and
  regulators feed a reaction; outputs are produced.
- **Moves along a route:** set membership hops are free. A complex is reached
  only by the reaction that forms it.
- **Excluded mid-route:** inert cofactors and free ubiquitin-family modifiers.
  They may still start or end a route.
- **Depth:** at most 12 reactions.

The variant without that exclusion is kept as `exp_routes_noblock.tsv`.
Predictions are canonical: build `20260928-1110_06ccb63`, solver `8d6d27c`.

| class | cases | right | wrong | accuracy |
|---|---|---|---|---|
| **a. in-pathway** | 752 | 570 | 182 | **75.8%** |
| b. only via other pathways | 76 | 21 | 55 | 27.6% |
| c. not in Reactome (no route ≤ 12; 3 readouts absent) | 21 | 2 | 19 | 9.5% |

- **74 of the 256 wrong experimental cases (29%) are outside the
  single-pathway model's reach** (b or c). That is more than the whole 50-case
  gap to MP-BioPath (643 vs 593); their hand edits added routes (specs/034
  §12d).
- **Mitotic G1/S is mostly scope:** 49 of its 89 cases are b or c, and 39 of
  its 54 errors.
- The class-b routes run mostly through FOXO-mediated transcription of
  cell-cycle genes (16), Cyclin D events (7), HIF regulation (7), Cyclin
  A/B1/B2 events (6), Fanconi Anemia (6) and PTEN localisation (6).
- **The fair score for the model is class a: 75.8%.** The 182 class-a errors are
  the ones to trace: TP53 75, PIP3/AKT 32, RAF 16, G1/S 15, WNT 13, HDR 12,
  Checkpoints 8.
- **Class b** is a test for the cross-pathway catalog. **Class c** is a
  curation-gap list (21 cases) for Adam.
- **Caveats:**
  - Reactome 97 is not the 2019 release that produced the ground truth.
  - The depth bound and the "no mid-route complex hop" rule are choices; the
    variant table shows the sensitivity to the ubiquitin/cofactor exclusion.
  - A route existing does not mean the curated mechanism produces the
    *expected sign*.

## Step 2a: blind drafts reconciled (2026-10-02)

- **Drafts:** two independent blind drafts, Fable (132 rows) and Sonnet
  (84 rows), each written without access to the benchmarks or to the other.
  Opus did not draft, having seen RAF's experimental expectations earlier.
- **Merged:** 167 distinct (gene, KO/80x, readout) expectations after
  normalising readout names; RAF/MAPK 77, PI3K/AKT 90.
  - 49 rows appear in both drafts: **42 agree** (86%), 7 are disputed.
  - **No dispute is a flat UP vs DOWN contradiction.** Six are definite vs
    uncertain (BRAF 80x on p-MEK/p-ERK, SOS1 KO/80x on RAS:GTP, PPP2CA KO,
    RICTOR KO). One is a real biology question: RAF1 KO on p-ERK, NONE
    (BRAF compensation) vs DOWN.
  - 118 rows come from one draft only.
- Files: `~/deltasignal-catalogs/analysis/041/bio041/reconciled.tsv`.
- **Adam's review** is collected on a private review page; verdicts are stored
  in its `verdicts` collection. **The model is not run on these expectations
  until the review is done**, so the review is blind to the model's
  predictions.
