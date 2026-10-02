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

## Step 3: topology motifs behind the 182 class-a experimental errors (2026-10-02)

**Method (Fable).** Every error was re-solved with the benchmark's own pins
(71 perturbations, 10 pathways; the harness reproduces 116 of 121 propagator
values exactly).
- 138 cases traced directly, 26 via a sibling readout, 18 by flags.
- Tables and traces: `~/deltasignal-catalogs/analysis/041/motif041/`
  (`motifs.tsv`, `counts.md`, `traces/`, scripts).
- Build `20260928-1110_06ccb63`, solver `8d6d27c`.

| motif | n | pathways | kind |
|---|---|---|---|
| **M4: a positive cycle with gain ≥ 1 rails; the TP53 hub self-loop** | **62** | 5: TP53 23, PIP3 25, HDR 8, ERBB2 3, WNT 3 | solver |
| **M3: the KO reaches a brake built from the pin, but not the recycled drive** | **27** | 3: RAF 16, PIP3 7, CCC 4 | solver |
| **M2: route severed in our network** (dissociation-sink copy; unexpanded produced set; pin-form choice) | **15** | 3: TP53 10, CCC 4, S 1 | generator |
| **M1: missing negative route** (a consumer represented only as a producer; branch competition; binder/repressor with a positive sink) | **14** | 4: WNT 8, HDR 4, S 1, G1 1 | generator |
| M9a: route-class artefact (see the correction below) | 39 | TP53 30, G1 8, ERBB2 1 | not a model error |
| M10: curated route followed, experiment disagrees (paralog redundancy, context) | 13 | G1 5, TP53 8 | not a model error |
| M9b: readout mapping / ground truth | 5 | | not a model error |
| M7: within 0.85–0.88 of the threshold | 4 | TP53 | not a model error |
| M5 (container pinned), M6 (AND) | 3 | 1 each | single pathway, no rule |

**The traced mechanisms:**
- **M4.** PIP3 sits in a 40-node cycle (PIP3 → PDPK1:AKT → PI3K). HDR's
  237-node resection machine has CHEK1 phosphorylation feeding its own input.
  ERBB2 has a 243-node cycle. WNT: a 1-of-19 ligand KO rails nuclear CTNNB1
  through 48 LRP6-phosphorylation outputs. TP53: the p-S15,S20 tetramer sits on
  a gain ≈ 1 self-loop, so every modulator moves every target 1.3–69x, and ATM
  KO even inverts.
  - Curator errors in the same pathways share the two shapes: PIP3 42 of 42
    wrong, HDR 74 of 78, ERBB2 59, WNT 115 of 178 railed, TP53 158 of 224.
- **M3.**
  - BRAF KO reads 63x UP: a drug/RKIP-bound RAF brake falls to 0 while the
    RAF-dimer set and its recycled scaffolds are OR-regenerated.
  - PDPK1 KO: a complex assembled from the pin inhibits the AKT step.
  - MDM2 KO: p14ARF:MDM2:MDM4.
  - specs/022's test does not fire because the pin is a catalyst/co-input,
    and x^(1−w) at x = 0 is still 0. Curator: RAF 37 of 37, PIP3 PDPK1 14,
    CCC MDM2 6.
- **M2.**
  - CDKN2A: the released TP53 tetramer copy dissociates to a dead end (the
    specs/015 sink).
  - CCNE1: a produced set is never expanded to the member that acts.
  - ATM: root_cycle pins the free p-ATM form, not the DSB-bound chain.
- **M1.**
  - WNT: the free CTNNB1 root has no depletion in-edge, and every other
    CTNNB1 copy is a *product* of destruction-complex reactions. So APC KO
    → 0 and 80x → 100, signs inverted. Curator APC/AMER1: 70 wrong, the same.
  - HDR: HR vs SSA branch competition.
  - p27 and RB1: binders/repressors with only a positive sink.
  - This is specs/034 §12d's "brake" class on our own networks.

**Correction to step 1.** The route search allowed free up-then-down set
hops. 39 errors have no route except such a walk, and nothing in the pathway
consumes where they stop: Ac-TP53 a dead end, PTEN reached through the set
"TP53,MDM2,MDM4,FOXO4,PTEN", p27 phosphorylated into a dead end.
- They belong in class b or c, not a.
- Out-of-reach errors are therefore at least 74 + 39 = 113 of 256 (44%).
  Counting the 13 M10 and 5 M9b cases, which are not model errors either,
  it is 131 (51%).
- **The four motifs in ≥ 3 pathways (M4, M3, M2, M1) account for 118 of the
  182** and are the work. M4 and M3 are solver shapes; M2 and M1 are generator
  shapes.

## Step 3b: M2, M3 and M4 checked against Reactome (2026-10-02)

Reactome 97 (DBInfo 2026-06-17). The build's reaction set matches it exactly
in all 8 pathways.
- **Every traced edge is a curated relation**, or a declared derived class
  (assembly, dissociation, depletion, set_member). There are 0 invented
  relations.
- Tables: `~/deltasignal-catalogs/analysis/041/verify041/` (`verify.tsv`,
  `trace_edges.tsv`, scripts).

| motif | curated as-is | generator artefact | mixed |
|---|---|---|---|
| M2 (15) | 4 (CCC ATM: pin-form/scope) | 11 (TP53 CDKN2A 10: tetramer copies joined only along precedingEvent; S CCNE1 1: output set dissolved to leaves) | 0 |
| M3 (27) | 0 | 4 (CCC MDM2: `_leaf_uuid` counts depletion edges as reachability and splits the pin) | 23 (RAF 16: curated drug and PEBP1 brakes, plus a variant-cap bundle making an OR set 28 AND inputs; PIP3 PDPK1 7: a curated SARS-CoV-2 brake) |
| M4 (62) | 34 (TP53 23, HDR 8, WNT 3: the loops are Reactome's own) | 0 | 28 (PIP3 25, ERBB2 3: curated cycles that close through the adopted depletion/set-pool representation) |

**Generator fixes identified:**
1. **Per-position uuid copies are joined only along precedingEvent.**
   - The WNT and TP53 cases.
   - `_get_or_create_entity_uuid`; reactions with no following event leave
     orphan copies.
2. **Output sets are dissolved to leaves,** orphaning the member complex that
   acts (CCNE1).
3. **`_leaf_uuid` treats depletion edges as reachability** (MDM2).
4. **The variant-cap bundle turns an OR set into AND inputs** (RAF;
   `_complex_variant_leafsets`, MAX_VARIANTS = 512).

**Scope / data (not topology):**
- Curated brakes from things absent in the assayed cell: drug-bound RAF, and
  a SARS-CoV-2 complex on PDPK1.
- CCC ATM: pin-form and scope.

**Solver:** the M4 loops in TP53, HDR and WNT are Reactome's own. No M4 loop
is a weld, so a loop-gain rule is legitimately a solver question there.
