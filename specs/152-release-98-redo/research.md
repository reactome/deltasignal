# specs/152 research: the release-98 move, redone

## G1 — the 97 identity gate: PASSES (2026-11-13)

Build `20261113-1235_baf608b_hgncfix97`: the 92 catalog pathways on Release
97, generator at LNG PR #132's head `baf608b` (the reference-lookup fix),
through scripts/catalog.sh (seed recorded, 92/92 complete, and the build log
shows the machinery alive again: `balance_dropped`, `output_slots_bound`,
`slot_binding_ties` all present).

Against the frozen catalog `20261026-1349_8b6527e`, specs/065 method
(probe_065a1/content_cmp.py — nodes labelled by their full nodes.csv content,
edge multisets under that labelling): **zero differing pathways of 92**, and
zero row-count mismatches across every side-table CSV of every pathway. The
189 references that newly resolve under the fixed predicate (IG V segments,
UBE2O, SORCS3, CD21 — present in 7 of our pathways) change no signature
decision in our scope, so no Adam adjudication is needed: the fix is a no-op
on Release 97, exactly as gated.

Note on infrastructure: the 97 graph was served from a port-free container
(`reactome-neo4j-97-trace`, image Release97, needs
`-e NEO4J_dbms_jvm_additional=-XX:-UseContainerSupport`), reached by
container IP via `NEO4J_URL` — the original 97 container's port publishing
did not survive the release swap. `NEO4J_URL` is not an LNG_* flag, so it is
not in BUILD.json's env record; the DBInfo release check in BUILD.json is the
authority on which graph a build read.

LNG #132 itself is CI-green with the G1 result recorded on the PR; **merge
needs a human review approval (repo rule), waiting on Adam.** G2 builds from
the PR head `baf608b` — the exact commit G1 validated — which catalog.sh
records in BUILD.json.

## G2 — the 98 build, properly: PASSES (2026-11-13)

Build `20261113-1311_baf608b_r98fix` through scripts/catalog.sh: 92/92
complete; BUILD.json records the env (PYTHONHASHSEED=0, LNG_EMIT_ONE_SIDED=1,
LNG_DIAGRAM_DIR=~/reactome-diagrams/98) and DBInfo release 98 (2026-09-18,
checksum 525378806530). The machinery is alive: `balance_dropped` > 0 in 38
pathway logs, `output_slots_bound` present in 60.

The scale of the earlier defect, in one number: this sound 98 build is
**209,876 nodes / 493,651 edges** — +0.2% over the frozen 97 catalog
(209,480 / 492,034), which is real curation growth. The defective build
(`20261113-0316_8b6527e_r98`) had 230,396 nodes and 332,885 edges — nodes
+10% from cross-product enumeration while edges fell 32%, because the dead
identity signatures also changed how copies share edges. Both distortions are
gone in the sound build.

## G3 — the replication, re-run on the sound build: PASSES EXACTLY (2026-11-13)

Arm `r98fix_frozen` (frozen solver, results under
`builds/20261113-1311_baf608b_r98fix/results/12b5db3/`):

| convention | 97 frozen | defective 98 (specs/150 R1) | sound 98 |
|---|---|---|---|
| curator held-out, all-rows (n=19,000) | 0.8977 / F1 0.8687 | 0.8946 | **0.8977 / F1 0.8687 — identical** |
| curator valid-only, 23,511 shared keys | 0.9014 | 0.8989 | **0.9014 — zero class changes** |
| experimental (scored 849) | 643 | 647 | **643 — identical** |

Zero predictions change on any axis: the per-pathway drift table is empty.
The catalog itself genuinely differs from 97 in 12 of 92 pathways (small
curation edits: NOTCH1 −1 node, RUNX1 +6, antibody pathways reshaped slightly,
+2-node additions in S_Phase, Mitotic G1, Cell Cycle Checkpoints, senescence,
Cell junction organization), and none of it touches a scored route.

Two conclusions, both owed to the trace:
- **The entire specs/150 drift was the defect** — the −0.32pp, the
  WNT/IL-4/13/PKN concentration, and also the +4 experimental "improvement"
  (647, "above MP-BioPath") that the defective cross-products manufactured.
  On a sound build, Release 98 neither helps nor hurts the benchmarks: the
  frozen record carries over digit for digit.
- The replication claim is now the strongest it can be: the rulebook's
  answers are invariant under a year of re-curation of these 92 pathways.

## G4 — parity: PASSES, zero differences (2026-11-13)

On the sound build (`specs/152-release-98-redo/g4/`, bench dumps via the
specs/140/141/142 harnesses, engine replays in the dev container through a RW
gate mount): pins 877/877 (root 806, root_cycle 58, none 13), readouts 983
curator + 272 korkut with zero differences, rails 24,356 folds with zero
differences and 3,922 railed — the frozen record's exact rail count.

## G5 (part 1) — the confidence refit: identical to the shipped table

`g5/confidence_98.json` (fit on the G3 results, basis from the g4 readout
dump): validation reproduces the frozen table digit for digit (held-out
0.9079 / ECE 0.0499 / AUC 0.6075 [0.528, 0.679]; experimental 0.7609 / ECE
0.1491), and the table CONTENT is bit-identical to `src/api/confidence.json`
(the 97 table) apart from the provenance fields — as it must be, given G3's
zero prediction changes. Provenance: `20261113-1311_baf608b_r98fix` /
`12b5db3`.

## G5 (part 2) — arrows, leads, packaging: COMPLETE (2026-11-13)

- **H1 PASSES identically for the third time**: korkut untouched class
  (n = 2,715), arm `arrows98fix_w`/`arrows98fix_k` on composites
  `20261113-arrows98fix-10c` (16 arrows, 34 bridges, 10 centres — matching
  R3): balanced 0.3834 vs 0.3333 (+0.0501, 95% [+0.0212, +0.0786]), the same
  156 differing answers, the same 37 right / 108 false / 11 inverted
  confusion (down-call precision 0.24). The p-AKT mechanism is stable across
  release 97, the defective 98 build and the sound 98 build. H2 neutral
  (curator 6 fixed / 6 broken of 4,687 untouched; experimental 0 of 20).
- `g5/leads_98.json`: the 10-centre/16-arrow table with sound provenance.
- Serving catalog `builds/20261113-r98fix-plus-leads` packaged (specs/149
  pattern: the sound 92-build + ctx-arrows_ composites, BUILD.json noting
  the composite build; its build_id is the measured
  `20261113-1311_baf608b_r98fix`, so the served scope comparison holds).

## The serving switch, re-executed (2026-11-13)

Adam's go-ahead stands from this spec's Why. The switch:
`src/api/confidence.json` and `leads.json` <- the g5 tables;
`src/api/scope.json` -> provenance `20261113-1311_baf608b_r98fix` /
`12b5db3` with the FROZEN numbers (89.77% / macro-F1 0.8687 / 19,000;
643/849 = 75.74% / 0.6638), because the sound 98 build reproduces them
exactly — the statement now says so in words. test_readout_protocol
provenance assertions and the CLAUDE.md scope-suite count (251 on the
10-centre table) move with it. Serving: PATHWAY_CATALOG ->
`builds/20261113-r98fix-plus-leads`.

## Amendment 1 (2026-11-20): G5 re-run on composites with every source pathway built

specs/149 amendment 1 (code review) found that `20261113-arrows98fix-10c`, the G5 composites, wrote the bridges of
the 11 arrows whose source pathway is outside the 92 but not the source pathways themselves, so 21 of 34 bridges
dangled. Adam chose the rebuild (2026-11-20). Rebuilt with the fixed builder, the same inputs and
`--sources 20261113-1345_baf608b_r98src2`: `builds/20261120-arrows98fix-10c-src`, **0 of 34 bridges dangling**,
counted independently. G5 was then re-run exactly as registered: `run_arm.sh` at `206e3bbb`, arms `arrows98src_k`
(korkut) and `arrows98src_w` (curator + experimental, `--limit 10`), scored with
`specs/148-diagram-arrow-bridges/build/score_148.py` against controls `r98fix_korkut_h` and `r98fix_frozen_h`,
re-run at the same commit on `20261113-1311_baf608b_r98fix`. The controls are identical to the recorded
`r98fix_korkut` (cff961d) and `r98fix_frozen` (12b5db3): 0 fold differences on 11,626 korkut, 24,100 curator and 849
experimental rows.

| | recorded (dangling composites) | rebuilt composites |
|---|---|---|
| H1 korkut untouched, n | 2,715 | 2,715 |
| balanced, arm vs always-no-change | 0.3834 vs 0.3333 | **0.3834 vs 0.3333** |
| difference, 95% by-drug-set | +0.0501 [+0.0212, +0.0786] | **+0.0501 [+0.0212, +0.0786]: PASSES** |
| differing answers; right / false / inverted downs | 156; 37 / 108 / 11 | 156; 37 / 108 / 11 |
| H2 curator untouched (4,687): fixed / broken | 6 / 6 | **6 / 11** |
| H2 experimental untouched (20) | 0 changed | 0 changed |

- **H1 is unchanged to the digit.** The 11 arrows that are now built move no korkut untouched answer, fold or call,
  in any centre: PIP3 0 of 786, RUNX2 0 of 360, Cell junctions 0 of 281, PTK6 0 of 195, ROBO 0 of 56. Mitotic G1's
  75 moved answers are the same p-AKT ones. specs/148's addendum claim ("the 12 new arrows change no korkut answer")
  is now true as a measurement, though the reason it gave was not one.
- The bridges are live. In PIP3 the PTEN bridge changes 250 of the 273 scored korkut readouts PIP3 already routes.
  F keeps PIP3's own answer for all of them, so H1 does not see them.
- **H2 gains five breaks, all on the new arrows.** RB1 overexpression in Mitotic G1 reads DOWN at three readouts
  (0.73–0.83×) through the senescence → INK4 arrows. RB1 knockout moves six more readouts slightly (1.01–1.12×, still
  NORMAL). RAC1 in ROBO reads DOWN under knockout (0.75×) and UP under overexpression (20.75×) through PAK. The
  curators say NORMAL in every case. H2 was registered "net-positive on each axis", so it now fails negatively on
  curator (−5), where before it failed neutrally. The curator set is SELECT. The client serves such answers
  only as labelled leads.
- Per-arrow outcomes and the served statements: specs/149 amendment 2.
