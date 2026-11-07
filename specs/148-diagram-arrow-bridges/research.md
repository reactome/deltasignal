# specs/148 research: the drawn arrows, measured

## The census (the build's first finding)

Of 115 release-97 `pathway→entity` arrows: 53 name source pathways outside the
92-pathway catalog (42 of those would resolve only through whole top-level
networks — the whole-Reactome sizing, SC-002); of the 62 in-catalog arrows,
**35 point at entities the networks already produce internally** (boundary
expansion and encapsulated sub-pathways — the generator absorbs most of what
the arrows draw; the specs/147 verdict again), 15 have no same-entity pair, 6
no produced source copy, and **6 bridge**: p-AKT from PIP3 into ERBB2 /
Mitotic G1 / Intrinsic Apoptosis; active ATM from DSB Response into Cell
Cycle Checkpoints; p-ABL1 and HES/HEY into RUNX1/2. 14 bridge edges, six
centres — small, but every one a famous handoff, including the single bridge
specs/132's tracing had validated as right (p-AKT).

## The measurement (arm `arrows97s`/`arrows97k` on build `20261107-arrows97-6c`, solver `ece7026`)

Scored under the registered F discipline (the pathway's own frozen answer
wherever it routes; the bridged answer only for control-untouched readouts):

- **H1 PASSES.** korkut, control-untouched readouts in the six centres
  (n = 1,553): bridged balanced accuracy **0.4056 vs always-no-change 0.3333**
  (+0.0723, 95% by-drug-set bootstrap [+0.0328, +0.1101]).
- **H2 fails as worded, neutrally.** Curator: 12 of 3,250 untouched cases
  change, 6 fixed / 6 broken, accuracy unchanged (0.8934). Experimental: none
  of the 31 untouched cases change. "Net-positive on each axis" is not met;
  no harm either.
- Run notes: `DS_CONTEXT_BRIDGE_ORDER=first` (amendment 1 — `last` trips the
  bench's composite-wide carrier-twin capacity guard, a census/solver
  mismatch, not a bridge finding) and the curator/experimental axes ran as a
  declared subset (`--limit 6`), since several full-run arm guards (specs/095
  copy pairs) legitimately fire on a six-pathway build.

## The traces

- **Every korkut fix is the drawn mechanism working**: all 37 fixed cases are
  AKT readouts (AKTpS473 15, AKTpT308 13, total AKT 9) under PI3K-inhibitor
  (ZS)-containing drug sets, in ERBB2 (19), Mitotic G1 (15) and Apoptosis (3)
  — the p-AKT arrow carrying a PI3K inhibitor's effect into pathways where
  AKT was an unreachable readout.
- **The cost rides the same bridge**: 108 false change calls (also mostly
  AKT-adjacent, ZS sets) — the PIP3 switch's overcalls (specs/136) imported
  with the signal. Raw accuracy on the untouched class falls while balanced
  rises: the bridge adds ranking information, not reliability. Measured
  reliability of a bridged "down" call here: 37 right of 156 made (24%).

## What this supports, and what it does not

The registered product shape stands up: bridged answers belong in a
**separate, labelled class** ("via a connected pathway — curator-drawn"),
with their own measured rates served the specs/146 way; on this evidence such
a call is a lead, not a conclusion, and the label must say so (a 24%-right
"down" is worth showing exactly because the alternative is a 70%-wrong
confident "no change" — but only with its number attached). It does NOT
support silently merging bridged answers into the main class — the same
verdict specs/133 reached, now with the mechanism resolved case-level.
Whether to serve the class is Adam's call; the whole-Reactome version (53
skipped arrows) needs top-level-network imports and is sized, not built.

## Provenance

Arrows: release-97 diagrams (the catalog's release), census committed;
member closures via `set_members_132.py`; build
`20261107-arrows97-6c` (context_rule `diagram_arrows`, whole sources, origin
table); arm runs under `run_arm.sh` (ARM.json beside each); control: frozen
`05dc867` (curator/experimental) and `a490cdc/k134_abund_a3` (korkut);
scoring `build/score_148.py`, seed 148. Registered at 6d9b02a before
measurement; amendment 1 at the guard failure, before any scoring.
