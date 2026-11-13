# specs/152: the release-98 move, redone on a sound build

## Why

specs/151 found the specs/150 Release-98 catalog defective: Reactome 98
removed every HGNC referenceGene row, the frozen generator's reference lookup
returned None for every protein, and balancing / output-slot binding / the
identity signatures were silently off build-wide. The 98 build also ran
without PYTHONHASHSEED=0 and without LNG_DIAGRAM_DIR (so it read the 97
diagram mirror), because it was invoked directly rather than through
scripts/catalog.sh. The serving switch was reverted the same day (PR #248).

The generator fix is logic-network-generator PR #132 (issue #131): both
reference queries now select a human ReferenceGeneProduct with a
referenceGene (same ids as the HGNC chain on 97; RNAs and non-human proteins
unchanged; 189 gene products whose HGNC rows were missing from release slices
now resolve, touching 7 of our 125 pathways), plus a loud refusal when a
build resolves zero references over >= 10 accessioned entities.

Adam (2026-11-13): "yes make a fix in the generator to make it so that we can
move to 98."

## Protocol (registered before any number is produced)

- **G1 — the 97 identity gate (the fix's acceptance test).** Build the 92
  catalog pathways on Release 97 with the fixed generator through
  scripts/catalog.sh (variant hgncfix97). Compare against the frozen catalog
  `20261026-1349_8b6527e` with the specs/065 content-labelling comparison
  (probe_065a1/content_cmp.py: nodes labelled by their nodes.csv content,
  edge multisets under that labelling). PASS = zero differing pathways.
  If only pathways among the 7 that touch the 189 newly-resolving references
  differ, each difference is traced to a named entity; adoption of such a
  difference is Adam's. Any difference outside those 7 fails the fix.
- **G2 — the 98 build, properly.** scripts/catalog.sh build on Release 98
  with `--env LNG_DIAGRAM_DIR=~/reactome-diagrams/98`, seed recorded in
  BUILD.json, DBInfo release checked = 98. The build must be complete (92/92)
  and the build log must show the balancer alive (`balance_dropped` > 0
  somewhere, `output_slots_bound` present).
- **G3 — the R1 replication, re-run.** Exactly specs/150 R1: frozen solver
  (tag freeze-2026-10-26 defaults) on the G2 catalog, both axes, same-
  convention comparisons against the frozen 97 record, the pre-registered
  band: curator held-out valid-only within ±2.0pp of 0.9079. The specs/151
  drift pathways (WNT, IL-4/13, RHO-PKNs) are re-read specifically: drift
  there should shrink if it was the defect, and what remains is re-curation.
- **G4 — parity.** specs/150 R2's three harnesses (pins, readouts, rails) on
  the G2 catalog: zero differences required, RW gate mounts as documented.
- **G5 — arrows, confidence, packaging.** specs/150 R3-R5 re-run: the 33
  source pathways built via catalog.sh on 98 (variant r98src2), the arrows
  pipeline and H1 gate, the confidence refit, leads table, packaged serving
  catalog. The serving switch is then re-executed (Adam's standing go-ahead,
  this spec's Why) with the specs/150 switch checklist: scope numbers under
  the verified conventions, provenance tests updated, live smoke against the
  recorded folds.

## Non-goals

No solver change; no rule changes; the freeze stands. Only the generator's
reference lookup moves, and only through LNG #132's gate.
