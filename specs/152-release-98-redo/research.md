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
