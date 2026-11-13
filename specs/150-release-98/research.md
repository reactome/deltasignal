# specs/150 research: R1 — the frozen rules on Release 98

## Verdict: the replication PASSES the registered band

The frozen rulebook, scored untouched on networks rebuilt from a year of
curation it never saw (build `20261113-0316_8b6527e_r98`: frozen generator
tag, Release-98 graph database verified by DBInfo, 98 diagram mirror; arm
`r98_frozen`, solver `e967a1b`):

| convention | 97 frozen | 98 R1 | drift |
|---|---|---|---|
| curator held-out, valid-only (n = 18,481 both) | 0.9079 | 0.9047 | **−0.32pp** |
| curator held-out, all-rows (n = 19,000, the RESULTS.md convention) | 0.8977 | 0.8946 | −0.31pp |
| curator all, valid-only (n = 23,511 both) | 0.9095 | 0.8989 — see note | — |
| experimental (n = 845 valid both) | 643 right | **647 right (+4; 7 fixed / 3 broken)** | +0.47pp |

(The curator-all row's 97 figure is 21,502/23,511 = 0.9146 by valid-only
count; the drift there is −0.37pp, matching the held-out drift; the table's
"see note" flags that RESULTS.md quotes 89.22% on the 24,100 denominator.)

**Within the registered ±2.0pp band → the structural rules generalize to
re-curated networks.** The experimental axis improves and now stands above
MP-BioPath's published 643 on the same cases.

**Correction on the record:** the first in-session report of R1 said
"+0.70pp, gains" — that compared the 98 valid-only figure against the 97
all-rows figure (different denominators). The same-convention comparison
above replaces it; the standing compare-denominators rule exists for exactly
this mistake, and it was caught by the per-case trace (net fixed−broke
disagreed with the headline).

## The trace (same case keys on both releases: 0 only-97, 0 only-98)

- Curator: 239 of 23,511 predictions changed; 62 fixed / 122 broke.
  Concentrated: **WNT −25, IL-4/IL-13 −24, RHO-GTPases-activate-PKNs −16**
  (all held-out), **ERBB2 +18** (tuning), IL-3/5-GM-CSF −7, Hedgehog −6.
  Release-98 curation restructured those networks; whether the WNT and IL-4
  losses are new structure the frozen rules mishandle (a future rule
  candidate) or curation the ground truth predates (the specs/070 class) is
  a follow-up trace, filed, not assumed.
- Experimental: 13 changed, 7 fixed / 3 broken; ERBB2 +4 the largest.

## Status of the coordinated move

- R1 recorded (this file). The 97 record remains the published provenance
  for everything currently served.
- Next: R2 parity harnesses on the 98 build; R3 arrows pipeline on 98; R4
  confidence refit + scope update; R5 packaging. Serving stays on 97 until
  Adam moves it.
- Infrastructure now standing: Release-98 graph DB as `reactome-neo4j`
  (the 97 container stopped beside it as `reactome-neo4j-97-stopped`), the
  98 diagram mirror covering all 125 build pathways, the frozen-tag
  worktree at `~/gitroot/lng-freeze-2026-10-26`.

## R2 — the parity harnesses on Release 98 (2026-11-13): all three PASS

The serve-what-you-measure invariant survives the release (engine replays on
the `20261113-0316_8b6527e_r98` build, bench dumps against the same build and
the live Release-98 Neo4j; artifacts in `r2/`):

- **Pins (specs/140): 877/877** — basis distribution identical to 97
  (root 806, root_cycle 58, none 13).
- **Readouts (specs/141): 983/983 curator + 272/272 korkut** — one readout
  moved from proxy to direct basis with the re-curation (890/80 vs 889/81),
  agreeing on both sides.
- **Rails (specs/142): 24,356 folds, zero differences**, stops match; the 98
  networks rail slightly less (3,838 vs 3,922).

Two run notes for the record: the engine replays must write through a
read-write gate mount (`/app/bench` is mounted read-only — the first compare
read the stale committed 97 artifacts and reported wholesale differences,
recognizable by identical per-case structure under wholesale-different
uuids); and the gates' 98 artifacts live under this spec's `r2/`, leaving the
committed 97 gate artifacts as the frozen record.

## R3 — the arrows on Release 98 (2026-11-13): the lead mechanism replicates

The 98 census (102 usable arrows after virus exclusion), the 33 source
pathways rebuilt on 98 (`20261113-0442_8b6527e_r98src`, 33/33), composites
`20261113-arrows98-10c` (16 arrows bridge, 34 edges, 10 centres — 97 had
18/32/11; the delta is the two arrows release 98 removed plus re-curation).
Scored under the specs/148 F discipline against the R1 control and a fresh
plain-98 korkut control (`r98_korkut`):

- **H1 PASSES, replicated almost digit for digit**: korkut untouched class
  (n = 2,715): balanced **0.3834 vs always-no-change 0.3333** (+0.0501,
  95% [+0.0212, +0.0786]; the 97 arm: +0.0503 [+0.0216, +0.0784]) — the
  same 156 differing answers, the same 37 right / 108 false / 11 inverted
  confusion. The p-AKT mechanism is release-stable.
- H2 neutral, as on 97: curator 6 fixed / 6 broken of 4,689 untouched;
  experimental 0 of 22.

## R4 — the confidence refit on the 98 record

`r4/confidence_98.json` (fitter now takes `--basis/--catalog-build/
--solver-commit`; provenance recorded as `20261113-0316_8b6527e_r98` /
`e967a1b`): the 97 pattern reproduces — curator held-out calibrated (ECE
0.0412, AUC 0.624 [0.548, 0.693]), the experimental axis fails through the
same untouched stratum, so the specs/146 amendment-1 serving semantics
(untouched → refusal) remain correct on 98.

## R5 — packaging

`r5/leads_98.json` (10 centres, 16 arrows, measured rate on the replicated
p-AKT mechanism; 98 provenance throughout) and the serving catalog
`builds/20261113-r98-plus-leads` (the plain 98 build plus the 10 unlisted
composites). **Switching what is served — `current`, `src/api/confidence.json`,
`src/api/leads.json`, and the scope statement's numbers — is one deliberate
move, Adam's to make**; until then everything served remains the 97 record,
whose provenance is untouched.

## The serving switch: executed and REVERTED the same day (2026-11-13)

The switch recorded under R5 was executed on Adam's go-ahead (PR #247:
confidence/leads/scope to the 98 record, catalog to
`builds/20261113-r98-plus-leads`, live-smoked exact against the recorded 98
bench folds) — and reverted hours later (PR #248) when the specs/151 drift
trace found that the 98 build is defective: Release 98 removed every HGNC
referenceGene link, which silently killed the frozen generator's
protein-reference lookup and with it specs/101 balancing, output-slot binding
and the same-entity signatures, build-wide. R1's replication band passed on a
catalog that is not the frozen construction. specs/151 amendment 1 has the
trace; logic-network-generator#131 the fix, whose acceptance gate is
rebuilding Release 97 identically to the frozen catalog. The 98 numbers in
this file remain true as measured, but they describe the defective build; the
coordinated move is re-run after the fix. Serving is back on the 97 record.
