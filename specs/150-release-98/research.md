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
