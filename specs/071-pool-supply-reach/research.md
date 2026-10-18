# specs/071 research: results of arms psr and psrd

Solver `dc7f56a` (`DS_POOL_SUPPLY_REACH=material`), run 2026-10-18 with
`scripts/run_arm.sh`. Both logs print `Pool supply reach: material`.

| arm | build | results | regulated producers kept |
|---|---|---|---|
| psr | `20261017-1246_5e663d9` (canonical) | `results/dc7f56a/psr` | 0 |
| psrd | `20261018-0815_ddbd828_dsub` (specs/070) | `results/dc7f56a/psrd` | 108 |

## psr: the canonical build

The rule keeps no producer on the canonical build. No canonical pool has a
supply that it reaches only through a catalyst or regulator edge. Every
prediction is identical to `6fb4329` on both axes: held-out +0, tuning +0,
experimental +0, with and without RAF.

- Prediction 3 (no pathway beyond ±10): **holds**, trivially.
- The registered gate passes, but only because nothing moves. On these networks
  the rule is a no-op, so adopting it would change no number. The decision is
  Adam's.

## psrd: the faithful build

`holdout_report.py --set-aside RAF_MAP_kinase_cascade`:

| comparison | held-out | tuning | experimental |
|---|---|---|---|
| psrd vs dsub | +0 | **+110** (159 / 49, 1 pathway: TP53) | **+16** (21 / 5, p 0.0025; MDM2 54%) |
| psrd vs canonical | −171 (23 / 194) | −177 | −17 (22 / 39, p 0.04) |

All movement is in Transcriptional Regulation by TP53, the tuning set. Held-out
does not move, so the gap to canonical is unchanged: the DSB and RUNX1
mechanisms (specs/070 §5, #1 and #2) remain.

TP53 curator cases, dsub to psrd:

| | fixed | broken |
|---|---|---|
| MDM2 knockout | 49 | 0 |
| CHEK2 knockout | 47 | 0 |
| MDM4 overexpression / knockout | 44 / 1 | 0 |
| AKT1, AKT2, ATM, BRCA1 | 18 | 4 |
| CDKN2A knockout | 0 | 45 |

### Prediction 1 (≥ 80 of the 100 MDM2 cases return to the canonical class): **fails, 49 of 100**

The 49 that returned are all knockouts. None of the overexpressions returned.
MDM2 now reaches its readouts (98 of 102 predict a change, against 0 under
dsub). However, 46 overexpression cases read UP where canonical and the
curators read DOWN. That class is the same at iteration caps 500 and 2,000.

Trace (MDM2 overexpression pinned at 80, cap 2,000;
`~/deltasignal-catalogs/trace_071/`):

- The MDM2 pool (pool1) rails to 100×. The TP53 ubiquitination pool (pool2:
  TP53 tetramer → MDM2:TP53 → PolyUb-TP53 → USP7 → TP53 tetramer) goes the
  wrong way: PolyUb-TP53 reads 0.29× and free TP53 tetramer reads 1.63×.
- The catalyst of the return step, free USP7 (R-HSA-3215278), reads 100×. USP7
  is a pool4 carrier: it is bound into DAXX:(PolyUb,p-MDM2):USP7 and released
  by "USP7 deubiquitinates MDM2" (R-HSA-3215295). The carrier rule
  (`DS_CYCLE_CARRIERS=on`) drops the release steps listed in
  `pool_carriers.csv`. That table lists only one variant copy of the release
  (`8e405807`).
- A second variant copy (`33469e3e`) reads a variant of the DAXX complex
  (`bb11841f`) that is not a pool4 member. It reads 100× and it feeds USP7
  unfiltered. So free USP7 follows the flux through MDM2 deubiquitination, and
  MDM2 overexpression deubiquitinates TP53.
- A smaller leak: the MDM2:TP53 intermediates read flux (11× and 47×), and
  p14ARF sequestration (R-HSA-6804998 → R-HSA-6804996) releases free TP53
  tetramer from them at the same fold.

This is the open-carrier-loop class of specs/070 §5 #1 (a release read as
supply), reached here through a variant copy that the pool tables miss. It is
not a defect of the 071 rule. The 071 rule makes the path reachable.

### Prediction 2 (TP53 non-convergence not below 663): **holds, 765**

TP53 non-convergence rises from 663 to 765 of 1,836 valid cases. Every MDM2,
MDM4, CHEK2 and CDKN2A case is non-converged. The LM fallback skips the
component because it holds pool nodes (`nonconv_skipped_pools=1`).

The iteration cap decides some of these classes (specs/066's test):

- **MDM2 knockout**: converges at cap 2,000, with the same class as at 500.
- **MDM2 overexpression**: does not converge at either cap; the class is
  stable.
- **CDKN2A knockout**: 47 of 51 readouts change class between cap 500 (UP) and
  cap 2,000 (DOWN). Curators expect no change, so the case is wrong either way.
  These are the 45 cases psrd breaks. CDKN2A scores 102 / 49 / 4 correct on
  canonical / dsub / psrd.

## Reading

- On the faithful networks the rule does what it was built for: MDM2
  translation is a supply again, and MDM2 knockouts read correctly.
- It also exposes the next defect on the same path: a USP7 release copy that
  the pool tables miss.
- The TP53 tuning gain (+110) is not trustworthy as a number. The component
  does not converge, and at least the CDKN2A cases are decided by the cap.
- Next, in order (each to be registered before it is measured):
  1. Make the carrier rule cover every variant copy of a listed release
     reaction (match by `reaction_stid`, not uuid). This could be a solver fix
     or a generator fix.
  2. LM fallback for pool components (specs/066 skips them).
  3. specs/070 §5 #1, the DSB open carrier loops, which psrd did not touch.
