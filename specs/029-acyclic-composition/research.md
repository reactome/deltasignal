# 029 — Composition edges that bridge severed routes without closing cycles

## Why

In the corrected failure anatomy (specs/026 corrections), 438 held-out misses
have a route in Reactome that our network severs. **390 of them break at a
component → complex step**: 208 where the complex node exists but the
component has no edge into it, and 182 where the complex is a nested
sub-complex the generator never creates.

`LNG_COMPOSITION_EDGES=1` adds component → containing-complex edges; as a class it loses:
−358 to −1,039 held-out outside IFN α/β. The traced cause is that
composition edges re-close the loops that the specs/018 fix opened.

**The rule (`DS_COMPOSITION_FILTER=acyclic`, bench side).** On the composition
build, drop each composition edge that would close a cycle: added one at a
time in sorted order, each checked against the network plus the edges kept so
far. That drops **1,114 of 5,651**; 4,537 remain. Tested:
- a bridging edge is kept and a closing edge is dropped;
- two edges that close a cycle only together are caught.

## Pre-registration (committed before any arm runs)

**Builds:** `20260926-0124_d4f4f64_ctrl` and `20260926-0126_d4f4f64_comp`
(specs/026). All arms are at the same commit, with current defaults
(`DS_KO_AGG=mean`).

| arm | build | settings |
|---|---|---|
| `ctrl` | ctrl | — |
| `comp_acyclic` | comp | `--bench DS_COMPOSITION_FILTER=acyclic` |
| `comp_acyclic_limit_novel` | comp | `--bench DS_COMPOSITION_FILTER=acyclic --server DS_COMPOSITION_MODE=limit_novel` |

- **Decides:** curator held-out in-release net against `ctrl`, McNemar p,
  pathways moved, perturbations.
- **Also reported:** net outside IFN α/β; the experimental axis; how many of
  the 438 severed routes become correct.
- **Adopt** only if held-out net > +15 with p < 0.05, **positive outside IFN
  α/β**, and experimental no worse than −15.

## Result: NOT adopted, and cycles are not the explanation

Arms at `f749b30`, paired against `ctrl` on the ctrl build:

| arm | held-out net | p | outside IFN α/β | IFN α/β | tuning | experimental | severed routes now correct |
|---|---|---|---|---|---|---|---|
| `comp_acyclic` | −666 (406 / 1,072) | 2e-69 | −810 | +144 | −72 | −18 | 156 of 438 |
| `comp_acyclic_limit_novel` | −11 (120 / 131) | 0.53 | **−111** (20 / 131) | +100 | +8 | −1 | 106 of 438 |

- **Filtering out cycle-closing edges helps, but not enough.** Compared with the
  unfiltered specs/026 arms, the loss outside IFN α/β shrinks: −1,039 to −810
  in the default mode, and −358 to −111 with `limit_novel`. It stays well
  below zero, so the rule fails.
- **Traced: DSB Repair still loses 108 with no cycle closed.** Of its 114
  breaks under `limit_novel`, 105 are false changes (NORMAL → DOWN), from
  knockdowns of PARP1, PARP2, FEN1, XRCC1, LIG3 and POLQ (15 readouts each).
  These are the MMEJ / alt-NHEJ proteins. Composition edges carry each
  knockdown into containing complexes that the other DSB readouts depend on,
  so readouts the curators call unchanged go DOWN.
- **This corrects the specs/026 and 018 inference.** The composition loss
  outside IFN is **over-coupling through complex membership**, not only
  re-welded loops. It is the same class as the specs/023 over-coupling finding,
  where false change is the largest error type.
- **Where the bridge stands.** A composition edge bridges a real severed
  route, as in IFN α/β, but a component's fold should not automatically
  propagate to every complex that happens to contain it. What distinguishes
  the two cases is not yet identified.


## Corrections (third review of PR #74)

- **The premise does not hold outside IFN α/β.** Composition edges route all
  200 of IFN α/β's severed held-out routes, but only **46 of the 238 outside
  it**, and 12 once the acyclic filter is applied. So outside IFN α/β, −111 is
  almost pure cost, not a real bridge outweighed by coupling. Of "156 / 106 of
  438 now correct", 150 / 100 are in IFN α/β and 6 are outside it.
- **The shrink compared with specs/026 is confounded.** The specs/026 arms
  used `DS_KO_AGG=max` and these arms use `mean`; the two controls alone
  differ by held-out +65. How much of the smaller loss is the filter, and how
  much is the readout rule, is not separable without an unfiltered arm under
  `mean`.
- **"Cycles are not the explanation" was an overclaim.** No new cycle is
  closed, but the damage flows through existing loops. PARP1's knockdown
  reaches 1,606 nodes instead of 35 by entering strongly connected components
  of 236 and 109 nodes, and 11 of its 13 broken readouts are reachable only
  through them. Whether it is `limit`'s one-way lowering or a loop collapse is
  not traced to a line.
- **Concentration.** DSB's −108 comes from 6 knockdowns across the same 15
  readouts, so McNemar is not meaningful there. The IFN α/β gain comes from 4
  genes and 25 readouts.
- **Convergence:** ctrl 1,649 of 1,725 solves converged; comp_acyclic 1,691 of
  1,725.
- **Order dependence.** 80 of the 1,114 dropped edges are chosen only by
  sorted uuid order, so which ones drop is label-dependent in the specs/013
  sense.
