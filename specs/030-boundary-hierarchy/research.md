# 030 — Build the nested complexes the generator never builds (step 1 of Adam's plan)

## Adam's plan, 2026-09-26

*"fix the 252 where it never builds a complex. And then when that is fixed
figure out how to fix when the protein has no edge to the complex in the 138
cases. And then try to fix the missing edges to reaction for the 12 cases.
once those are built try to work with the deltasignal to make the pinning
better."*

## What the 252 actually are (traced before any code)

**Traced: JAK1 KD in Interferon α/β.** JAK1 is pinned and flows through the real
reactions (receptor binding → IFNAR1 recruitment → JAK activation → STAT
phosphorylation → release of p-STAT2:p-STAT1) to **ISGF3 [cytosol]**. That is
produced, and consumed by nothing. The next complex, **ISGF3:KPNA1:KPNB1**, is a
root: produced by nothing, and translocated to the nucleus. Reactome has no
reaction binding ISGF3 to KPNA1 and KPNB1; Reactome joins them only by
hasComponent, through the nested **ISGF3:KPNA1**, which no reaction uses.

**How the generator decomposes root complexes.** It goes straight to base leaf
proteins (`get_terminal_components`: STAT1, STAT2, IRF9, KPNA1, KPNB1). It skips
every intermediate complex, so the ISGF3 node the pathway produces is never
joined to the root that contains it.

**Across all 438 severed routes.** A route that prefers reactions and uses
composition only where it must shows:
- 178 need **no** composition step at all: Reactome connects them by reactions
  alone, so the break is some other generator defect (next steps);
- the composition steps needed go mostly into ROOT complexes: 264 of 364.

## The rule (logic-network-generator `LNG_BOUNDARY_HIERARCHY`, PR #97; default ON since ba624ce, `=0` for the flat mode)

Decompose a root complex one hasComponent level at a time:
- a component that already exists as a node, and is not downstream of the
  root (the specs/018 rule), is **joined** to the complex and not descended
  into;
- a nested complex that does not exist is **built** as a node and decomposed
  in turn;
- anything else falls back to its terminal leaves, as before.

Generator tests (tests/test_boundary_hierarchy.py):
- flat mode (`=0`) unchanged;
- the nested complex is built and the produced species joined;
- a copy downstream of the root is not reused;
- the hierarchy is the default;
- after the review, below: root-copy preference, two joins cannot jointly
  close a cycle, and a nested complex is built per root.

## Pre-registration (committed before any build or arm)

**Builds** at generator `ae8016c`: `hier` (`--env LNG_BOUNDARY_HIERARCHY=1`)
and `hierctrl` (no flag), so both share one generator commit.

**Arms:** `hierctrl` and `hier`, current defaults, through
`scripts/run_arm.sh --catalog`.

**Decides:** curator held-out in-release net, McNemar p, pathways moved,
perturbations. **Also reported:** net outside IFN α/β; the experimental axis;
how many of the 438 severed routes become correct; and the change in
strongly connected component sizes, since new assembly edges must not weld
loops.

**Adopt** only if held-out net > +15 with p < 0.05, positive outside IFN α/β,
and experimental no worse than −15.

## Result: ADOPTED

Builds `20260926-0834_ae8016c_hierctrl` and `20260926-0835_ae8016c_hier`,
bench `e6bc117`, current defaults.

| | held-out net | p | notes |
|---|---|---|---|
| **hier vs hierctrl** | **+228** (238 / 10) | 9e-58 | 8 pathways moved (+5 / −2), 38 perturbations |
| Interferon α/β | +200 | | every one of its 200 severed held-out cases |
| **outside IFN α/β** | **+28** (38 / 10) | 6e-5 | DSB Repair +27 |
| tuning | +6 | | TP53 |
| experimental | +2 (2 / 0) | | |

- **Severed routes.** 216 of the 438 severed held-out routes are now correct,
  against 0 in the control.
- **Convergence** is unchanged: 1,673 against 1,671 of 1,725 solves.
- **Structure.** The hier build has 948 more nodes (the nested complexes it
  builds) and 630 fewer edges: joining an existing component replaces several
  leaf edges with one.
- **Cycles.** Cyclic nodes rise by 157 of 8,076 (+2%) in 5 pathways: DNA
  Repair +85, DSB +53, HDR +15, RAF/MAP +6, and one −2. That is the known
  residual of the generator's reachability snapshot, which is not updated as
  assembly edges are added. DSB still gains +27.

It meets every pre-registered condition. Next: make
`LNG_BOUNDARY_HIERARCHY=1` the generator default (LNG PR #97), rebuild the
canonical catalog, and re-baseline; then Adam's step 2.


## Review of PR #97: the gain outside IFN α/β was mostly a weld

The adversarial review reproduced every number above, then showed how much
less the part outside IFN α/β means than it looked:

- **Two joins welded a new loop.** Each passed the specs/018 downstream test
  against the pre-loop snapshot, and neither saw the other.
  - p-RPA heterotrimer joins a nested complex under LIG1:ERCC1:…, which
    catalyses Completion of SSA.
  - RAD9:HUS1:RAD1, an output of Completion of SSA, joins a root that
    produces p-RPA.
  - Together they put readout 5686663 inside a new 162-node SCC. 16 of DSB's
    29 fixes read out there, and so do both of its breaks. **The DSB "+27" is
    not independent of the weld.** The cycle rise, +157 nodes, is this
    mechanism, not the 4-edge specs/018 residual.
- **What remains outside IFN α/β.** Excluding DSB, the other 6 pathways net
  **+1** (9 / 8).
- **"216 of 438 severed routes"** is IFN α/β's 200 plus DSB's 16 through the
  weld. **0 of 222** were repaired in the other 17 pathways.
- **Concentration.** IFN α/β's +200 comes from 8 perturbations (JAK1, IFNAR2,
  PTPN11, SOCS1, both directions) over 25 readouts. It is clean: IFN gains no
  cyclic nodes. The McNemar p-values assume independent cases.
- **Nested complexes shared across roots.** A shared nested complex was
  decomposed against whichever root came first: no welds, but 4 missed joins
  in S Phase, and uuid-order dependent.
- **Only the first eligible copy is joined.** That is also Adam's step 2:
  PDGFB's pinned ROOT copy fed only the processing reaction, while the PDGF
  A/B heterodimer's assembly edges came from a produced PDGFB copy.

**Fixed (LNG 8563de9, hierarchy mode only; flat mode byte-identical):**
- the downstream test is updated as each join is emitted;
- root complexes are processed in a deterministic order (stid, then uuid);
- nested complexes are built per root;
- when joining, a ROOT copy (in-degree 0) is preferred over a produced copy.

Tests: two joins cannot jointly close a cycle; root-copy preference;
per-root nesting. Both of the first two are mutation-checked.

## Pre-registration: re-measure (committed before the arm runs)

- **Setup:** the `hier2` build at LNG 8563de9 (hierarchy default on), against
  the same `hierctrl` build and arm (flat mode unchanged).
- **Adopt** only if held-out net > +15 with p < 0.05, the net **outside IFN α/β
  is positive**, and experimental is no worse than −15.
- **Also reported:** the cyclic-node change (it must not rise because of these
  joins); IFN α/β; DSB; the other pathways; and severed routes repaired
  outside IFN α/β.
