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

## The rule (logic-network-generator `LNG_BOUNDARY_HIERARCHY=1`, PR #97; off by default)

Decompose a root complex one hasComponent level at a time:
- a component that already exists as a node, and is not downstream of the
  root (the specs/018 rule), is **joined** to the complex and not descended
  into;
- a nested complex that does not exist is **built** as a node and decomposed
  in turn;
- anything else falls back to its terminal leaves, as before.

Three generator tests: default flat behaviour unchanged; nested complex built
and produced species joined; a copy downstream of the root not reused.

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
