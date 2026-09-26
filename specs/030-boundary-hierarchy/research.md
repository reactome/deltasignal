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
