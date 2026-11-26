# specs/174: why each Lun 2019 half-A case passes or fails, and whether K562's curated regulators act

**Status:** registered 2026-11-26, before any per-case class was counted.
**Requested by:** Adam, 2026-11-26 ("do you think it would be helpful to look through and figure out why each one of
the failures failed and why each one that passed passed", then "ok"). Analysis only: no rule, no adoption.

## Why

specs/173 found that on Lun 2019 half A a route carries some signal (enrichment 1.44 [1.13, 1.82]) but the frozen model
is no better than reachability (balanced accuracy 0.379 against 0.372). specs/169 found no route signal at all on K562
CRISPRi. Before building anything that learns which routes transmit (a candidate specs/175), we need to know, class by
class, why passes pass and failures fail: background responsiveness, a missing or wrong curated edge, an edge silent in
this cell, wrong direction, or the readout itself.

## Data and roles

- **Lun 2019 half A** (`lun2019_a`, FIT): the 3,127 measured P ∩ A pairs scored by specs/173 (code `bad27ef0`, same
  measured call, same frozen top114 predictions, sha256 `0c0ec033…`). **Half B is never read, summarised or counted.**
- **Replogle 2022 K562 split A** (`replogle2022_a`, SELECT): specs/169's truth and predictions. Split B stays sealed.
- Per-case tables are written outside the repo (`~/deltasignal-catalogs/analysis/174/`). The repo gets counts and
  classes, as for specs/173.

## Registered classes (Lun half A)

**Readout type**, from the P ∩ A cases of each readout (pairs with a route: `reach` UP):
- **background**: changes in ≥ 25% of the pairs without a route;
- **route-specific**: not background, at least 10 pairs with a route, changes in ≥ 20% of them and at least twice as
  often as without a route;
- **quiet**: the rest.

**Case class** (`frozen` and `measured` in {up, down, not}; a conflict counts as not):

| class | route | model | measured |
|---|---|---|---|
| R1 right change | yes | up/down | same direction |
| R2 wrong direction | yes | up/down | opposite direction |
| R3 called, silent | yes | up/down | not |
| R4 missed on route | yes | not | up/down |
| R5 right no-change | yes | not | not |
| N1 changed, no route | no | — | up/down |
| N2 no change, no route | no | — | not |

Each class is crossed with readout type. A pass is R1, R5 or N2; a failure is R2, R3, R4 or N1.

**Registered questions** (descriptive; no gate):
1. What share of R1 sits on background readouts? If more than half, most of the model's right change calls are
   readouts that move under any over-expressed kinase, and they say little about routes.
2. How concentrated are the failures? The share of R2+R3 in the five readouts with the most of them.
3. For each failure class, the dominant readout × perturbed-protein groups (largest first, until 50% of the class is
   covered, at most eight groups), each traced through the network: the shortest signed path from the pinned protein
   to the readout in the scored pathway, and the model's fold at the readout. Each traced group gets one mechanism
   label: **missing negative** (biology damps the readout, the network has no such edge), **inactive in cell** (the
   route is curated but the readout does not move), **sign** (the network's sign is opposite), **readout**
   (the measurement cannot show the predicted change, e.g. a total-protein readout read as a modified form), **other**
   (stated). The label is a hypothesis from the trace and the literature, recorded as such.
4. Passes: for R1 on route-specific readouts, the same tracing for the largest groups, to state which curated routes
   transmit.

Readout nodes for tracing are matched from the catalog's own tables (`leaf_refs.csv`: UniProt accession and phosphosite
position), not Neo4j; a trace whose matched readout node differs from the scored one is reported as such.

## Registered check (K562 split A)

The pairs where the knocked-down gene's protein is in a regulator (positive or negative) of a reaction producing the
target's protein, in the served Release 98 networks (`logic_network.csv` regulator edges, `leaf_refs.csv` accessions).
For each: measured change, knockdown efficiency, and the knocked-down gene's overall effect (share of its scored targets
that change). Question: do the regulator pairs that fail come from knockdowns that move nothing at all (the regulator is
not active, or the knockdown not effective, in K562), or from knockdowns that act elsewhere but not on the curated target?

## Not in scope

No rule, no fitting, no adoption. Anything proposed from this is registered separately and judged once on Lun half B.
