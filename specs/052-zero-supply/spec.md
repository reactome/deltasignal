# specs/052: zero supply — a protein that is not there is in no complex

**Status:** **not adopted** (2026-10-13; research.md R1). Curator held-out
−20 (p 1.1e-5), tuning −15, experimental +4 (n.s.). `DS_ZERO_SUPPLY_MODE`
stays `off`. Registered 2026-10-13, before any arm. The code and its unit
tests were written alongside this text. Nothing has been run on the
catalog.

**Origin:** specs/021, the decryptM trace (class 2), and Adam's go-ahead
(2026-10-13) to work on the largest traced defect.

## The defect

Take gefitinib on Signaling by EGFR, readout EGFR Y1092.
- Both EGFR root forms are pinned to 0.
- 413 of the network's 3,083 nodes still rail to 100×.
- The railed component is a cycle curated in Reactome (verified in Neo4j,
  release 97). EGF-like ligands:p-6Y-EGFR:CBL:GRB2 (R-HSA-182928) runs
  through 183058, 182994, 183002, 8951490 and 183084, which regenerates
  R-HSA-182928.
  - Every member contains EGFR.
  - Its co-inputs (CIN85:endophilin, Beta-Pix:CDC42:GTP, EPS15L1, EPN1)
    de-repress to the 10× depletion ceiling once the EGFR complexes that
    consume them fall away.
- The cycle then sustains itself. Complexes made of EGFR exist at up to
  100× in a cell with no EGFR.

That breaks conservation of mass, not a parameter choice. It is the
specs/014 positive-cycle knife-edge, set off by a knockout, and it is not
specific to EGFR.

Knockout cases predicted UP but expected otherwise:
- curator, canonical `1f58a1b`: 167 of 596 knockout→UP calls (WNT 33,
  RAF 31, MET 20);
- experimental: 44 of 75 (Mitotic G1 22, RAF 13).

How many of those are this mechanism is unknown. The arm measures it.

## The rule: `DS_ZERO_SUPPLY_MODE=hold` (default `off`)

**Data.**
- The generator's `nodes.csv`: `node_kind` and `member_leaves` per node.
- `member_leaves` is exact per variant. The base-complex `containment.csv`
  would over-include other variants' leaves.
- A bundle without `nodes.csv` reports `hold: no node table` and solves as
  `off`.

**Definitions, for one solve.**
- A leaf L is a *candidate* if it is in `member_leaves` of a node pinned to
  0 in this solve (an observation of 0 that clears the confidence gate).
  Only a pin can remove a supply.
- A node *carries* L if L is in its `member_leaves`, whatever its kind. A
  reaction node carries L if it has an `input` edge from an L-carrier.
- An edge *carries* L if its source and target both carry L and its type
  is `input`, `output`, `assembly`, `dissociation`, `set_member`,
  `variant_pool` or `variant_split`.
  - `catalyst`, `regulator` and `depletion` edges never carry: an enzyme is
    not built into its product.
- An *L-source* is a non-reaction L-carrier with no carrying in-edge.
- *Supplied* means reachable along carrying edges from an L-source that is
  not pinned to 0.
  - A node pinned above 0 is itself a supplied source.
  - A node pinned to 0 neither is supplied nor passes supply on.
- A node *requires* L if it carries L and its kind is `simple_entity`,
  `simple_complex`, `set_variant` or `dissociation_sink`. Pools
  (`set_pool`, `variant_pool`) are OR nodes, and reactions are not
  entities, so neither requires L.

**Hold.**
- Every unpinned node that requires a candidate leaf L and is not supplied
  with L is held at 0, the way cofactors and drugs are held at baseline.
- An observation that pins wins.

**Reporting.**
- The solve reports `zero_supply_rule` (off / hold / hold: no node table),
  `zero_supply_leaves` (candidates with no supply) and
  `zero_supply_held`.
- The bench prints them, and refuses a `hold` run in which nothing was
  held.

**What it does not do.**
- It never raises a value.
- It never touches a supplied node.
- It never acts on a partial knockdown (a pin above 0) or on
  overexpression.
- A protein the network also makes elsewhere (a translation output, any
  source not pinned) stays supplied, so root pinning keeps its meaning.

## Predictions (fixed now)

1. The traced case reads DOWN: gefitinib, EGFR Y1092, Signaling by EGFR.
2. Signaling by EGFR, held-out decryptM, kinobeads map: 18 of 29 scored
   cases are predicted UP. 13 are measured DOWN, and they move off UP. The
   5 measured unchanged are not predicted.
3. An overexpression or partial-knockdown solve is bit-identical to `off`.
   This is asserted on the catalog: 0 changed predictions among
   overexpression cases.
4. Production code defaults are unchanged until adoption.

## Arm and gate

- Arm `zs`: `--server DS_ZERO_SUPPLY_MODE=hold` on the canonical build
  `20261012-1401_01206ca`, by `scripts/run_arm.sh`. The control is the
  canonical scoring `1f58a1b`, re-run at the arm's commit if the code
  differs in anything but this rule.
- Reported:
  - curator tuning (FIT) and held-out (SELECT), and experimental (FIT);
  - net, McNemar p and macro-F1;
  - pathways moved and distinct readouts in the discordant cases;
  - the knockout→UP subgroup.
- **By-construction cases are reported apart.** A case whose readout nodes
  all require the knocked-out leaf reads DOWN because of the rule, not the
  propagation. Gains are given with and without those cases.
- decryptM 1 (FIT) is re-scored and reported. It motivated the rule, so it
  is not evidence for it.
- **Adopt** if:
  - curator held-out net ≥ 0 and experimental net ≥ 0;
  - neither axis is significantly negative (McNemar p < 0.05 against);
  - the gain does not come only from by-construction cases.
- Otherwise it stays off. Every arm counts in the specs/051 SELECT tally.
- **Traced either way:**
  - the largest-gaining pathway;
  - any pathway that loses more than 5 cases.

## Non-goals

- No generator change.
- No change to pinning protocol (root_cycle stays).
- No partial-supply bound (a cap at external supply). That would be a
  different rule.
