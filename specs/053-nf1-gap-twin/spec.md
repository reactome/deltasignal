# specs/053: a closed recycling copy is joined to its fed identity twin

**Status:** registered 2026-10-13, before any generator code was written or
any arm run. Trace: research.md T1–T5. Choice of generator over solver:
Adam, 2026-10-13 (research.md T5, option 2).

## Why

An NF1 knockout cannot reach RAS (research.md T1–T3). The RAS GAP is a copy
of SPRED dimer:NF1 whose only producer is its own release step
(R-HSA-5658231), while the copy that R-HSA-5658438 makes from NF1 is a
different uuid with the identical `nodes.csv` identity. Reactome holds one
SPRED dimer:NF1. Nothing curated links "SPRED dimer binds NF1" to "RAS GAPs
bind RAS:GTP" (R-HSA-5658435 has no precedingEvent), so the generator's
position-aware uuids keep them apart.

## The rule (generator, `LNG_CLOSED_TWIN_JOIN`, default off)

After the logic network is built and before any table is written:

1. **Closed copy.** An entity node X (not a reaction node) is *closed* when
   it has at least one in-edge, every in-edge is an `output` edge, and
   every producing reaction has, among its `input` edges, a node whose
   `member_leaves` strictly contain X's (X is released from a complex it
   was bound into). Leaves are `member_leaves`, or the node's own entity
   when it has none.
2. **Fed twin.** A node W of the same pathway with a `nodes.csv` identity
   equal to X's in every column but the uuid (`node_kind`,
   `diagram_entity_id`, `compartment`, `member_leaves`, `source_sets`,
   `chosen_members`), that is not itself closed and has an in-edge.
3. **Join.** X is merged into W: every edge to or from X is re-pointed to
   W, exact duplicate edges are dropped, and X disappears from every
   output table. With several fed twins, X joins the one with the smallest
   uuid after sorting by entity key, so the build does not depend on
   iteration order; the others are left as they are.
4. All side tables (`nodes.csv`, `pools.csv`, `pool_carriers.csv`,
   `carrier_loops.csv`, `containment.csv`, mappings) are computed on the
   joined network. A joined copy has an outside producer, so it is no
   longer closed for specs/049, and it is a carrier with a non-release
   producer for specs/039.
5. `BUILD.json` / the generator log records the number of joins per
   pathway.

Not done: no edge is added between copies that are both fed (that is the
UUID-silo bridge, closed at 1:3.9 against), and a closed copy with no fed
twin is left alone (specs/049 `closed_twin` handles root copies).

**Measured scope on the canonical build** (research.md T4, DS-side census
on `20261012-1401_01206ca`, so an expectation, not a check of the code):
1,752 closed nodes, 690 with a fed twin, in 10+ pathways (largest: Mitotic
G1 / S-phase family R-HSA-73894 164, R-HSA-5693532 116, R-HSA-983169 95).
Cases whose pins reach a fed twin but not its closed copy: 556 curator (472
correct today), 51 experimental (32 correct). RAF: 5 joins, the 17 NF1
cases.

## Arms

Generator build `C` = canonical generator commit + this rule, `LNG_CLOSED_TWIN_JOIN=1`,
same Reactome release. Solver: code defaults at HEAD of registration.

| arm | build | solver | compared against |
|---|---|---|---|
| `cj0` | canonical `20261012-1401_01206ca` | defaults | — (control; must reproduce `1f58a1b`) |
| `cj` | `C` | defaults | `cj0` |

Both by `scripts/run_arm.sh`, curator and experimental axes.

## Predictions

- **J-P1 (structural).** The build's join count is within 10% of 690, and
  no pathway outside the census list changes edges.
- **J-P2.** NF1 knockout: at least 9 of the 11 cases (6 curator, 5
  experimental) read UP, through the GAP: one case traced to show the GAP
  carrier reading NF1's value.
- **J-P3.** NF1 overexpression: at least 4 of the 6 curator cases read
  DOWN.
- **J-P4.** No prediction changes in a pathway with no join.
- **J-P5.** No direction is predicted for the other ~540 touched cases.
  Joins add producers to recycling copies inside cyclic components, the
  family that hurt in specs/018 (welds, held-out −185); TP53 and DSB repair
  are named in advance as where a loss would come from.

## Decision rule

The faithfulness gate of specs/049 amendment 5, against `cj0`:
- curator held-out net > 0 at McNemar p < 0.01;
- experimental net loss not significant (p ≥ 0.05);
- no pathway losing more than 10 cases on either axis;
- J-P1 and J-P4 holding;
- every experimental loss traced to a named mechanism.

Held-out is SELECT (specs/051), so passing chooses the rule; it is not a
test result. If the gate passes, `LNG_CLOSED_TWIN_JOIN=1` is recommended as
the generator default, and the change is Adam's call. If it fails, nothing
is adopted and the NF1 cases stay wrong, recorded as such.
