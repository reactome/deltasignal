# specs/044 research: F7, sets inside root complexes

## Evidence (2026-10-03)

**What it is.** A set component of a root complex was flattened into its
members, and each member was emitted as a required (`and`) assembly input of
the complex.
- Code: `_decompose_hier` → `get_terminal_components`, logic-network-generator
  `src/logic_network_generator.py` around line 1986.
- Curation: "VAV1 Rho/Rac effectors:GDP" (R-HSA-114543) is GDP plus "VAV1
  effectors", a CandidateSet of RHOA, RHOG, RAC1, RAC2 and CDC42.
- Solve, through the API, GPVI: a knockout of any ONE of the five gives the
  complex 0, the same as knocking out all five.

**Scope** (build 20260928-1110_06ccb63, `~/deltasignal-catalogs/analysis/044/f7_scope.py`):
3,828 of 10,165 assembly edges reach their complex only through a set. 557
complex nodes, in 81 pathways, need two or more alternatives at once.

**Benchmark exposure** (`f7_cases.py`, `f7_routes.py`, solver 8d6d27c):

| Axis | Valid cases that perturb a set alternative | Accuracy on them | Accuracy on all cases |
|---|---|---|---|
| Curator | 3,920 (100 genes, 53 pathways) | 88.2% | 85.8% |
| Experimental | 206 (20 genes, 10 pathways) | 78.6% | 70.2% |

Cases that are correct and whose readout is reachable **only** through the F7
edges:

| Axis | Knockout → DOWN | Overexpression → UP |
|---|---|---|
| Experimental | 38 | 46 |
| Curator | 311 | 348 |

Wrong cases that are reachable only through the F7 edges:
- Curator, 74 that a looser rule could fix: knockout predicted DOWN where the
  curators expect NONE (35), and overexpression predicted UP where they expect
  NONE (39). Also 40 knockouts predicted DOWN where the curators expected DOWN
  but the prediction was wrong in another way. 79 in all.
- Experimental: 12 wrong cases in all.

**Interpretation.** The curators and the experiments both behave as if losing
one alternative of these sets lowers the complex. The solver's set-pool rule
agrees: `DS_SET_POOL_MODE=product` (specs/033) multiplies member folds. The
literal "any member suffices" (`max`) was refuted for set members (specs/038:
held-out −120, experimental −57). So F7 is fixed **structurally**:
- the set becomes one node, fed by its members over `set_member` / OR edges,
  and that node is an AND input of the complex;
- the semantics stay with the solver's set-pool rule;
- no "any member suffices" arm is run.

OR between **different reactions** that produce one entity is a different
question: issue #90.

## Change

logic-network-generator `feat/f7-set-nodes` (83a71a2), in a separate worktree
so catalog builds are not disturbed. In `_decompose_hier`, a set component
becomes a set node:
- a member the network already has is joined;
- a member that is a complex is built level by level;
- a member that is a set becomes a nested set node;
- anything else becomes a leaf.

Modifier-isoform sets stay atomic. Everything is in `boundary_edges.csv`.
Verified: GPVI's VAV1 complex now reads GDP + "VAV1 effectors"
(set_member/OR: RHOA, RHOG, RAC1, RAC2, CDC42). Curated edges are unchanged.

## Pre-registration (committed before the arm runs)

- **Arms:**
  - `f7ctrl`: catalog from generator main 1491276, built as variant
    20261003-1058_1491276_f7ctrl.
  - `f7set`: the same generator plus 83a71a2.
  - Both are scored by the same solver commit, with code defaults.
- **Prediction: nearly identical.** Under `product`, a set node's fold is the
  product of its members' folds, and the complex ANDs it with the other
  components. That is the same product as before. Four things can still move a
  prediction:
  1. the 100× cap is applied at the set node as well;
  2. the per-input sensitivity transform is applied at one more hop;
  3. a member that is a complex is now built hierarchically instead of
     flattened (the same product under multiplication);
  4. **a member the pathway already produces is now JOINED to that node, which
     creates a new route**, like specs/030's joins. This is the likeliest source
     of change.
- **Expected:**
  - fewer than 50 changed predictions on the curator axis (of 24,100) and
    fewer than 5 on the experimental axis (of 849);
  - any change traced to one of the four mechanisms above, with the
    discordant cases traced;
  - held-out net change within ±15, the regeneration noise floor.
- **Decision rule:**
  - adopt on faithfulness if held-out and experimental are each within the
    noise floor;
  - a gain needs McNemar p < 0.05 and more than one pathway, or it is not a
    gain;
  - a loss beyond the noise floor on either axis blocks adoption until it is
    traced.
- **Report:**
  - both axes, curator held-out and all;
  - number of predictions changed;
  - how many pathways and distinct readouts the changes span;
  - McNemar p.

## Result (2026-10-03)

Builds:
- `f7ctrl`: 20261003-1058_1491276_f7ctrl, 72,322 nodes / 214,197 edges.
- `f7set`: 20261003-1113_83a71a2_f7set, 73,570 nodes / 215,731 edges.

Solver a952354, code defaults (Gauss-Seidel sweep):

| | Net | Fixed / broken | p | Pathways |
|---|---|---|---|---|
| Curator held-out | +8 | 26 / 18 | 0.29 | 6 |
| Curator tuning | +88 | 105 / 17 | | TP53 +96, Mitotic G1 −8 |
| Experimental | −6 | 4 / 10 | 0.18 | Mitotic G1, TP53 |

**The TP53 swing is node-label noise, not F7.** Two builds with identical
networks (by stable id, curated rows in the same order) but freshly minted
uuids differ by 158 TP53 predictions:
- canonical 20260928 against `f7ctrl`: TP53 −145, held-out +2;
- the cause is that the Gauss-Seidel sweep order inside a cycle follows label
  order (specs/013), and TP53's cycle is the catalog's largest.

So any number for the tuning pathways is a coin flip across rebuilds.
Re-scored with `DS_SCC_SWEEP=jacobi`, which does not depend on labels, the
TP53 change disappears:

| (Jacobi) | Net | Fixed / broken | p | Where |
|---|---|---|---|---|
| Curator held-out | +5 | 12 / 7 | 0.36 | NODAL +8, DDX58 +2 |
| Curator tuning | −9 | 6 / 15 | 0.078 | Mitotic G1 −8, DSB repair −5 |
| Experimental | −5 | 0 / 5 | 0.062 | Mitotic G1 only (MYC, RB1 and RBL1 knockouts) |

**Trace.** Mitotic G1's boundary changes are almost all the same produced
sources rewired through set nodes (cyclin D, CDK4/6, CDKN1A/B/C, cyclin E/A
families), which is identical under `product`. One edge is new (mechanism 4):
Cyclin E:p-T160-CDK2, which the pathway produces, now joins the set
"CCNA:p-T160-CDK2, CCNE:p-T160-CDK2" that it is a member of. That is faithful
to the curation, and it adds a route into the p27↔CDK2 loop (the M1 / M4
motifs of specs/041).

**Against the pre-registration:**
- held-out (+5 / +8) and experimental (−5 / −6) are within the noise floor;
- experimental changes: 5, against "fewer than 5" predicted;
- all changes are traced: the rewiring, plus one new join into a curated loop.

The decision rule permits adoption on faithfulness. The residual loss is in
one tuning pathway and runs through a curated loop.

**Incidental:** `jacobi` versus Gauss-Seidel on `f7ctrl` gives curator all
85.22% vs 85.20%, and experimental 69.47% vs 69.70% (−2 cases). This is the
measurement to start from for making the sweep label-independent by default
(specs/013, code review F14).
