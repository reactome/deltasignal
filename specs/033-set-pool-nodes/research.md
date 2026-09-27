# specs/033 — One node per set-valued catalyst or regulator

**Status:** pre-registered 2026-09-26, before any code or arm.

## The defect

A Reactome catalyst or regulator can be an **EntitySet**, meaning any one of
these members plays the role. An example is "RAF activating kinases"
(R-HSA-5672734), the catalyst of "Phosphorylation of RAF": JAK2, SRC, BRAF,
RAF1, ARAF, CaMKII and three phosphorylated MAP2K dimers.

**What the generator does now** (`append_regulators`, subunit/member
decomposition, a design chosen 2026-07-15 because bundling lost single-gene
propagation):
- It attaches **every member as its own edge on every reaction copy**, marked
  AND for catalysts and positive regulators and OR for negative regulators.
- Nothing is copied per member.

**So the reaction node computes**
`inputs × JAK2 × SRC × … × MAP2K dimer`:
- two members up multiply;
- one member at 0 kills the reaction;
- a member made downstream of the reaction makes the reaction require its own
  product.

**Traced in RAF/MAP kinase (specs/032):** the three MAP2K dimers are made only
by "Dissociation of RAS:RAF complex", downstream of RAF phosphorylation. The
component self-requires, does not converge, and ends at 0 under every RAS/RAF
perturbation (experimental 28.6%, against MP-BioPath's 93.9%).

**Inventory (build `20260926-1221_590301c`, groups = reaction copy × set):**

| role | groups | current member flag |
|---|---|---|
| catalyst | 1,881 | AND |
| positive regulator | 56 | AND |
| negative regulator | 146 | OR |

- 63 pathways have at least one catalyst group.
- These edges close only 201 of 8,064 cyclic nodes (RAF 84).

**What the ground truth says** when the perturbed gene is itself a member of a
catalyst set:

| axis | cases | accuracy | expected change |
|---|---|---|---|
| curator | 1,030 | 70.7% (KO) / 63.9% (UP) | 68% |
| experimental | 53 | 57.7% (KO) / 29.6% (UP) | 96–100% |

For comparison, other genes score 85–87% on curator, with 38% expected change.
**So a single member is expected to matter.** A plain mean over members (which
makes one member of nine nearly invisible) contradicts the ground truth. The
failures are mostly missed changes (curator 125 DOWN read as NORMAL, 133 UP
read as NORMAL), 208 of them in Mitotic G1.

## The design

**Generator (`LNG_SET_POOL=1`):** for a bare EntitySet catalyst or regulator
(not a modifier-isoform set):
- emit **one pool node per set per pathway**, with the set's stable id, which
  is the entity the diagram draws;
- add a `set_member` edge (pos, OR) from each member to the pool, emitted once;
- add **one** edge from the pool to each reaction copy, in the original role
  and flag.

Every member reaches every reaction the set serves, exactly as before, so
reachability and cyclic structure are unchanged. Complexes are untouched:
under complex-as-node they are already one node.

**Solver (`DS_SET_POOL_MODE`):** a node whose incoming edges are all
`set_member` is a pool, and its fold is combined from its members' folds.

| mode | combination | one member KO | one member up 80x | several members up |
|---|---|---|---|---|
| `product` | ∏ folds, capped as hill_sat is | 0 | 80x | multiply |
| `extreme` | the member with the largest \|log fold\| (a tie goes to the lower fold) | 0 | 80x | the largest |
| `geomean` | (∏ folds)^(1/n) | 0 | 80^(1/n) | damped |
| `mean` | arithmetic mean (the existing OR rule) | (n−1)/n | ≈80/n | diluted |

- `product` on pooled networks is close to today's arithmetic for catalysts
  and positive regulators.
- For negative regulators the pool is a single inhibitor term, where today
  each member divides separately.
- Default: `product` until adopted.

## Pre-registration

**Arms,** all through `scripts/run_arm.sh`, on one pooled build plus the
canonical build:
- `base`: the canonical build `20260926-1221_590301c`, code defaults;
- `pool_product`, `pool_extreme`, `pool_geomean`, `pool_mean`: the pooled
  build with each mode.

**Choosing the mode (no held-out data used):**
- the winner is the mode with the best **tuning-split** curator macro-F1 plus
  experimental macro-F1 (each on its own case set; summed);
- ties go to the simpler mode, in the order product < extreme < geomean < mean.

**Adopt** the pooled generator plus the chosen mode against `base` only if all
of these hold:
- curator held-out net ≥ −15 (the regeneration noise floor) and not
  significantly negative (p ≥ 0.05, or net ≥ 0);
- experimental net ≥ 0;
- tuning curator net ≥ −15.

It is adopted on faithfulness (one node per curated participant) as long as
it is not harmful. That is the same basis as specs/020's variant sharing.

**Predictions:**
- RAF experimental rises under `extreme` and `geomean` (the loop stops
  requiring its own product);
- `pool_product` moves few predictions against `base`, apart from the TP53
  uuid draw;
- `pool_mean` loses set-member cases on both axes.

**Also reported:**
- per-pathway net on both axes;
- concentration;
- convergence counts;
- the set-member-case accuracy table above for each arm;
- cyclic node count (it must be unchanged);
- node and edge counts.
