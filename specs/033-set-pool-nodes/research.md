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

## Pre-registration amendment (after review, before any arm)

An adversarial review of both code changes found three things the
pre-registration did not state. They are fixed or specified here, **before
any arm has run**:

1. **Readouts and pins are unchanged.** A pool node carries the set's stId and
   its members, and the benchmark indexed every node under both.
   - A set-valued readout that is also a catalyst set (p-T,Y MAPK dimers,
     p-AKT, GSK3, p-MAPK8/9/10) would have resolved to the pool instead of
     the members' mean.
   - Every member readout would have gained the pool in its mean.
   - Both benchmark loaders now skip `set_pool` rows (test in
     `bench/test_analysis_numbers.py`), so the arms measure the model, not a
     change of readout.
2. **Depletion.**
   - A pooled phosphatase or ubiquitin-ligase set now depletes its substrate
     through its pool: one edge instead of one per member. Under `product` that
     is the same fold; under `extreme` or `geomean` the depletion strength
     follows the mode. This is part of the design and is measured with it.
   - The "a catalyst does not deplete itself" guards compared stIds, and a
     pool's stId is the set's, so a pool would have depleted its own member.
     This happens in 4 E3-ligase sets and would have closed new cycles.
   - Fixed in LNG 81499ef: a pool never depletes one of its members.
3. **The cycle gate** is restated as **cyclic nodes excluding pool nodes must be
   unchanged**. A pool on an existing loop (RAF's is the motivating one) is
   itself cyclic.
   - On the first pooled build: 8,078 cyclic nodes, of which 14 are pools, so
     8,064 non-pool, equal to the base.
   - It is re-checked on the rebuilt catalog.

**Known approximations of `product` against today:**
- **Caps:** the pool caps the member product at 100x before the reaction
  multiplies. Two members at 20x read 100x, not 400x. The difference appears
  only above the cap.
- **Negative-regulator sets become one divide term.** It is equal to the old
  per-member divides below the cap. Under `DS_INHIBITOR_OR=1` (not the
  default) min over members becomes 1/product.
- **The self-inhibitor rule** (specs/022) can now flag a pooled inhibitor set
  whose containment includes the reaction's input, where before it flagged
  only the containing member.

**Build:** the arms run on a pooled build regenerated at LNG 81499ef, not on
the first pooled build.

## Result — the pooled structure ADOPTABLE with `product` (2026-09-26)

**Setup:**
- Base: canonical `20260926-1221_590301c`. Pooled: `20260926-2229_81499ef_setpool2`
  (299 pools, 1,911 set_member edges, 64 depletion edges from pools).
- Non-pool cyclic nodes: 8,064, equal to the base. **Gate 3 passes.**
- All arms at solver 039ce74, through `run_arm.sh`. Pins are identical across
  arms (1,269 nodes over 864 perturbations).

**Choosing the mode** (tuning curator macro-F1 plus experimental macro-F1, as
pre-registered):

| mode | tuning curator mF1 | experimental mF1 | sum |
|---|---|---|---|
| **product** | 0.7558 | 0.5933 | **1.3491** |
| extreme | 0.7557 | 0.5933 | 1.3490 |
| mean | 0.7401 | 0.5952 | 1.3353 |
| geomean | 0.7341 | 0.5933 | 1.3274 |

**`product` against base:**

| split | net | fixed / broke | p | gate |
|---|---|---|---|---|
| curator held-out | +6 | 10 / 4 | 0.18 | ≥ −15, not significantly negative: **pass** |
| curator tuning | +10 | 10 / 0 | 0.002 | ≥ −15: **pass** |
| experimental | +3 | 3 / 0 | 0.25 | ≥ 0: **pass** |

Held-out for the other modes is reported only; none of it was used to
choose: extreme +9 (p = 0.035); geomean −20 (TP53 −79, IFN-γ −17); mean −6
(TP53 −68).

**Verdict:** the pooled structure plus `product` passes every gate. It is
**adoptable on faithfulness**: one node per curated catalyst or regulator,
the set node the diagram draws, and 10,950 fewer edges, with no measurable
harm.

**Predictions:**
- `pool_product` moves few predictions: **held** (24 of 23,511 curator cases).
- RAF experimental rises under extreme or geomean: **only partly**:

| arm | RAF experimental (of 49) | RAF curator (of 84) | experimental predictions DOWN / NORMAL / UP |
|---|---|---|---|
| base | 14 | 40 | 34 / 12 / 3 |
| product | 17 | 43 | 31 / 12 / 6 |
| extreme | 17 | 43 | 31 / 12 / 6 |
| geomean | 19 | 46 | 29 / 12 / 8 |
| mean | 22 | 48 | 26 / 12 / 11 |

MP-BioPath scores RAF experimental 46 of 49. **The collapse has a second
cause.** "MAP2Ks and MAPKs bind to the activated RAF complex" takes F-actin,
CNKSR2 and Ca2+ as inputs, and in this component they are produced only by
"Dissociation of RAS:RAF complex", downstream. The input is wired to the
recycled, downstream copy. This is the same wrong-copy pattern as specs/030's
joins, and it is traced next.

## Note added at merge (2026-09-27)

Under the hierarchy (specs/030), a root complex with a set component now joins
the set's POOL node (one assembly edge: pool → complex) instead of the member
leaves. So `DS_SET_POOL_MODE` also shapes that assembly. Under `product` the
result is numerically close to the member join. It is included in every pooled
build measured here and in specs/038.
