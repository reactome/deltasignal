# specs/049: pools that share intermediates, and carrier loops (pre-registration)

**Goal.** Finish specs/048 lever B: fix the RAF knockout→UP class without
breaking the RAF overexpression cases. **Experimental is the primary axis.**
The probes behind this spec are in `research.md` (P1 to P5). No arm has run.

Control throughout: canonical build `20261010-1415_c0e1b2a`, solver main
(`f38d4ba`) at code defaults, `root_cycle` pins.

## What the probes established

1. **Shared intermediates.** The MAPK1 and MAPK3 pools share 246 of 254
   nodes, the capped RAF:MAP2K:MAPKs complexes, and are both dropped (P1).
   - Keeping shared nodes in both pools (`LNG_POOL_SHARED=1`), with the
     solver writing each shared node once from all its pools
     (`update_pools!`), fixes the 9 ARAF/BRAF/RAF1 knockout cases.
   - It breaks the 9 overexpression cases. RAF net is 0 (P2).
2. **The overexpression break is an undecided carrier loop.** The 14-3-3
   (YWHAB dimer) cycle has no producer outside itself (P3).
   - Its value is set by the iteration cap: 0.55× at 100 sweeps, 0.047× at
     500, 0.0077× at 2,000.
   - On canonical networks the MAPK loop rails and hides it.
   - `find_pools` rejects such sets by design as "carrier loops": one
     core-only signature, so there is nothing to interconvert.
3. **The cap decides no canonical prediction** (P4, arm `iters2000`: 0 of
   23,511 curator, 0 of 845 experimental).
4. **Census** (P5): 885 distinct carrier loops in 54 pathways, 279 of them
   closed. Only two closed loops are undecided under any benchmarked
   perturbation: 14-3-3 in RAF, and one 2-form loop in Fanconi.

## The rule: a carrier loop's free form reads its supply, not its own release

**Derivation.** A carrier loop is one protein R bound and released unchanged:
E → E:X → … → E. Mass is conserved, so free E plus all bound forms equals
the total amount of E.
- **Assumption:** the bound forms hold a negligible share of the total
  (unsaturated). This is the convention specs/039 amendment 2 already uses:
  pool intermediates have no baseline share.
- **Consequence:** free E equals total E, and total E is set by E's supply,
  not by how fast the loop turns.
- **Supply** means E's producers outside the loop: synthesis, or another
  reaction outside the loop. If there are none (a closed loop), the supply is
  constant and free E reads baseline.
- **This is the rule `DS_CYCLE_CARRIERS` already applies** to an enzyme
  inside a modification pool: "the free form reads only its other
  producers" (`src/core/cycles.jl`, `update_pool!`). Here the same rule is
  applied to carrier loops that are not part of a pool.
- **Bound forms are not touched.** They keep reading their reactions, so a
  partner's signal still passes through them (RAF up → RAF:14-3-3 up).

**Not specs/035.**
- specs/035 held every node "released unchanged" at baseline, decided per
  node from the graph. It held PIP3 itself, a second messenger, and lost 122
  experimental cases.
- Here a carrier must be a protein R. Its loop comes from Reactome's own
  reaction participants (`find_pools`'s R-steps), and only its core-only free
  form is affected.
- specs/035 recorded this narrower rule ("require a protein or complex
  entity") as a new pre-registration, not run until now.

## Design

**Generator** (`logic-network-generator`, branch
`049-shared-intermediate-pools`):
- **New table `carrier_loops.csv`**, always written, even when empty. The
  logic network itself is unchanged.
- One row per (loop, free-form node, release reaction node). The columns are
  `loop_id`, `r`, `closed`, `free_uuid`, `free_stid` and `release_uuid`.
- **Loop:** a strongly connected R-form set that `find_pools` rejects as a
  carrier loop, taken on its first pass, de-duplicated across R.
- **Free form:** a loop node whose proteins are exactly the loop's core.
- **Release:** a loop step reaction node with an output edge into the free
  form.
- **Closed:** no form of the loop has an entity-level in-edge from outside
  the loop's step reactions. The in-edges counted are `output`,
  `variant_split`, `variant_pool` and `set_member`.
- **Exclusions,** each counted in the stats:
  - a free form that any pool manages (as a state, intermediate, step node or
    carrier) is not listed (`carrier_loop_conflicts`);
  - a free form claimed by two loops is listed once.

**Solver:** `DS_CARRIER_LOOP_MODE=off` (default) | `closed` | `all`.
- **`closed`:** free forms of closed loops have their release in-edges
  dropped. Those are by definition all of their producers, so each becomes a
  root and reads baseline unless an observation pins it.
- **`all`:** the same for every loop. A fed free form then reads only its
  producers outside the loop.
- **Pins still win.** A gene's root_cycle pin set is computed by the bench
  from the network files, so it is not affected.
- The solve reports `carrier_loop_rule` and `carrier_loop_forms`. A bundle
  without the table reports `<mode>: no carrier table`, and the bench refuses
  that arm.
- `off` must leave the solve unchanged, including under `LNG_POOL_SHARED`
  builds.

**Shared intermediates** (already on the branches, P2): `LNG_POOL_SHARED=1`,
and `update_pools!` in the solver. Builds that use it also need
`LNG_POOL_JOIN_MEMBERS=1`, `LNG_POOL_REGEN=protein` and `LNG_BIND_STOICH=1`
(specs/048).

## Builds and arms

**Builds.** Each is the full 92-pathway catalog from the generator branch,
filed under `~/deltasignal-catalogs/builds/`.
- **K:** every flag at its default. This is the canonical networks plus
  `carrier_loops.csv`.
- **S:** `LNG_POOL_JOIN_MEMBERS=1`, `LNG_POOL_REGEN=protein`,
  `LNG_POOL_SHARED=1`, plus `carrier_loops.csv`.

**Arms.**

| arm | build | solver | compared against |
|---|---|---|---|
| `kctl` | K | `DS_CARRIER_LOOP_MODE=off` | canonical results (identity check) |
| `kclosed` | K | `closed` | `kctl` |
| `kall` | K | `all` | `kctl` |
| `sctl` | S | `off` | `kctl` |
| `sclosed` | S | `closed` | `kctl` (this is lever B) |

## Structural checks (before any score is read)

1. Both builds complete for all 92 pathways.
2. `kctl` reproduces the canonical predictions exactly, on both axes.
3. Census figures for K are recorded before scoring:
   - loops, closed loops, listed free forms and conflicts, per pathway;
   - RAF must list the YWHAB dimer (R-HSA-2028645) as a closed free form.
4. On S, RAF has a MAPK1 pool and a MAPK3 pool, and 0 shared nodes are
   dropped.
5. Coverage (valid cases, `Pinned:`) is within 1% of `kctl` on both axes.

## Predictions

**`kclosed`, against `kctl`:**
- **K-P1.** RAF experimental is unchanged. The MAPK loop rails over 14-3-3
  whatever 14-3-3 reads.
- **K-P2.** Changes are confined to pathways with a listed closed free form.
  Experimental net is within ±5.
- **K-P3.** Curator held-out net is within ±20. No direction is predicted.
  specs/035's +42 held-out came from a broader set in HDR and ERBB2 and is
  not evidence for this one.

**`kall`, against `kctl`:**
- **K-P4.** More pathways move than under `kclosed`. No direction is
  predicted on either axis. This arm is the record of whether "closed" is
  the right boundary.

**`sctl`, against `kctl`:** this reproduces P2 at catalog scale.
- **S-P1.** The 9 ARAF/BRAF/RAF1 knockout→UP cases become correct.
- **S-P2.** The 9 ARAF/BRAF/RAF1 overexpression cases become wrong. RAF net
  is within ±2.

**`sclosed`, against `kctl` (lever B):**
- **S-P3.** The 9 knockout cases stay fixed.
- **S-P4.** At least 7 of the 9 overexpression cases are correct (p-MAPK
  readouts UP).
- **S-P5.** RAF experimental is ≥ 35/49 (canonical 28/49). Experimental net
  is ≥ +7.
- **S-P6.** The NRAS knockout cases stay wrong (P2: 2.2–2.4×). That is
  another mechanism.
- **S-P7.** Curator held-out: no direction is predicted. MAPK cycles exist in
  many signalling pathways.

## Decision rule

Each candidate (`kclosed`, `kall`, `sclosed`) is decided against `kctl` on
its own.
- **Adopt** if:
  - experimental net > 0 (for `kclosed` and `kall`, ≥ 0 together with a
    held-out net > 0);
  - curator held-out net ≥ −15;
  - no single pathway loses more than 10 cases on either axis;
  - the structural checks hold.
- **If more than one is adopted,** only `sclosed` is a candidate for
  canonical. Its `closed` rule must also pass alone (`kclosed`). A rule
  adopted only as part of a combination is recorded as such.
- **Any loss beyond a floor is traced to a case first.**
- **Report:** both axes, held-out and tuning, McNemar p, the pathways and
  distinct readouts moved, convergence (and `DS_MAX_ITERS` drift for RAF,
  re-run with `drift.jl`), and each prediction marked met or missed.

## Not in scope

- Saturated carriers, where bound forms hold a real share. That would need a
  baseline bound fraction, a new parameter.
- The undecided non-carrier cycles of P5: DSB repair and response (MDC1,
  RNF8, RNF168) and EGFR (PLCG1, SH3GL1, ARHGEF7, EPS15).
- The leftover specs/047 classes (PTEN mRNA, the Mitotic G1 E2F silo, TP53
  DAXX, PIP3 THEM4) and NRAS knockout.

## Amendment 1 (2026-10-11, before any arm ran): what was built

Two details of the design changed in implementation. Neither touches a
prediction or the decision rule.

- **The table lists in-edges, not release reactions.**
  - Column `loop_source_uuid` replaces `release_uuid`. It holds every node
    inside the loop with an entity-level in-edge into the free form.
  - Usually that is a loop step reaction node. Under a variant cap it can be
    a plain node that hands the free form on by `variant_split`, which is
    what the solver must drop.
  - In a closed loop that is every in-edge of the free form.
  - Generator `465de02`: `find_pools(carrier_loops=…)`.
- **Free forms that pool carriers also hold are excluded.** As in the design,
  a free form that a pool manages is excluded. Pool carriers count as managed
  too, so the RAS GAPs in RAF stay with the RAS pools.
  - On the canonical RAF networks this excludes 12 of the 17 loops.
  - 5 free forms are listed, the YWHAB dimer (R-HSA-2028645) among them.
- **The solver drops the listed edges** before the solve
  (`src/core/carrier_loops.jl`, `apply_carrier_loops`).
  - Tests: `test/test_carrier_loops.jl` (106 assertions), and the
    shared-pool fixture in `test/test_cycles.jl` (24).
  - The bench sends the table with a POSTed network. It refuses a
    `closed`/`all` arm whose build has no table, and an arm in which no solve
    changed a free form.
