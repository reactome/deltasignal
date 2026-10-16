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

## Amendment 2 (2026-10-11, before any arm ran): adversarial review of both branches

Two review agents read the diffs. Fixed before any build was kept (the first
K build was stopped and discarded):

**Generator** (`937191d`, `8116f24`):
1. **Under `LNG_POOL_JOIN_MEMBERS`, a loop released through `variant_split`
   read as fed, and its free form was not listed.** The plain node a step
   writes was outside `inside`. This is the RAF shape, and build S needs it.
   The plain node now counts as inside, with tests for both shapes.
2. **The B1 regeneration exemption** (a pool node an input feeds) applied at
   the defaults. It is now gated on `LNG_POOL_JOIN_MEMBERS`, so the defaults
   are main's. Tested.
3. **One row per (free form, source),** closed if any of its loops is. Free
   forms and conflicts are counted once.
4. **Every positive entity in-edge counts as supply** for "closed", not only
   output and the variant/set edges.

**Solver:**
1. **`update_pools!` swallowed a NaN** from a pool. The final residual check
   then missed it. It now returns it.
2. **Order of rules.** The edge drop now runs before the specs/035 conserved
   set and the silo bridges are resolved, so they all see one network.
3. **Only the listed activating edge is dropped.** A depletion or inhibition
   between the same pair stays.
4. **A carrier-table row with a missing value is an error,** and a malformed
   JSON row is an ArgumentError.
5. **A node that is both shared and a pool state** is left to its pool. The
   generator does not emit that case.
6. **Docs reworded.** A closed free form reads baseline unless a silo bridge
   or a boundary edge still feeds it.

Tests: `test_carrier_loops.jl` 113, `test_pools.py` 38.

**Recorded, not changed: what a shared node reads.** It is written as the
mean of each pool's estimate (Σnum/Σden).
- **Fixture:** in the review's `twin` fixture with S's supply at 80×, pool 1's
  estimate of the shared complex is 80× and pool 2's is 8.18×. The complex is
  written as 44×. The binding reaction alone would give about 16×.
- **Why:** each pool re-evaluates its pool-fed inputs with the other
  substrate near baseline (amendment 5).
- **What main does:** on the same fixture, main does not converge, and gives
  8.18× or 80× depending on pool order. The branch converges and does not
  depend on order.
- **Effect:** the modified forms and the readout are the same as main's. Only
  the complex moves.
- **This is a property of the registered design, not a fix.** If a shared
  complex decides an `sctl`/`sclosed` case, the case is traced before the
  result is read.

## Amendment 3 (2026-10-11, registered before these arms ran): carrier loops with drugs held

**Why.** P7 traced ARAF knockout on build S. With 14-3-3 held by `closed`,
the knockout zeroes a drug-bound RAF complex (R-HSA-9653109) that inhibits
R-HSA-5672978, and the step de-represses to 4.6×. `DS_DRUG_MODE=inert`
(specs/032) holds that complex at baseline. specs/032 measured it without
conserved carriers, where the RAF loop collapsed (experimental −5), so this
is a new arm, not a re-run.

**No code change.** Both flags exist. Solver commit: HEAD at registration.

| arm | build | solver | compared against | question |
|---|---|---|---|---|
| `kdrug` | K | `DS_DRUG_MODE=inert` | `kctl` | specs/032 on today's canonical networks (attribution) |
| `kcd` | K | `closed` + `inert` | `kctl` | the two solver rules without shared pools |
| `scd` | S | `closed` + `inert` | `kctl` | the lever: both rules with shared pools |

**Predictions.**
- **D-P1.** `scd`: at least 7 of the 9 ARAF/BRAF/RAF1 knockout cases are
  correct (p-MAPK readouts DOWN), through the BRAF/RAF1 path, not 14-3-3
  drain. One fixed case is traced to confirm that 14-3-3 reads baseline.
- **D-P2.** `scd`: at least 7 of the 9 ARAF/BRAF/RAF1 overexpression cases
  stay correct.
- **D-P3.** `scd`: RAF experimental ≥ 35/49, and experimental net ≥ +5
  (the −1 of `sclosed` plus the knockouts).
- **D-P4.** `kcd`: RAF knockouts stay wrong. On K, P6 traced the ARAF
  knockout readout to the MAPK rail, which shared pools address and the
  drug rule does not. RAF within ±2 of `kctl`.
- **D-P5.** `kdrug`: experimental ≤ 0 and RAF not improved, as in specs/032.
- **D-P6.** Curator held-out: `kcd` and `scd` within 30 of `kclosed` (+81)
  and `sclosed` (+61) respectively. No direction is predicted for `kdrug`.
- **D-P7.** The ERCC1, ERCC4 and CCND1 drain losses of `kclosed` persist in
  `kcd` and `scd`; drugs are not on those paths.

**Decision rule.** The same gates as before, each arm against `kctl`:
experimental net > 0, curator held-out net ≥ −15, no pathway losing more
than 10 cases on either axis, structural checks holding (coverage within 1%,
`Drugs: rule inert` with node-solves held > 0, the carrier line as in P6/P7).
- `scd` passing makes three rules a candidate together: `LNG_POOL_SHARED`
  (with its two companion flags), `closed` and `inert`. Each is recorded as
  adopted only in combination unless it passes alone; none has.
- A canonical change from `scd` needs the generator flags too, which is
  Adam's merge (LNG #111).
- Any loss beyond a floor is traced to a case first.

## Amendment 4 (2026-10-12, registered before this arm ran): a closed carrier is not drained by binding either

**Why.** P10 classified every case `closed` changes against Reactome. Its
four held-out losses that have a causal curated route are one mechanism:
PTPN11 knockout in Signaling by FGFR1–4. Under `closed`, free GRB2
(R-HSA-74686) keeps a depletion in-edge from the CBL feedback complex
R-HSA-5654267 (Activated FGFR:p-FRS2:p-PTPN11:p-CBL:GRB2). The knockout
zeroes that complex, GRB2 de-represses to the 10× cap, and RAS:GTP reads UP.

The `closed` rule rests on bound forms holding a negligible share of a
conserved carrier: free E equals total E, and only the supply sets it. On
the same assumption, binding into a complex cannot lower free E, so a
depletion edge into a closed free form is the same error as the release
edge `closed` already drops. Keeping one and dropping the other was an
oversight in the original rule, not a measured choice.

**The rule.** `DS_CARRIER_LOOP_MODE=closed_supply` does what `closed`
does, and also drops every depletion in-edge (negative, `edge_type`
`depletion`) into a free form of a closed loop, from any source. Inhibitions
are not touched; there are none into closed free forms on build K.
`closed` and `all` are unchanged. The solve reports the dropped depletion
edges in `carrier_loop_edges_dropped`, as it does the release edges.

**Size, counted on build K before the arm:** 204 closed free forms; 1,024
listed release edges; 24 depletion edges into 18 closed free forms in 11
pathways (MHC class I antigen presentation 5, TGF-β receptor 4, FGFR1–4
2 each (GRB2 and SOS1, R-HSA-112199), ERBB2 2, DNA Damage Bypass 2, L1CAM 1,
Mitotic G2-G2/M 1, DSB Response 1).

| arm | build | solver | compared against | question |
|---|---|---|---|---|
| `ksup` | K | `closed_supply` | `kctl` (gate), `kclosed` (attribution) | does the depletion drop remove the PTPN11 losses without new ones? |

**Predictions.**
- **S-P1.** Structural: against `kclosed`, no prediction changes outside
  the 11 pathways listed above. `carrier_loop_edges_dropped` exceeds
  `kclosed`'s by exactly the depletion edges present in each solved network.
- **S-P2.** The four PTPN11 knockout cases in FGFR1–4 that `kclosed`
  broke are correct under `ksup`, and free GRB2 reads 1.0× in the traced
  FGFR1 case.
- **S-P3.** Curator held-out against `kctl` ≥ +85 (`kclosed` +81 plus the
  four). No direction is predicted for MHC class I, TGF-β or ERBB2.
- **S-P4.** Experimental equals `kclosed` (−1). None of the 11 pathways is
  one of the ten experimental pathways except as listed by the bench; any
  experimental change is traced.

**Decision rule.** Two gates, both registered now and both reported.
1. **The registered gate as before**, against `kctl`: experimental net > 0,
   curator held-out net ≥ −15, no pathway losing more than 10 cases on
   either axis, the structural checks above.
2. **A route-qualified gate** (new; P9, P10 and Adam's 2026-10-11 point that
   the ground truth itself must be questioned). Every changed case is
   classified with `route_class.py` against its own pathway and against the
   union of the 92. A case whose expected answer is UP or DOWN and which has
   **no route at all** (class `none`, even leniently) in both its own
   pathway and the union is reported separately and not counted: the
   network cannot be blamed for not moving a readout Reactome does not
   connect. Every other case counts, including `recycled` and `lenient`
   routes. The gate is then the same thresholds on the counted cases.

If gate 1 passes, `closed_supply` is recommended for adoption. If only gate
2 passes, it is recommended with the excluded cases listed by name, and
the decision is Adam's. If neither passes, it is not adopted. Any loss
beyond a floor is traced to a case first.

## Amendment 5 (2026-10-12, registered before this arm ran): a closed carrier reads its own gene

**Why.** P11 traced the AURKA losses, and the ten SMAD7 losses `kclosed`
already had, to one gap. "A closed free form reads its supply" must mean
that when the carrier's own gene is perturbed, the free form follows it. But
the rule holds the free form at baseline unless it is itself pinned, and the
bench pins the gene's root copy, a different node of the same entity. On
build K, 11 of 204 closed free forms have such a root copy (GRB2 in
FGFR1–4, SMAD7 and SMURF2 in TGF-β, AURKA, ubiquitin in ERBB2, ANK1 in L1CAM,
CALR and TAPBP in MHC class I).

**The rule.** `DS_CARRIER_LOOP_MODE=closed_twin` does what `closed_supply`
does, and gives each closed free form an activating OR in-edge from every
root copy of the same entity: a node with the same `reactome_id` and no
in-edge in the network, that is not itself a listed free form. The free
form then reads its root copy (the mean, if several), which is baseline
unless the gene is perturbed. The solve reports the edges added in a new
diagnostic, `carrier_loop_twin_edges`.

| arm | build | solver | compared against | question |
|---|---|---|---|---|
| `ktwin` | K | `closed_twin` | `kctl` (gate), `ksup` (attribution) | does the gene reach its own closed carrier? |

**Predictions.**
- **T-P1.** Structural: against `ksup`, no prediction changes outside the
  seven pathways with a root twin (FGFR1–4, TGF-β receptor, Mitotic
  G2-G2/M, ERBB2, L1CAM, MHC class I). `carrier_loop_twin_edges` > 0 in
  those pathways' solves and 0 elsewhere.
- **T-P2.** The two AURKA tuning cases (p-S95-PHLDA1) are correct again,
  with free AURKA reading the pinned value.
- **T-P3.** At least 6 of the 10 SMAD7 held-out cases that `kclosed`
  broke are correct again.
- **T-P4.** Curator held-out ≥ `ksup` (+124 against `kctl`).
- **T-P5.** Experimental within 1 of `ksup` (−2); none of its changed cases
  has a twin carrier on its path except possibly ubiquitin in ERBB2.

**Decision rule.** Gates 1 and 2 of amendment 4 are computed and reported.
Neither can pass if T-P5 holds, so a third gate is registered, and it is the
one the decision uses. It is the standard the project has adopted under
before: vn7 (specs/046, experimental −10, not significant) and cycle balance
(specs/039 amendment 5, experimental −2, not significant).
3. **Faithfulness gate.**
   - Curator held-out net > 0 at McNemar p < 0.01.
   - Experimental net loss not significant (p ≥ 0.05).
   - No pathway losing more than 10 cases on either axis.
   - Structural checks holding.
   - Every experimental loss traced to a named mechanism: either the case has
     no curated route (class `none` or `lenient` in its pathway and in the union),
     or it is a mechanism already recorded in this spec (RAF: P7/P8).

If gate 3 passes, `closed_twin` is recommended for adoption as the default,
with the losses listed by name; the change of default is Adam's call. If it
fails, nothing is adopted.

## Amendment 6 (2026-10-12): `closed_twin` is the default

Gate 3 passed (research.md P12): curator held-out +122 (p 2e-23),
experimental −2 (p 0.73), no pathway losing more than 10, structural checks
holding, and every experimental loss named. Adam accepted the
recommendation on 2026-10-12, so `DS_CARRIER_LOOP_MODE` now defaults to
`closed_twin`. `off` stays available for the A/B and to reproduce earlier
numbers.

**What the default needs.** `carrier_loops.csv` is written by the generator
from LNG #111 onwards. A bundle without it reports
`closed_twin: no carrier table` and solves exactly as `off` (tested). The
bench refuses it, as it does a build with no `pathogens.csv`. So:
- builds before LNG #111, the canonical `20261010-1415_c0e1b2a` included,
  are re-run with `--server DS_CARRIER_LOOP_MODE=off`;
- the canonical numbers move only when LNG #111 is merged, a catalog is built
  from generator main, and that build is benchmarked under the default. Until
  then, `ktwin` on build K is the measurement of record. K was built at
  8116f24 with the generator defaults, and `kctl` on it reproduced the
  canonical results exactly.

**Known losses carried into the default** (P12): CALR knockout and
overexpression → 983416 (2 tuning cases, a recycled route where the expected
answer is no change); ERCC1/ERCC4 in HDR (the drain recorded in P9); the
SMAD7 and PTPN11 → SPRY2 cases, where an enzyme knockout should raise its
substrate and the network has no edge that says so. That is the next lever,
and it is registered separately before it runs.

**Adversarial review of amendments 4–6 (2026-10-12)**, before merge. Three
confirmed defects, all fixed:
1. Influence scores (API and CLI) were computed on the network as given, so
   a closed free form was credited through the loop and depletion edges the
   solve had dropped, and its root twin got nothing for the edge it solved
   through. Both now score `solved_carrier_network(network)`. On a fixture,
   the twin's influence rises by exactly 1 and the release complex's falls by 1.
2. The bench refused any run where no free form changed, or no twin edge was
   added. Under the default, that refused valid subset runs: 56 of 92 pathways
   on K have no closed loop, and only 9 have a root twin. The bench now
   counts what each network makes possible, by the amendment 5 rule and
   independently of the solver. It refuses only when closed loops are listed
   and nothing changed, or when the twin edges added differ from the twins
   available. Checked on all 92 pathways of K: 204 closed free forms, 16 twins
   by the bench's count and 16 by the solver's, with no pathway differing.
3. `docs/API.md` still gave `off` as the default, and CLAUDE.md omitted
   `carrier_loop_twin_edges`.

One plausible defect, also fixed: an inline network can send
`reactome_id: ""`, which made every id-less root a twin of every id-less
closed free form. An empty id now names no entity. No shipped build has
such a node.

One accepted limit: with several twins, the OR mean sums in uuid-sorted
order. That moves values only at ulp level, inside the specs/013 `pred_ui`
tolerance.

**Canonical rebuild (2026-10-12).** LNG #111 was merged as 01206ca, a tree
identical to 8116f24, which built K. Build `20261012-1401_01206ca`, solver
1f58a1b, default config:
- every curator and experimental prediction equals `ktwin`'s; `pred_ui`
  differs by at most 0.01, which is relabelling;
- curator held-out in release 89.29% / 0.8612;
- experimental 589/849.

docs/RESULTS.md carries the box.

## Amendment 7 (2026-10-15, after the week review, before any code): a carrier shared by several pools excludes every pool's release steps

The week's adversarial review found, and reproduced on the canonical build,
that the specs/039 carrier rule (`DS_CYCLE_CARRIERS=on`, default) excludes
only the CURRENT pool's release steps when it writes a carrier. Each of the
12 RAS GAP free forms in RAF (R-HSA-5673001) is a carrier of all four RAS
isoform pools and has no producer except those pools' four release
reactions. Each pool in turn writes the GAP from the mean of the OTHER three
pools' release flux, which re-closes the gain-1 release loop the rule exists
to cut. The value depends on which pool writes last, the writes disagree
within every sweep, and the RAF solve does not converge (pool1 base state
pinned at 10: GAPs 4.03×, converged false; the documented rule gives 1×).
These are the only multi-pool carriers in the canonical build.

**Rule (in place, under `DS_CYCLE_CARRIERS=on`):** a carrier reads the mean
fold of its producers other than the release steps of EVERY pool that lists
it as a carrier (1 if none), as the rule's own description says. Every pool
then computes the same value, so write order no longer matters.
`cycle_carriers_held` counts distinct carriers (fixed separately, PR #123).

**Arm `carrierfix`:** the code at this amendment (on top of specs/061
amendment 3), against `rtwin58g4` (same code without it), RAF reported
separately since only RAF has multi-pool carriers.

**Prediction:** only RAF predictions move (no other pathway has a
multi-pool carrier); RAF convergence improves; the GAPs read about 1× under
RAS perturbations. Which way RAF's scored cases move I do not know; RAF is
set aside from rule decisions (specs/054), so the decision is on
faithfulness (the documented rule) provided nothing outside RAF moves.

### Amendment 7 result (2026-10-16)

Isolated by an arm with only this change reverted (`9ff8335`, local branch,
`results/9ff8335/iso_nocarrier`) against `76d7306` (`rtwin58g6`):

| | net | fixed / broke | p |
|---|---|---|---|
| curator held-out | 0 | 0 / 0 | |
| curator tuning | +21 | 21 / 0 (all RAF) | < 1e-4 |
| experimental | +14 | 14 / 0 (all RAF) | 1e-4 |

Nothing outside RAF moves, as predicted. 36 curator and 22 experimental RAF
cases go from unconverged to converged. RAS overexpressions (HRAS, NRAS,
KRAS) now read UP and KRAS knockouts DOWN, as the curators and experiments
expect; under the old carrier write the RAS overexpressions read DOWN
(specs/055 noted them as wrong). RAF is set aside from rule decisions
(specs/054), so this is adopted on faithfulness: it implements the rule's
documented behaviour, it converges, and nothing outside RAF changes.
