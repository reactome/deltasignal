# specs/049: pools that share intermediates — probe record

Continues specs/048 lever B (Amendment 1). Probes only: no arm has run, and
nothing here is pre-registered. Canonical build `20261010-1415_c0e1b2a`,
solver main `f38d4ba`.

## P1. Which pools collide (RAF, R-HSA-5673001)

Generator branch `049-shared-intermediate-pools` = main plus the specs/048 B1
and B2 flags (`LNG_POOL_JOIN_MEMBERS=1`, `LNG_POOL_REGEN=protein`), ported
unchanged from `wip/048-mapk-pools`.
Script: `~/deltasignal-catalogs/analysis/049/shared_pools_probe.py`.

- `find_pools` finds 6 candidate pools in the first pass and keeps 4. It
  drops 246 shared nodes, which reproduces specs/048.
- The two dropped pools are MAPK3 (states R-HSA-109842 / R-HSA-59284) and
  MAPK1 (R-HSA-109853 / R-HSA-59282). Each has 2 states, 3 intermediates,
  2 paths and 254 nodes.
- **They share 246 of their 254 nodes: 243 step reaction copies and the 3
  intermediates. They share no state.** Both cross the same capped
  RAF:MAP2K:MAPKs complexes. The bundle node stands for the curated set
  MAPKs (R-HSA-169291), not for one protein.
- The four RAS GDP/GTP pools share nothing.

## P2. Keep shared intermediates in every pool that claims them

- **Generator:** `LNG_POOL_SHARED=1`. A node claimed by two pools is removed
  only if it is a state of either; intermediates and step copies stay in
  both.
- **Solver:** `update_pools!` (src/core/cycles.jl). Pools of one component
  update in order, as before. A node that two of them write is written once,
  after all of them. Its value is the baseline-flow-weighted mean over every
  path of every pool through it, which is the single-pool rule applied
  across pools.
  - Each pool's flows sum to 1, so the pools weigh equally. That matches the
    `mean` a set pool uses.
  - With no shared node the code path is the old one.
  - The solve reports `cycle_shared_nodes`.
- **RAF build:** `~/deltasignal-catalogs/analysis/049/raf_b12s`. It has 6
  pools (12 states, 54 intermediates), and 0 shared nodes are dropped.

**RAF experimental, the probe against canonical:** 28/49 against 28/49, net 0.

| class | n | canonical | probe |
|---|---|---|---|
| ARAF/BRAF/RAF1 knockout → p-MAPK readouts | 9 | 45–70× UP (wrong) | 0.007–0.014× DOWN (**fixed**) |
| ARAF/BRAF/RAF1 overexpression → p-MAPK readouts | 9 | 31–45× UP (right) | 0.66× DOWN (**broken**) |
| NRAS knockout → p-MAPK | 3 | 36–54× | 2.2–2.4× (still wrong) |
| HRAS/KRAS/NRAS overexpression | 15 | right | right (3–13× in place of 13–100×) |

Converged: 8 of 14 RAF solves (the canonical build converges in 215 of 244
experimental solves).

## P3. Why RAF overexpression reads DOWN: the 14-3-3 cycle has no steady state

Traced RAF1 overexpression (80×) to readout 1268261.
Scripts: `trace_dir.jl`, `node_inedges.jl`, `leak_scan.jl`; a temporary
`DS_POOL_DEBUG` print, not committed.

1. **The MAPK pool's phosphorylation path has drive 0.44.** Its steps are
   [0.044, 10, 1, 1]. The first step is R-HSA-5672972 (MAPKs bind the
   activated RAF complex).
2. **That step's partner is low.** R-HSA-5672718::pool (p21 RAS:GTP:activated
   RAF dimers) reads 0.0275×.
3. **The partner is fed by a cycle that has collapsed.** The cycle is
   R-HSA-2028645 (YWHAB dimer) → 5672954 → 5672696 → 5672957 (PP2A
   dephosphorylates KSR1) → 5672692 → 5672958 → YWHAB dimer, together with
   5672960 (YWHAB dissociates from dephosphorylated RAF). The 14-3-3 dimer is
   a carrier: nothing outside the cycle produces it, and every step releases
   it unchanged.
4. **The collapse is internal.** 1,018 nodes sit at 0.03–0.06×. No negative
   edge from a source above 1.2× enters that set, and only one positive input
   below 0.8× does: the MAPK pool's own step.
5. **The value is set by the iteration cap, not by a fixed point.** On the
   canonical network, YWHAB dimer under RAF1 overexpression reads:

   | `DS_MAX_ITERS` | YWHAB dimer | readout 5674385 |
   |---|---|---|
   | 100 | 0.55× | 63.1× |
   | 500 (default) | 0.047× | 63.1× |
   | 2000 | 0.0077× | 64.4× |

   The canonical prediction does not move only because the MAPK loop rails,
   whatever its input. It is right for the wrong reason.

**Conclusion.**
- Shared-intermediate pools do what specs/048 asked of them: knockout of the
  activating kinase now lowers p-MAPK.
- They also remove the rail that hid an undecided upstream cycle. The
  overexpression cases then read whatever 500 sweeps leave of the 14-3-3
  carrier.
- **Lever B cannot be scored until that cycle has a defined steady state.**
- The 14-3-3 cycle is a binding carrier, not a modification. `find_pools`
  rejects it by design: its core-only forms carry one signature, a
  "carrier loop" (specs/039 amendment 3).

## P4. How many canonical predictions the iteration cap decides: none

Arm `iters2000` (`--server DS_MAX_ITERS=2000`) on the canonical build,
solver `a2413d2`. That commit is this branch, which is byte-identical to main
without shared pool nodes. Results are in
`builds/20261010-1415_c0e1b2a/results/a2413d2/iters2000`.

- **No prediction changes on either axis.**
  - Curator: 0 of 23,511 valid cases; macro-F1 moves in the fourth decimal.
  - Experimental: 0 of 845.
- More solves converge: curator 1,658 of 1,725 and experimental 219 of 244,
  against 215 of 244 at the default.
- **Reading:** cycles like 14-3-3 are undecided in value but not in any
  current call. The rails they feed saturate first.
- So the iteration cap is not a lever on canonical. It becomes one only when
  something removes a rail, as the MAPK pools do.

## P5. Carrier-loop census, and which loops the iteration cap decides

Canonical build, generator `049-shared-intermediate-pools` with every flag at
its default. Scripts are in `~/deltasignal-catalogs/analysis/049/`:
- `carrier_census.py` captures each set `find_pools` rejects as a carrier loop
  on its first pass.
- `drift.jl` solves every valid benchmark perturbation at `max_iters` 500 and
  2000. It pins the benchmark's own `gene_uuids` at 0 or 0.8, and lists nodes
  that differ by more than 2× between the two caps.

**Census.**
- The 17,515 `carrier_loops` in the stats count one per protein R. They are
  **885 distinct loops in 54 pathways**.
- Most are small: the median loop has 3 forms and the largest has 81.
- **279 are closed**: no form of the loop has an in-edge from outside the
  loop's own step reactions. The other 606 are fed.
- Per experimental pathway (loops / closed):

  | pathway | loops | closed |
  |---|---|---|
  | RAF | 17 | 17 |
  | ERBB2 | 53 | 2 |
  | PIP3 | 4 | 3 |
  | HDR | 30 | 14 |
  | Cell Cycle Checkpoints | 190 | 28 |
  | Mitotic G1, S Phase, Prophase | 6, 4, 1 | 1, 1, 1 |
  | WNT, TP53 | 0 | 0 |

**Drift.**
- **Experimental:** 5 of 193 solves have any node that moves between the caps.
  - The three RAF overexpressions (1,634 nodes each). Of these, 9 are the forms
    of the closed YWHAB loop (R 48889) and the rest are downstream of it. Every
    node falls (median 0.049× → 0.008×), and the solves converge only at 2000.
  - PIP3, TP53 knockout and overexpression: the same 32 nodes in both, so the
    drift does not depend on the perturbation. 3 of them are forms of a fed
    loop (R 61460).
- **Curator:** 18 of 1,699 solves.
  - The 5 above.
  - Fanconi: FANCA and FANCC overexpression, 13 nodes, a closed 2-form loop
    (R 144344).
  - Ten more drifts involve no closed carrier loop:
    - DSB Repair and DSB Response: MDC1, RNF8 and RNF168 overexpression, 281
      to 4,159 nodes. Fed loops at most.
    - EGFR: PLCG1, SH3GL1, ARHGEF7 and EPS15, 21 to 1,738 nodes. Fed loops at
      most.
- None of these changes a prediction (P4).

**Reading.**
- Closed carrier loops are common (279), but almost all are decided: at
  steady state their gain stays at 1 under every benchmarked perturbation.
- Only two are undecided under any benchmarked perturbation:
  - 14-3-3 in RAF, which blocks lever B;
  - one 2-form loop in Fanconi.
- The large undecided sets in DSB and EGFR are other cycles, not carrier loops.
- **So a conservation rule for closed carrier loops would cover 279 loops by
  construction, but it is measurable on one experimental class (RAF
  overexpression under shared pools).** A pre-registration must therefore
  state its effect on canonical predictions in advance. Expected: none.

## Next

1. Design for conserved carrier cycles (14-3-3, scaffolds), which
   specs/035 refuted as a blanket "hold at baseline". The question is
   whether a carrier cycle can be solved by conservation, as a modification
   pool is, without holding signal carriers. P5 scopes it to closed loops
   (279), of which only 14-3-3 is known to matter.
