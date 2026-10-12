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

## P6. Arms `kctl` and `kclosed` (build K, solver `3c64f21`)

Build K is `20261011-1609_8116f24_carrierK` (generator `8116f24`, flags at
their defaults). Results are in `builds/<K>/results/3c64f21/{kctl,kclosed}`.

**Structural checks:**
1. Both builds are complete, 92 of 92 pathways each.
2. **`kctl` reproduces the canonical predictions exactly:** 0 of 24,101
   curator rows and 0 of 850 experimental rows differ.
3. **K census:**
   - 2,083 rows and 1,087 free forms. 204 free forms are closed, in 36
     pathways.
   - K's pool tables match canonical, row for row in every pathway, 122 pools
     in all.
   - RAF lists 5 closed free forms, the YWHAB dimer R-HSA-2028645 among them.
4. **S:** RAF has 6 pools, MAPK3 (pool5: R-HSA-59284 / 109842) and MAPK1
   (pool6: R-HSA-59282 / 109853) among them, and 0 shared nodes are dropped.
5. Coverage is identical on both axes.

**`kclosed` against `kctl`:**

| split | net | fixed / broke | p | pathways |
|---|---|---|---|---|
| curator held-out | **+81** | 102 / 21 | 5.8e-14 | 10 (+9 / −1) |
| curator tuning | +88 | 98 / 10 | 2.7e-19 | 4 (+3 / −1) |
| **experimental** | **−1** | 3 / 4 | 1 | 4 (+1 / −3) |

- Curator, best pathways: HDR +64, ERBB2 +25, NER +18, PDGF +18. Worst:
  TGF-β −9, RAF −3.
- Converged: curator 1,646 of 1,725, experimental 214 of 244 (`kctl`: 215).
- Changed: 3,902 curator and 1,194 experimental free-form solves.

**Predictions:**
- **K-P1 missed:** RAF experimental lost 1 case (ARAF KO).
- **K-P2 met:** experimental within ±5.
- **K-P3 missed:** held-out +81, outside ±20.

**Decision: not adopted on its own.** The rule requires experimental ≥ 0
together with a positive held-out; experimental is −1.

**Every experimental change traced** (`diff_walk.jl`: the solve under `off`
and under `closed`, walking back from the readout).

In each case, under `off` a closed carrier drained to about zero when the
perturbation cut the loop's flow. That drain was the only route from the gene
to the readout.

- **ERCC1 KO and ERCC4 KO → R-HSA-5686483 (HDR), broke.**
  - The pin set holds ERCC1:XPF (R-HSA-109943), a closed free form. The pin
    wins, as designed.
  - Under `off`, RPA (R-HSA-68462 and phospho-RPA R-HSA-5685169, both closed
    free forms) drained to 2e-6×, and the readout followed.
  - Under `closed`, RPA holds at baseline, and nothing else connects ERCC1
    to the readout, so it reads 1.0.
- **BRCA2 KO (HDR), fixed** by the same mechanism in reverse: RPA no longer
  drains.
- **ARAF KO → R-HSA-5672718 (RAF), broke.**
  - Under `off`, the 14-3-3 dimer drained to 0.001×, and the readout followed
    (0.02×, correct).
  - Under `closed`, 14-3-3 holds at baseline. The readout then reads
    R-HSA-5672734::pool at 12×, which is the railed MAPK loop, the specs/048
    knockout→UP class.
  - This is the class `sctl`/`sclosed` target.
- **CCND1 overexpression → R-HSA-157443 (Mitotic G1), broke.**
  - Under `off`, the closed free form R-HSA-8848007 drained to 2e-6×.
  - Under `closed`, it holds at baseline. The readout nodes it fed move from
    about 0 to 1.0, and the case aggregate falls from 43× to 1.0.
- **EGFR KO and PTK6 KO (ERBB2), fixed:** ratios 0.87 → 0.69 and 0.88 → 0.62.

**Reading.**
- On canonical networks, carrier drain is a route the propagator relies on.
  A partner's knockout empties the carrier, and every complex downstream of
  the carrier follows.
- That is right for the wrong reason where the true route runs through the
  partner (ERCC1:XPF), and simply wrong where it does not (BRCA2).
- Held-out says the conserved reading is better overall (+81, 102 fixed
  against 21 broken). Experimental is balanced (3 fixed against 4 broken).
- `kall` and the S arms are still running.

## P7. Arms `kall`, `sctl` and `sclosed`; decision

Every arm is scored against `kctl`.
- `kall` used solver `3c64f21`.
- `sctl` and `sclosed` used `4b1ecf7`. That commit adds only research.md to
  `3c64f21`, so the solver code is identical.
- Build S is `20261011-1639_8116f24_sharedS`.

| arm | curator held-out | curator tuning | experimental | worst pathway |
|---|---|---|---|---|
| `kclosed` | +81 (102/21), p 6e-14 | +88 | −1 (3/4) | TGF-β −9 |
| `kall` | **−55** (180/235), p 0.008 | +101 | −7 (8/15) | **IFN α/β −168** |
| `sctl` | −15 (13/28), p 0.03 | −49 | −6 (9/15) | **HDR −49** |
| `sclosed` | +61 (110/49), p 1.5e-6 | +44 | −1 (2/3) | TGF-β −9, Class I MHC −9 |

RAF experimental reads 28/49 under `kctl`, `sctl` and `sclosed` alike.

**Predictions:**
- **K-P4 met:** more pathways move under `kall`: 28, against 14 under
  `kclosed`.
- **S-P1 met:** under `sctl` the 9 ARAF/BRAF/RAF1 knockout cases are fixed.
- **S-P2 met:** under `sctl` the 9 overexpression cases break. RAF net is 0.
- **S-P3 missed:** under `sclosed` all 9 knockout cases are wrong again.
- **S-P4 met only trivially:** the 9 overexpression cases are right because
  nothing in RAF changed against `kctl`.
- **S-P5 missed:** RAF stays at 28/49, and experimental net is −1.
- **S-P6 met:** the NRAS knockout cases stay wrong.

**Decision: no arm is adopted.**
- `kclosed` and `sclosed` fail on experimental, −1 each.
- `kall` and `sctl` fail on every gate.

### Why lever B failed: the knockout fix was also a carrier drain

Traced ARAF knockout to readout 1268261 on build S (`diff_walk.jl`,
`node_in2.jl`).

**Under `sctl`, the readout reads 0.009×, correct, but for the wrong reason.**
- ARAF knockout drains the 14-3-3 carrier loop (P3).
- With it go every RAS:RAF complex and every step copy of the MAPK pool's
  phosphorylation path: R-HSA-5672978 copies read 0.009×.
- So the specs/049 P2 "knockouts fixed" was the same drain, not the removal
  of ARAF from the path.

**Under `sclosed`, 14-3-3 holds at baseline, and the readout reads 2.1×,
wrong.**
- The MAPK pool's path now runs through the BRAF/RAF1 complexes, which is
  biologically fair: they are redundant.
- But 20 of the 30 step copies of R-HSA-5672978 read **4.6×**.
- Each copy has a negative regulator, R-HSA-9653109: a drug (R-ALL-9657579)
  bound to a RAF complex variant that contains ARAF.
- ARAF knockout zeroes that drug complex. Under `divide`, a zero inhibitor
  de-represses the step about 10-fold: 0.46 × 10 = 4.6.

This is the RAF drug-complex class (specs/032, `DS_DRUG_MODE=inert`, not
adopted on its own: experimental −5).

**Reading:**
- The RAF knockout→UP class has two layers.
  1. Carrier drain, which `closed` removes.
  2. Drug-complex de-repression, which specs/032 addresses.
- Any arm that fixes RAF knockouts without addressing both is fixing them
  through 14-3-3 drain.
- Shared-intermediate pools (`sctl`) also break HDR, 6 overexpression cases
  and 49 tuning cases. That is not traced here.

### What is live

- **`closed` is a curator gain on both canonical and shared builds**
  (held-out +81 and +61, p ≤ 1.5e-6, with 12 or 11 pathways gaining).
- Its experimental cost is exactly the drain-dependent cases (ERCC1, ERCC4,
  CCND1), one fix short of the gate.
- **A candidate next registration:** `closed` together with
  `DS_DRUG_MODE=inert`, on build S.
  - The traces predict the RAF knockouts would then read DOWN through the
    real path.
  - specs/032's drug rule lost RAF cases on its own. That was before carrier
    loops were conserved, so it is a new arm, not a re-run.
  - It must be registered before it runs.

## P8. Arms `kdrug`, `kcd` and `scd` (amendment 3; solver `cf31ed5`)

`cf31ed5` adds only the amendment to `eea8f8c`; the solver is that of
`3c64f21`. Every arm is scored against `kctl`
(`bench/analysis/arm_compare.py`; case diffs by
`~/deltasignal-catalogs/analysis/049/exp_diff.py`).

**Structural checks hold.** Coverage is identical on both axes (18,481 /
5,030 curator, 845 experimental). `Drugs: rule inert` holds 15,086 curator
and 6,648 experimental node-solves. The carrier line reads `closed` on `kcd`
and `scd`, with 3,444 shared pool nodes on `scd`.

| arm | curator held-out | curator tuning | experimental | RAF exp | worst pathway |
|---|---|---|---|---|---|
| `kdrug` | +8 (8/0), p 0.008 | +3 | −5 (5/10) | 21/49 | RAF −7 (exp) |
| `kcd` | **+89** (110/21), p 9e-16 | +101 | **+1** (14/13), p 1 | 27/49 | TGF-β −9 (cur) |
| `scd` | +69 (118/49), p 1e-7 | +42 | −6 (10/16) | 21/49 | RAF −7 (exp) |

**Predictions:**
- **D-P1 partly met.** `scd` fixes 6 of the 9 knockouts (ARAF, BRAF at
  three readouts each), short of 7. They are real paths: converged, the
  14-3-3 dimer (R-HSA-2028645) reads 1.0×, and p-MAPK reads 0.47–0.63×.
  RAF1 knockout still reads 3.45× UP.
- **D-P2 missed.** `scd` breaks 13 RAS overexpression cases (HRAS, NRAS,
  KRAS): the readouts fall from 13–100× to 0.97–1.14×.
- **D-P3 missed.** RAF 21/49; experimental −6.
- **D-P4 missed.** `kcd` fixes 9 RAF knockouts (ARAF, BRAF, NRAS), but by
  collapse: they read 1e-6 to 4e-7, and HRAS/NRAS up read 0 in the same
  arm (10 broken). The RAF loop sits at baseline unperturbed and falls to
  zero under any perturbation; HRAS up does not converge. RAF net −1.
- **D-P5 met.** `kdrug` is −5 experimental, RAF 28 → 21: the specs/032
  overexpression collapse replicated on today's networks.
- **D-P6 met.** `kcd` +89 against `kclosed` +81; `scd` +69 against
  `sclosed` +61.
- **D-P7 met.** ERCC1, ERCC4 (HDR) and CCND1 (Mitotic G1) are broken in
  `kcd` and `scd` as in `kclosed`.

**Outside RAF, `kcd` gains +2 experimental:** EGFR and PTK6 knockouts in
ERBB2 (0.87 → 0.69, 0.88 → 0.62, the carrier rule) and two CDK4 knockouts in
Mitotic G1 (19.7 → 0.84, the drug rule), against BRCA2 fixed and ERCC1,
ERCC4, CCND1 broken.

**Against the registered rule:**
- `kdrug` and `scd` fail (experimental −5, −6).
- **`kcd` passes every registered gate:** experimental +1 > 0, held-out +89,
  no pathway losing more than 9, structural checks hold. It passes as a
  combination; neither rule passes alone (`kclosed` −1, `kdrug` −5).
- **The margin is one case at p = 1,** and its RAF component is a loop
  collapse, wrong-for-the-wrong-reason in both directions (net −1). The gain
  does not depend on it: outside RAF the net is +2.

**Decision (Adam, 2026-10-11): `kcd` is not adopted,** although it passes
the registered gate. This overrides the registered rule, and is recorded as
an override:
- the experimental margin is one case at p = 1;
- RAF's predictions would come from a collapsed loop, wrong-for-the-wrong-
  reason in both directions;
- the curator gain is almost all `closed`'s (+81 alone), so the open
  question is its drain-dependent experimental losses, not the drug rule.
No default changes.

**What the drug rule does to RAS overexpression** (both builds): with the
drug complexes held at baseline, a RAS-up signal no longer reaches p-MAPK;
the readout reads ~1× on S and 0 on K. That is the specs/032 mechanism,
traced there (the RAF phosphorylation step's catalyst is made downstream of
it). Holding drugs fixes the knockouts and removes the overexpression route,
so `inert` cannot pass on RAF without that loop fixed first.

## P9. Do the cases `closed` loses have a curated route? (Neo4j, Release 97)

`closed` loses ERCC1 and ERCC4 knockouts in HDR, CCND1 overexpression in
Mitotic G1 and ARAF knockout in RAF; each was right under `off` only through
a carrier drained to zero. Before treating them as losses, each was checked
against the pathway as curated (`~/deltasignal-catalogs/analysis/049/
curated_route.py`: reaction → reaction through an output used as an input,
catalyst or regulator; `ROUTE_STRICT=1` lets only the gene's own forms climb
into complexes, so a released subunit cannot jump into any later complex that
contains it).

**ERCC1, ERCC4 knockout → R-HSA-5686483 (Holliday junction resolution).**
- ERCC1:XPF is used only in the SSA branch (R-HSA-9980006, 5686657,
  9980021, 5686663). The readout is in the HRR branch.
- The only curated link is machinery released by "Completion of SSA"
  (R-HSA-5686663) and used again in HRR: free RAD51 (strict), or MRN, ATM,
  BRCA1, EXO1/DNA2/BLM, KAT5, BRIP1 (lenient, via the D-loop complex).
- That is a recycled-component hop, the same class as the RPA drain the
  network used under `off`, not a causal route. SSA and HRR are
  alternative branches.

**CCND1 overexpression → R-HSA-157443 (CCNE1, made by R-HSA-8961846).**
- Strict: **no curated route.**
- "Cyclin D:CDK4/6 phosphorylates RB1 and prevents RB1 binding to E2F1"
  (R-HSA-69227) outputs only monophosphorylated RB1; no free E2F1. The
  release is in the reaction's name, not its outputs.
- The lenient route runs through p-RBL2 into the DREAM complex, which
  inhibits E2F1 transcription: the opposite sign.
- Under `off` the case was right through the drained PTK6 product
  R-HSA-8848007, which is not on any curated route to CCNE1.

**ARAF knockout:** traced in P7; the curated route exists, and needs the drug
complex and the MAPK pools handled (P8, `scd`).

**Reading.**
- Three of `closed`'s four experimental losses were right under `off` for
  a reason Reactome does not contain. Under `closed` the network says what
  the curation says: no route, so no change.
- The experimental answer for those cases comes from outside the pathway's
  mechanism (cross-pathway effects, or biology not curated here). Getting
  them right from this pathway alone would be an artefact.
- This argues for `closed` on faithfulness, the ground on which specs/048
  adopted the bind fix. It is post hoc, so it is recorded as an argument for
  Adam's decision, not as a re-scoring.

## P10. Every case `kclosed` changes, classified by its curated route

Adam (2026-10-11): before counting a gain, check that getting the case right
is a good idea. The curators' answers may come from reasoning the pathway does
not contain, and experimental answers may need other pathways.

`~/deltasignal-catalogs/analysis/049/route_class.py` (Neo4j, Release 97)
classifies each case `kclosed` fixes or breaks against `kctl`:
- **causal:** a curated route in which no hop is a part released unchanged,
  and only the gene's own forms climb into complexes;
- **recycled:** a route only if released-unchanged parts (RPA, RAD51, MRN)
  may hop;
- **lenient:** a route only if any entity may climb into a complex that
  contains it;
- **none:** no route in the pathway.

It reproduces the P9 hand traces (ERCC1/ERCC4 recycled, CCND1 lenient, ARAF
causal). `ROUTE_UNION=1` searches the union of the 92 catalog pathways.

**Curator held-out (102 fixed, 21 broken):**

| expected | route | fixed | broken |
|---|---|---|---|
| no change | causal | 2 | 0 |
| no change | recycled / lenient / none | **77** | 0 |
| up or down | causal | 14 | 4 |
| up or down | recycled / lenient / none | 9 | 17 |

- **77 of the 102 fixes are "no change" where Reactome has no causal
  route.** Under `off`, a drained carrier moved these readouts; under
  `closed`, the network says what the curation says.
- **Where a causal route exists, `closed` is +14 / −4,** and every fixed case's
  route sign agrees with the expected direction.
- **Where the curator expects a move without a causal route, `closed` is
  +9 / −17.** Those answers come from reasoning beyond the pathway, and
  `off` matched them through drain.
- So the +81 is +77 faithful no-change, +12 on causal routes (with the 2
  causal no-change fixes) and −8 on curator calls the pathway does not
  support.

Tuning (99 fixed, 10 broken) has the same shape: 89 no-change fixes without
a causal route, causal +9 / −3 (the 3 are the RAF knockouts).

**The 4 causal held-out breaks are one mechanism: PTPN11 knockout in FGFR1–4,**
whose RAS:GTP readout (R-HSA-109783) goes from about 0 to 6.67×. Traced in
FGFR1:
- the free GRB2 (R-HSA-74686, a closed free form) keeps a **depletion**
  in-edge from the CBL feedback complex "Activated FGFR1:p-FRS2:p-PTPN11:
  p-CBL:GRB2" (R-HSA-5654267);
- PTPN11 knockout zeroes that complex, so GRB2 is de-repressed to the 10×
  depletion cap;
- GRB2:SOS1 → "Activated FGFR1:p-FRS:GRB2:SOS1" → RAS nucleotide exchange
  (R-HSA-5654392) then read 10×.
- Under `off` the same GRB2 drained to 3e-6×, so the readout read DOWN.

`closed` drops only the activating loop edge (amendment 2, fix 3). A
depletion into a closed free form by a complex that binds it is the
consumption that conservation says is negligible, so it is a candidate for
the same rule. That is a new rule and needs its own registration.

**Experimental (3 fixed, 4 broken), within the pathway and across all 92:**

| case | change | own pathway | all 92 pathways |
|---|---|---|---|
| EGFR KO (ERBB2) | fixed | causal, sign agrees | causal |
| ARAF KO (RAF) | broken | causal, sign agrees | causal |
| BRCA2 KO (HDR), expected no change | fixed | recycled | causal, either sign |
| PTK6 KO (ERBB2) | fixed | recycled | causal, either sign |
| ERCC1, ERCC4 KO (HDR) | broken | recycled | causal, either sign |
| CCND1 up (Mitotic G1) | broken | lenient | lenient |

- Across the catalog, ERCC1/ERCC4 do reach the readout through other
  pathways, but through cycles with an inhibition, so the route gives either
  sign and says nothing about direction.
- CCND1 → CCNE1 has no causal route anywhere in the catalog (P9: R-HSA-69227
  releases no E2F1).

**Reading.**
- `closed`'s curator gain is mostly the network ceasing to move readouts
  that Reactome does not connect to the gene. Its losses are mostly cases
  `off` got right only through drain.
- On cases with a causal route, `closed` is +15 / −5 on both axes together,
  and the losses are two mechanisms: RAF (drug complex and MAPK rail, P7/P8)
  and the GRB2 depletion above.
- This is recorded as evidence for a faithfulness decision on `closed`, not
  as a re-scoring. The registered numbers stand.

## P11. Arm `ksup` (amendment 4; solver `dfa8a04`, build K)

`DS_CARRIER_LOOP_MODE=closed_supply`. Structural lines: `Carrier loops: rule
closed_supply`, 3,902 / 1,194 free-form-solves changed (curator /
experimental), as `kclosed`; coverage identical (18,481 held-out, 845
experimental).

| comparison | curator held-out | tuning | experimental |
|---|---|---|---|
| `ksup` vs `kctl` | **+124** (143/19, p 1e-24) | +97 | **−2** (3/5, p 0.73) |
| `ksup` vs `kclosed` | +43 (49/6, p 2e-9) | +9 | −1 (0/1) |

Worst pathway against `kctl`: TGF-β receptor −9 (the SMAD7 cases `kclosed`
already lost); RAF −3. Against `kclosed` every pathway is ≥ 0 net.

**Predictions.**
- **S-P1 met.** Against `kclosed`, every changed case is in the 11 listed
  pathways (FGFR1–4, ERBB2, DNA Damage Bypass, Mitotic G2-G2/M on both axes).
- **S-P2 met.** The four PTPN11 KO → RAS:GTP cases are fixed (6.67× → 0.667×),
  as are the eight PTPN11 KO/UP no-change cases on the FGFR complex
  readouts (10× / 0.1× → 1.0×).
- **S-P3 met, by more than predicted** (+124 against ≥ +85). The extra +39 is
  no-change cases in FGFR1–4, ERBB2 and Mitotic G2-G2/M that the depletion
  had been moving.
- **S-P4 missed.** Experimental −2, one worse than `kclosed`: PTPN12 UP in
  ERBB2 (readout p-6Y-ERBB2 heterodimers, expected DOWN, 0.668× → 1.0×).
  Reactome has no route from PTPN12 to that readout in the pathway or in the
  union of the 92 (class `none`); the old DOWN came through the RNF41
  depletion now dropped.

**Gates (amendment 4).**
1. Registered gate: experimental −2. **Fails.**
2. Route-qualified gate: PTPN12 UP is excluded (no route anywhere), leaving
   experimental 3 fixed / 4 broken = −1. CCND1 is `lenient`, not `none`, so it
   counts. **Fails.**

**Decision: `closed_supply` not adopted**, under both registered gates.

**The six new curator losses, traced.**
- **PTPN11 KO → PP2A:p-SPRY2:GRB2 (FGFR1–4, expected UP, 5× → 0.5×).**
  `kclosed` was right for a wrong reason: the complex read 5× because free
  GRB2 sat at the 10× de-repression cap. With GRB2 at 1.0×, the complex
  follows PP2A:p-SPRY2, whose two producers are SRC phosphorylation
  (R-HSA-1295609, 1.0×) and "PTPN11 dephosphorylates SPRY2"
  (R-HSA-1549564, 0 under the knockout), averaged to 0.5×. The expected UP is
  substrate accumulation: the readout is the input that PTPN11's reaction
  consumes, and the network has no depletion edge for that consumption.
- **AURKA KO/UP → p-S95-PHLDA1 (Mitotic G2-G2/M, tuning).** AURKA has two
  nodes: the pinned root (`a290ed4c`) and the carrier loop's free form
  (`233bb72f`, closed). Under `closed` the only edge left into the free form
  was a depletion from SCF-FBXL7:AURKA, which carried the knockout through two
  sign inversions. Under `closed_supply` the free form reads baseline: **the
  gene's own free form is cut off from its perturbation.**
- **BTC KO (ERBB2, 3 tuning cases)**: small magnitude moves (0.82× → 0.90×)
  across the DOWN threshold; not traced further.
- RCHY1 KO/UP in DNA Damage Bypass: `lenient` routes only.

**The rule's own gap.** The second AURKA point is general. "A closed free
form reads its supply" must mean the gene's level when the carrier's gene is
perturbed, but the rule holds the free form at baseline unless it is pinned,
and the bench pins the root copy, a different node. On build K, 11 of 204
closed free forms have a root copy of the same entity: GRB2 in FGFR1–4,
SMAD7 and SMURF2 in TGF-β, AURKA, ubiquitin in ERBB2, ANK1 in L1CAM, CALR
and TAPBP in MHC class I. The ten TGF-β SMAD7 losses that `kclosed` already
had (SMAD7 KO/UP reading 1.0×) are this mechanism. It is the subject of
amendment 5.
