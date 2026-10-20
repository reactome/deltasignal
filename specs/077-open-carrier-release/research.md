# specs/077 research: results of arms ocf and ocfd

Solver `882bef1`, run 2026-10-19 with `scripts/run_arm.sh`, `DS_CARRIER_LOOP_OPEN=fed`.

| arm | build | against | results |
|---|---|---|---|
| ocf | canonical `20261017-1246_5e663d9` | `results/bf4b688/rst` | `results/882bef1/ocf` |
| ocfd | faithful `20261018-0815_ddbd828_dsub` (+ specs/071, 072 rules) | `results/bf4b688/rstd` | `results/882bef1/ocfd` |

The logs report 13,160 dropped release edge-solves (ocf, curator) and 20,478 (ocfd, curator).

## ocf: canonical (gated)

`holdout_report.py --set-aside RAF_MAP_kinase_cascade`:

| split | net | fixed / broken | p |
|---|---|---|---|
| curator held-out | **+35** | 82 / 47 | 0.0026 |
| curator tuning | +16 | 36 / 20 | 0.044 |
| experimental (RAF set aside) | **−5** | 8 / 13 | 0.38 |

Curator held-out in release goes from 89.91% to 90.10%, and macro-F1 from 0.8706 to 0.8741.

By pathway (curator): MET +20, EGFR +16, IL-2 family +8, HDR +8, Interferon
α/β 0 (2 / 2), **intrinsic apoptosis −26** (4 / 30), and nothing else below −2.
Experimental: Mitotic G1 −5, ERBB2 −2, S phase −1, Cell Cycle Checkpoints +2,
HDR +1. Converged curator cases fall from 23,284 to 23,244 (65 lost, 25 gained).

**The registered gate fails** on its third condition: intrinsic apoptosis is
−26 (below −10). Held-out passes (+35 ≥ 0), and experimental passes at its limit (−5 ≥ −5).

## ocfd: faithful build

- Curator net **+217** (312 fixed, 95 broken).
- DSB Repair **+182** (182 / 0): every MMEJ case recovered.
- MET +20, IL-2 +16 and EGFR +16 as on canonical.
- Losses: intrinsic apoptosis −26, RUNX2 −21 (CBFB and ESR1 knockouts) and TP53 −5.
- Experimental is −5, the same cases as `ocf`.

## Predictions

1. In `ocfd`, DSB Repair gains at least 150: **holds**, +182.
2. Interferon α/β moves by at most 10 cases: **holds**. It is 0 net (2 / 2) on both builds, against −168 under `all`.
3. In `ocf`, fewer than 3 pathways lose more than 10 cases: **holds**. One does (intrinsic apoptosis).

## The apoptosis loss (traced)

Every changed case is a BAD, YWHAE, BCL2, BAX or AKT1 perturbation. All of them now read **exactly 1.0×** at
every readout (caspase-3, caspase-7, apoptosome and tBID forms). Under `rst` these perturbations moved the readouts.
All 30 broken cases are BAD knockout and overexpression, YWHAE knockout and overexpression, and AKT1 overexpression
(6 readouts each). The 4 fixed cases are BAX knockout and overexpression, where the curators say NORMAL.
BCL2 is wrong under both: inverted under `rst`, NORMAL under `ocf`.

The cut edge is "BAD displaces tBID from BCL-2 sequestration" (R-HSA-139897):
BAD + tBID:BCL-2 → BAD:BCL-2 + MYS-BID(61-195), where MYS-BID(61-195) is tBID at the mitochondrion. The generator lists MYS-BID
(R-HSA-114266) as the free form of an open loop:

- Its outside producer is "Translocation of tBID to mitochondria" (R-HSA-139920).
- Its loop releases are R-HSA-139897 and "tBID activates BAK protein" (R-HSA-139895).

`fed` drops both releases, so BAD no longer reaches tBID. BAD, its 14-3-3 sequesterer
YWHAE, and AKT1, which phosphorylates BAD, are all cut off.

This release is not the RBBP8 kind. The conservation argument ("bound forms are
negligible, so free = total") is wrong for a sequestration, where the bound share is
the point, and a displacer changes the free share. The curators agree:
BAD knockout reads DOWN and overexpression UP.

Structurally, the releasing reaction has a second non-cofactor input, the
displacer (BAD). Among the edges `fed` drops (`count_displacers.py`, this directory):

| build | single-input release | release with a displacer |
|---|---|---|
| canonical | 682 | 23 |
| faithful | 969 | 49 (24 in R-HSA-170834) |

The MMEJ release of RBBP8 (R-HSA-5687675) and the other DSB releases are single-input.
In intrinsic apoptosis only R-HSA-139897 has a displacer. The other releases are
caspase-3 and caspase-7 from DIABLO:XIAP, the apoptosome from DIABLO:XIAP, and
tBID from BAK and BAX. They are single-input, and their loss is not traced here: the
changed cases go through BAD.

## Reading

- `fed` does what it was built for: DSB +182 on the faithful build, and
  Interferon α/β untouched. The canonical held-out net is +35 (p 0.0026), most of it MET,
  EGFR and IL-2.
- It fails the gate on one traced mechanism: a release driven by a displacer is a real
  supply change, and the rule removed it.
- **Not adopted as registered.** A candidate amendment (not run): exempt releases
  whose reaction has a second non-cofactor input. It is post hoc, motivated by this
  arm, and needs Adam's go-ahead and its own pre-registered arm.

## Amendment 1: arms ocs and ocsd (`fed_single`)

Solver `24e3c1f`, run 2026-10-19 with `scripts/run_arm.sh`. Post hoc: the rule was
written after `ocf` failed, so the held-out figure below is selected on held-out.

| arm | build | against | results |
|---|---|---|---|
| ocs | canonical `20261017-1246_5e663d9` | `results/bf4b688/rst` | `results/24e3c1f/ocs` |
| ocsd | faithful `20261018-0815_ddbd828_dsub` (+ specs/071, 072 rules) | `results/bf4b688/rstd` | `results/24e3c1f/ocsd` |

Dropped release edge-solves: 12,692 (ocs curator) against 13,160 under `fed`.

### ocs: canonical (gated)

| split | net | fixed / broken | p |
|---|---|---|---|
| curator held-out | **+65** | 84 / 19 | < 1e-4 |
| curator tuning | +20 | 34 / 14 | 0.0055 |
| experimental (RAF set aside; RAF does not move) | **−5** | 6 / 11 | 0.33 |

| | accuracy | macro-F1 |
|---|---|---|
| held-out, in release | 89.91% → **90.26%** | 0.8706 → **0.8757** |
| tuning, in release | 79.14% → 79.53% | 0.7789 → 0.7822 |
| experimental | 71.38% → 70.79% | 0.6320 → 0.6355 |

By pathway (curator, 17 moved): MET +20 (30 / 10; CBL −10), EGFR +16, IL-2 family +8,
HDR +8, intrinsic apoptosis **+4** (against −26 under `fed`), worst VEGF −2.
Against `ocf`, `ocs` is +34 (42 / 8): apoptosis +30 (every BAD, YWHAE and AKT1 case
returns), Mitotic G1 +6, EGFR 0 (6 / 6), Cell Cycle Checkpoints −2 (MDM2 knockout).

**The registered gate passes**: held-out +65 ≥ 0; experimental −5 ≥ −5, at its limit;
no pathway below −10 (worst −2).

### ocsd: faithful build

- Curator net **+272** (311 / 39): held-out +257 (277 / 20), tuning +15.
  Curator all pathways 85.01% → 86.14%, macro-F1 0.8210 → 0.8333.
- DSB Repair **+182** (182 / 0), MET +20, IL-2 +16, EGFR +16.
- RUNX2 **0** (−21 under `fed`): its CBFB and ESR1 losses were displacer releases too.
- Losses: TP53 −5 (AKT1, AKT2, CHEK2 knockouts), VEGF −2.
- Experimental −5, the same 17 cases as `ocs`.

### Predictions

1. `ocs`: intrinsic apoptosis within ±5 of `rst`: **holds**, +4.
2. `ocsd`: DSB Repair gains at least 150: **holds**, +182.
3. `ocs`: curator held-out at least +20: **holds**, +65.
4. `ocsd`: RUNX2 −21 not predicted to change: **wrong**, it returns to 0. Not traced
   beyond that; the exemption keeps whatever release cut it.

### The experimental loss (partly traced)

The 11 broken cases: Mitotic G1 6, S phase 2, ERBB2 3 (ERBB3 knockout, PTK6 knockout
and overexpression). 6 fixed: Mitotic G1 MYC and RB1 knockouts (3), S phase CCNE1
overexpression, ERBB2 PTK6 overexpression, HDR CHEK1 knockout.

- Mitotic G1 knockouts of CCND1, CCNE1 and CDK4 at CCNA1/CCNA2 (4 cases) are a
  threshold shift, not a lost route: `rst` reads 0.83–0.84× (DOWN), `ocs` 0.865–0.875×
  (NORMAL). The same readouts read 14× (UP) under `fed`, so the amendment improves on `fed` here.
- CCNE1 overexpression at CCNA1/CCNA2 (2 cases) is a lost route: 20× under `rst`,
  1.08× under `ocs` and `fed`. `fed_single` drops two releases in Mitotic G1: the SCF(Skp2)
  ligase released after p27/p21 degradation (R-HSA-187574) and Cyclin E:p-T160-CDK2
  released after phosphorylating RB1 (R-HSA-188390, single non-cofactor input). The free
  kinase keeps CAK phosphorylation (R-HSA-188350) as its outside producer. Why that
  producer does not carry CCNE1 overexpression is traced below.
- S phase and ERBB2 losses are traced below.

### Reading

- `fed_single` passes the registered gate, with two caveats. It is post hoc. Its
  experimental net is −5, exactly the limit, with 2 cases a lost route that is not
  fully traced.
- On the faithful build it closes most of the DSB Repair gap and removes the RUNX2 loss.
- Adoption is Adam's decision.

### CCNE1 overexpression in Mitotic G1 (traced)

`g1trace.jl` and `g1c.jl` (this directory) solve R-HSA-453279 on the canonical build
with both CCNE1 roots pinned at 80, as the benchmark does. They reproduce the arms
exactly: CCNA1 reads 20.06× under `off` and 1.0807× under `fed_single`.

- **Under `off`, the 20× is the artefact the rule removes.** CCNE1 overexpression
  fills Cyclin E/A:CDK2:p-p27/p21:SCF(Skp2):CKS1B (R-HSA-187565). Degrading the p27/p21
  (R-HSA-187574) releases the SCF(Skp2):Cks1 ligase (R-HSA-187547), and that release is
  read as new ligase: free SCF reads **88.9×**. The ligase then multiplies back into the
  binding step (R-HSA-187552). Cyclin A:p-T160-CDK2 (R-HSA-187952) reads 100×, and so do
  the CCNA-containing forms that make up most of the CCNA1 readout set. The curated answer (UP) was
  reached through a cyclin E overexpression making 89× of a ubiquitin ligase.
- **Under `fed_single`, the free ligase reads 0.273×, not baseline.** Its two remaining
  producers are Cks1 association (R-HSA-187545) and the p130 ubiquitination release
  (R-HSA-1363331, a separate cycle the table does not list). Cks1 association reads
  0.273× because free CKS1B (R-HSA-187533) is drained. It has two depletion in-edges, from the
  substrate complex R-HSA-187565 (7.0× under `fed_single`) and from the p130
  complex. These are complexes that contain CKS1B depleting free CKS1B. Under `off`, CKS1B is
  drained to the 0.1× floor too, but the release masked it.
- The edges are depletion edges the generator derived, not a curated reaction: R-HSA-187565 is formed
  from SCF(Skp2):Cks1, not from free CKS1B.

**The generator's ubiquitin pass emits them** (`_emit_substrate_depletion_edges`,
LNG `5e663d9`; the specs/057 class, confirmed against Neo4j Release 97). The pass treats
any input leaf whose gene is also in the output as the "substrate", and adds a depletion edge from the
reaction's catalyst to every network copy of that leaf. Two Mitotic G1
ubiquitinations hit CKS1B:

| reaction | catalyst (= its own input) | output | CKS1B in the output |
|---|---|---|---|
| R-HSA-187575 "Ubiquitination of phospho-p27/p21" | R-HSA-187565 (cyclin E/A:CDK2:p27/p21:SCF(Skp2):CKS1B) | the same complex + 3xubiquitin | unchanged, inside the Ub complex |
| R-HSA-1363331 "Ubiquitination of p130 (RBL2) by SCF (Skp2)" | R-HSA-1363326 (p-RBL2:SCF(Skp2):Cks1) | CUL1:SKP1:SKP2:CKS1B + Ub-p-RBL2 | unchanged, the ligase released |

The real substrates are p27/p21 and p130. CKS1B, SKP2, CUL1 and SKP1 are the ligase,
carried through unchanged. The rules of specs/057 (LNG PR #113, closed, not adopted)
treat the two reactions differently:
- `LNG_UB_SUBSTRATE_CHANGED`: a leaf whose own stId is in the output is carried. That
  drops both drains.
- With `LNG_UB_SUBSTRATE_UBX` (amendment 1): such a leaf is still the substrate if it sits
  inside an output complex that has Ub as a direct component. That drops the p130 drain
  (R-HSA-1363331) and keeps the p27/p21 one (R-HSA-187575), where CKS1B sits inside
  the Ub complex.

specs/057 measured these rules before `fed_single` existed, when released carriers
still read their release as supply and so masked the drain. Its result
(held-out −39 / −15, not adopted) does not say how they behave together with `fed_single`.

So `fed_single` is right about the ligase. Conservation says free SCF ≈ total SCF, which
cyclin E does not set. The lost case comes from removing a wrong route that happened
to give the right answer. What it exposes is a second defect: a derived drain on a
subunit by the complex that contains it. `closed_supply` drops such drains into
closed free forms (specs/049 amendment 4). `fed_single` drops releases into open free
forms but not the drains on their supply. No rule is proposed here.

### The other five experimental losses (traced)

`trace5.jl` (this directory; one line per case: pathway, gene, direction, readout, …,
pinned uuids, readout uuids) reproduces every arm value.

| case | readout | `off` | `fed_single` | mechanism |
|---|---|---|---|---|
| S phase CCNE1 KO | R-HSA-187568 (p27/p21 complex with SCF + 3xUb) | 0.014× | 1.004× | SCF release |
| S phase CCNE1 KO | R-HSA-187920 (p-FZR1, p-RB1) | 0.346× | 1.226× | SCF release |
| ERBB2 ERBB3 KO | R-HSA-1306951 (NRG:ERBB3:ERBB2:GRB7) | 0.807× | 0.915× | threshold shift |
| ERBB2 PTK6 KO | R-HSA-109783 (RAS:GTP) | 0.681× | 1.000× | PTK6 release |
| ERBB2 PTK6 OE | R-HSA-109783 (RAS:GTP) | 55.5× | 1.000× | PTK6 release |

- **S phase** is the Mitotic G1 mechanism again. Under `off`, the CCNE1 knockout drops the
  SCF(Skp2) ligase's release, so free ligase reads 0.39× and the p27/p21 degradation
  chain falls to 0.019×. Under `fed_single` the free ligase reads 1.155×. That is above
  baseline, plausibly because the knockout empties the complex that drains CKS1B (the
  generator's ubiquitin pass, above). DOWN was reached by a cyclin E knockout setting the
  ligase's level.
- **PTK6** (2 cases): "ERBB2 promotes PTK6 autophosphorylation" (R-HSA-8848005) releases
  the p-ERBB2 heterodimers that PTK6 bound. Under `off` that release made PTK6 set the
  receptor's level: heterodimers at 67× under PTK6 overexpression, SHC1 complexes 100×,
  RAS:GTP 55×. In Reactome PTK6 is downstream of ERBB2. This is the same class as
  PLC-gamma setting the PDGF receptor's level, which specs/068 removed (adopted, +18).
  The measured direction may well be real biology, but this network has no curated
  route for it.
- **ERBB3 KO** is a threshold shift (0.81× → 0.92×). The releases from PTK6 and PLCG1
  steps into the ERBB2 heterodimers no longer carry the knockout.

**Summary of all 11 experimental losses:**
- 6 are the release artefact the rule removes: Mitotic G1 2, S phase 2, ERBB2 PTK6 2.
  Each one's right answer came through a released carrier read as new supply.
- 5 are threshold shifts across the DOWN boundary: Mitotic G1 4, ERBB2 ERBB3 1. None of them loses a route.
- None is a correct curated route that the rule cut.

## Adoption (2026-10-19)

Adam, 2026-10-19: "ok" to adopting `fed_single` after the eleven experimental losses were
traced. `DS_CARRIER_LOOP_OPEN=fed_single` is the code default; `off` restores the earlier
solver. Next, as specs/078: the generator's ubiquitin pass (the CKS1B drain) on top of
these defaults.
