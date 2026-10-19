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
