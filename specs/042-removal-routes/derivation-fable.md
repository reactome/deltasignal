# specs/042 — derivation (Fable, blind; 2026-10-02)

## 0. Verification: curated topology or ours? (Neo4j R97 vs build 20260928-1110_06ccb63)

For every entity in an M1 trace: its in-pathway roles in Neo4j, and the build's
uuid copies of that stId with their edges (`scratchpad/fable042/topo_check.py`).

| cases | Reactome | our build | class |
|---|---|---|---|
| WNT APC, AMER1 KO/80x (8) | CTNNB1 [cytosol] 448839 is **one** entity: input of 195304 (binds destruction complex) and 201669 (nuclear import), output of 201685 (release). Exit 2130282 (proteasome). No synthesis. | **3 copies**: 9aa158df root (in 0 → 195304 only); 094a7dc9 (← 201685 → 201669); 2c39329c dissociation dead end. AMER1 [cytosol] split the same way. | **generator artefact**: root split + dissociation sink |
| TP53 CDKN2A (10; M2/M1) | TP53 Tetramer 3209194 is one entity: output of 6804762, 3215310 (USP7), 6804996 (p14ARF); input of 21 reactions. | **4 copies**: b223452c acts; the USP7 and p14ARF products each dead-end in a dissociation. | **generator artefact** |
| HDR BLM KO → resolution (3) | HJ complex 5686228: input of 5686410 (BLM dissolution), 5693584 (cleavage), 9853389 (FIGNL1). | one copy, faithful | **curated as-is**; mixed with M4 (BLM in the EXO1,DNA2:BLM,WRN set pool rails the chain to 0) |
| HDR RTEL1 KO → D-loop resolution (2) | extended D-loop 5686104: input of 5686440 (MUS81), 5693539 (ligation), 5693589 (SDSA, RTEL1). | one copy | **curated as-is** |
| HDR PALB2 KO/80x, RAD51 KO → SSA (3) | resected end 5685162: input of 5693561 (RAD51/BRCA2 → HR), 5693580 (RAD52 → SSA), 5684882 (CHEK1, **returns** it via 5684887). | one copy | **curated as-is** |
| S-Phase CDKN1B 80x → p-FZR1/p-RB1 (1) | CCNA:CDK2 141608: input of 187934 (p27 binding), 187949 (CAK → active 187952, the readout's catalyst), 174164 (WEE1). | one copy | **curated as-is** |
| G1 RB1 KO → CCNA1 (1) | free RB1 68642: input of 9018017 (binds E2F1/2/3:DP set → 68644, **no consumer**) and 69227 (CDK4/6, "prevents RB1 binding to E2F"). CCNA1 is expressed from E2F1:DP 68653, a set *member*. | faithful; the E2F1 root assembles into set and member | **curated as-is**: dead-end sequestration |
| G1 CDKN1B KO → CCNE1 (1; M9b/M1) | p27 binds CCNE:CDK2 68374 → 68376, a dead end; active CDKs are curated **negative regulators** of 9018017. | faithful; readout also matches p27-containing complexes | **curated as-is** + readout mapping |

WNT and TP53/CDKN2A are copy splits. The other six groups are faithful and
lack a semantics Reactome never states as an edge: one state with two fates,
or a binding that parks its partner.

## 1. What the verified cases call for

**(a) Generator: rejoin.** A produced entity must land on the node the pathway
consumes: one node per stId per pathway wherever the stId is both produced and
consumed; dissociation edges point at the existing copy. Census in §5. This
re-closes cycles the uuids had cut (034 §11); (b) is what makes that safe.

**(b) Rejoin alone fails; specs/039 alone fails.** Rejoined, free CTNNB1's only
producer is the release: under APC KO the loop has zero supply and iterates to
0, still inverted. As a 039 **closed** pool (free ⇄ bound, φ₀ = 0.1): free =
1/(0.1·0 + 0.9) = **1.11, NORMAL**. The fate APC controls is degradation, which
a closed pool cannot see. The missing quantity is the **total**: supply over
removal.

**(c) Solver:** fates of one species compete; a parked partner is removed.

## 2. The rule

Species X = one reference entity (UniProt/ChEBI); carriers = nodes whose
containment leaves include X. Transitions = reactions consuming a carrier as
*input* and producing a carrier (transfers). Exits = reactions consuming a
carrier and producing none (Ub, cofactors ignored). A transition whose carrier
is **released unchanged** (035's test, followed along X-carriers only, curated
outputs only, not derived dissociation edges) is a *return*, not a fate.

- **R1 total.** T_X = s_X / Σ_i w_i·u_i over removal routes i (maximal
  transfer chains from the free state to an exit), bounded by the 011 clamp
  [0.1, 10]. u_i = product of step drives with X at baseline (039 amendment
  2); a return branch on the route multiplies the committed step by
  k·u_c/Σ_k u_k. w_i = 1/n (equal shares, 039's branched-pool assumption);
  n = 1 is exactly MP-BioPath's "remover ⊣ X". The free state (root and
  rejoined return copies) reads s_X·T; downstream reads ordinary propagation.
- **R2 competition.** At a state with n ≥ 2 **transforming** transitions,
  transition i reads inflow × n·u_i/Σ_j u_j. Returns take no share; the
  perturbed fate is capped at n-fold, never divided.
- **R3 sequestration** (arm B). A transition into a carrier with **no
  consumer** is an exit for every non-cofactor species it parks, drive = the
  partners' fold.

Baseline exact (all u = 1 ⇒ T = 1, shares 1). Off byte-identical:
`DS_REMOVAL_MODE`, `DS_FATE_COMPETITION`, `DS_REMOVAL_SEQUESTER`;
`LNG_REJOIN_COPIES` on the generator.

## 3. Not the failed rules

| | naive depletion (−14pp) | general consumption (NORM-F1 0.39→0.19) | here |
|---|---|---|---|
| lowers | input of every catalysed reaction | every consumed input, from co-inputs | **sibling fates / the total** |
| perturbed route | ×u·(1/u): flux conserved → NORMAL | same | **n·u/Σu ≥ 1**, or **T·u** |
| identity | stId (p-S "consumed") | stId | reference entity: a transfer is not a removal |
| enzyme binding | fires | fires | a return; excluded |

Both failures lowered the signal's carrier and conserved flux through the
perturbed step. R1/R2 never divide the perturbed fate. The Pi/Ub depletion
edges (phosphatase removes the phospho species; ligase the protein) are R1's
one-step case and are load-bearing (−90 to remove). 019 bridged severed sinks
(nothing bridged here); 035 held loop inputs (nothing held).

## 4. Composition

- **039**: closed pool = R1 without exits (T = s); π unchanged, scaled by T.
  R2 at a pool state is 039's branched split, with drives.
- **011**: same `H_dep` operator and clamp; `edge_type=removal`, deduplicated.
- **022/040**: drivers carrying X or self-fed from X are excluded; no closure.
- **033/038**: an exit consuming a set reads its OR-mean over virtual copies.

## 5. Detection and census (`scratchpad/fable042/census.py`, canonical build)

The generator (it alone has `referenceEntity`) ships `removals.csv` (species,
free-state uuids, exit, route steps, siblings) and the rejoin map; the solver
computes u_i, T and shares per solve.

| | count | pathways |
|---|---|---|
| root split (WNT shape) | 499 stIds | 50 |
| dead-end produced copies (TP53 shape) | 614 stIds / 923 copies | 48 |
| dissociation-sink copies | 9,738 | — |
| protein removal reactions (R1) | **97** on 243 species | 40 |
| branch states, any siblings | 2,064 | 88 |
| … ≥ 2 transforming fates (R2, lower bound) | 1,302 | 86 |
| dead-end bindings (R3) / partners | 1,257 / 1,670 | 89 |

WNT removals: CTNNB1 (2130282), DVL1/2/3, AXIN1/2, APC2; S-Phase: p27/p21
(187574). The R2 bound used an stId return test that wrongly excluded free RB1
(a derived dissociation) and the resected end (RPA re-entering through another
species); the species-aware count is a pre-arm deliverable.

## 6. Worked cases (pin 80 / 0; UP ≥ 1.15, DOWN < 0.85)

**WNT.** Rejoined free CTNNB1; route bind (APC/AMER1 in the complex) → CK1α →
GSK3 → {β-TrCP → 2130282 | PP2A → release 201685, drive u_r = WNT receptor};
n = 1, committed fraction 2u_ub/(u_ub + u_r). APC or AMER1 KO: u_bind = 0 ⇒
T → **10** (cap); free, nuclear and TCF targets **UP** (today 2e-13). APC 80x:
1/80 → **0.1 DOWN** (today 100). WNT1 KO: fraction 2 ⇒ **0.5 DOWN** (today ≈ 0,
still right). WNT 80x: 2/81 ⇒ **10 UP**.

**HDR.** Resected end, fates HR and SSA (CHEK1 is a return), n = 2. PALB2 80x
(HR step 5693620): SSA = 2/81 = **0.025 DOWN**, HR 1.98 UP. PALB2 or RAD51 KO:
SSA **2.0 UP**. RTEL1 KO at 5686104 (n = 3): MUS81 and ligation **1.5 UP**.
BLM KO at 5686228 (n = 2): resolution **2.0 UP**, once the set-pool product
stops zeroing the chain.

**p27 (S-Phase, R2).** CCNA:CDK2, fates p27-binding, CAK, WEE1 (n = 3). CDKN1B
80x: CAK fate 3/(80+1+1) = **0.037** ⇒ p-FZR1/p-RB1 **DOWN**. KO (p21 copy
stays 1): 3/(0.5+2) = **1.2 UP**.

**RB1 KO (G1, R3).** For species E2F1 the dead-end RB1 binding is an exit,
n = 1, u = RB1 fold ÷ active-CDK inhibitors: RB1 KO ⇒ T_E2F1 = **10**; the E2F1
root, E2F1:DP and CCNA1 expression **UP**. CDKN1B KO raises active CDK2 (R2),
lowering that exit's drive ⇒ CCNE1 **UP**.

## 7. What can go wrong

- **Rejoining re-closes artefact loops** the uuids had resolved (034 §11);
  safe only with R1. Arm rejoin+R1 against rejoin alone; record the
  cyclic-node delta per pathway (020).
- **R2 is the blanket risk**: 1,302+ states in 86 pathways; its guard is the
  transform test. If NORM-F1 falls by more than half the gross gain, refuse it
  as general consumption was refused.
- **R3**: 1,257 dead-end bindings, every partner a brake. Own arm.
- **Equal shares** understate a dominant exit (three ligases: one KO = 1.5x);
  direction holds for any shares. The n-fold cap forbids magnitude claims.
- **Before adoption**: species-aware census; pathways moved and distinct
  readouts; McNemar held-out, both axes; trace one gain and one loss per
  pathway with |Δ| ≥ 10; the specs/041 expected-biology review.
