# specs/039 problem 3 — derivation 3 (Fable, blind)

2026-09-28, from `problem3.md`, `research.md`, `derivation2-fable.md`,
`cycles.jl`, build `pools039d` (arm tables under `results/`), and R97.
`derivation3-opus.md` was not read. Numbers below come from a scratch copy of
`src` with the rule implemented (no tracked file touched), run on scratch
copies of the three bundles. φ₀ = 0.1; nothing fitted.

## What the arms actually show (re-traced)

- **α/β −56 is two equal-sized effects, not one.** Under `bal01` JAK1/IFNAR2
  80x read 1.000 (28 + 28 lost); SOCS1 KO/80x passed. Under `bal01a4` JAK1 and
  IFNAR2 are **fixed (100, +56)** and SOCS1 is **lost (1.000, −56)**: the hold
  did reach 912681 and made SOCS inert. "Unchanged at −56" is a coincidence of
  counts; research.md's "not reached" is wrong.
- **γ's +34 under `bal01` was the 100x cap, not conservation.** SOCS binds the
  p-JAK2 state (9% share). Exact tracking would give pp = 1.000 as in α/β; the
  state wants 80 × 10.9 = 870, is capped at 100, so the inhibitor reads 100,
  H = 0.01, and pp = 80 × 0.01 × 10 = 8. With the base state (82% share) bound
  in α/β, nothing hits the cap and the cancel is exact. Same mechanism, one
  artefact.
- JAK1 80x reaches the α/β receptor pool only at s = 1.25 (a 13-copy OR mean
  over ligand variants); the readouts reach 100 by downstream STAT loops. The
  pool's job is only not to clamp that 1.25 to 1.000.

## The rule

**A drive is a rate per unit of source. Nothing the pool itself makes may
enter its own drives as the pool's protein.** For every input of a step
(activator, catalyst, inhibitor, depletion) that is a mass-flow descendant of
the pool (states, intermediates, step copies) and is not a state or carrier,
its **drive-side value is re-evaluated from its own inputs with the pool at
baseline** (recursively over pool-fed ancestors; fixed order; a pool-fed
ancestor loop is iterated with the states pinned). What remains is what other
species contribute; what is removed is the pool's own concentration entering
its rate a second time. It differs from amendment 4 in one word: *re-evaluate*,
not *hold*. Amendment 4 discarded the SOCS fold along with the receptor fold.

- A:SOCS = AND(A, SOCS) reads **SOCS fold** in the drive (and the live product
  A·SOCS as a node). u_forward = 1/SOCS_fold: SOCS KO de-represses, SOCS 80x
  suppresses, A does not cancel itself.
- p-BRCA2:SEM1 ← p-CHEK1 ← ATR-active intermediate reads **CHEK1 × BRCA2
  folds**, not the intermediate's v/k_cat.

## Q1. Forms that leave and do not return

**Outside the pool**, as now (they are not in the R-graph SCC). They read the
live product of their inputs (A:SOCS = 0.846 × 0.5 = 0.423 under SOCS1 KO), a
dissociation sink reads its flux, degradation reads its input. They do not
drain π: the pool's total is s, on the standing assumption that turnover is
uniform across states. Making them *sink states* would need the baseline
share of the pool's turnover that leaves through each sink; with π leaking,
supply must equal Σ_sinks k_e·u_e·π_e, and that k_e is a second φ we do not
have. Declined: not pre-registrable without a fitted number. Stated cost: SOCS
80x does not shrink the receptor pool, it only inhibits the step.

## Q2. Complexes built from a pool form that regulate the pool

The double count is exactly **the R-dependence** of the regulator; the genuine
feedback is **everything else** in it. For A + SOCS → A:SOCS the sequestered
fraction is [SOCS]/(K + [SOCS]), first order in A: more receptor does not
raise the *fraction* inhibited. `divide` with x = [A:SOCS] makes the forward
rate zeroth-order in A (the cancel), or, for a catalyst built from a state,
second-order (Mitotic G2's CDK1 loop). The rule gives α/β and γ the **same**
answer (u = 1/SOCS_fold); their opposite signs under amendment 4 were "SOCS
removed" landing on a pathway whose gain needed SOCS (α/β) and one whose gain
did not (γ, s = 80 through the ligand). The one class where R-dependence is
real, product-form autocatalysis (B catalyses A → B), is linearised by the
rule; the pre-registration already gave up bistability, and the build has 0
such pools (177 source-form, 0 product-form).

## Q3. Inputs made from intermediates (HDR)

Both. The drive rule alone reproduces amendment 4's HDR numbers (early
readouts 8.99, D-loop resolution 0.113, p-CHEK1 0.000 where curators say no
change). That residual is the intermediates reading v/k_cat of an exit that
is RAD51 loading onto a DNA-bound machine. So **detection**: a multi-step path
is an enzyme cycle only if **every non-R entity a step consumes is output
again by a step of the same path** (the enzyme comes back out). Otherwise the
partner is being transformed, R is a stoichiometric passenger, and the
intermediates are the partner's stages: the path is dropped and counted; a
state left without an exit drops the pool. In HDR the DSB machine 5684128
enters at step 1 and leaves as 5685317 (with RAD51, BRCA2), and p-BRCA2:SEM1
never returns. HDR then reverts to the iteration, i.e. ctrl (302), the best
of any arm there.

## Q4. What the generator can detect, and how

- **Pool-fed step inputs** (drive rule): graph reachability from the states
  through activator edges to a step's inputs; the solver already counts it
  (`cycle_self_fed_inputs_held`). Build census: **6 of 40 pools** — α/β
  (regulator 912681, 13 slots), γ (873821), HDR (32 inputs), Insulin pool2
  (input Insulin:p-6Y-INSR; regulator GRB10:INSR, a SOCS-like sequestration),
  Mitotic G2 pool2 (catalyst p-CDK1:CCNB1 set, a positive loop), MET (Ub).
- **Machine paths**: stId-level reaction inputs vs outputs along the path in
  Neo4j (`input`/`output` of each step). Census on the 9 multi-step paths:
  **3 fail** — HDR p2 (5684128, 9763139), EGFR pool1 p2 (182917), Insulin
  pool1 p1 (74674, the ligand). RAF, ERBB2, intrinsic apoptosis pass.
- **Sinks**: already excluded by SCC membership; nothing to add.

## Q5. The three cases, and SOCS1 KO (rule implemented, φ₀ = 0.1)

| case | states (base / modified) | readout | verdict |
|---|---|---|---|
| α/β JAK1 80x | 1.25 / 1.25 (s = 1.25) | ISG20 100 | **UP** (bal01 1.0, a4 100) |
| α/β IFNAR2 80x | 1.25 / 1.25 | 100 | **UP** |
| γ IFNG 80x | 80 / 80 | 100 | **UP** (= bal01) |
| HDR RAD51 80x | pool dropped → ctrl | 100 | **UP** (drive rule alone: 8.99 / 0.113) |
| α/β SOCS1 KO | 0.846 / **1.692** (u = 2, OR-mean of SOCS1/3 copies) | 100 | UP (= curators) |
| α/β SOCS1 80x | 1.216 / 0.030 | 0.0 | DOWN (= curators) |
| γ SOCS1 KO | 1.048 / 0.524 / 1.048 (ring, second step de-repressed: J barely moves) | 873814 50.5, ISGs 100 | UP by downstream loops; the pool itself reads NORMAL |
| γ SOCS1 80x | 0.218 / 8.82 / 0.218 | 0.001 | DOWN |
| HDR SOCS1 | SOCS1 is not in the network | — | NORMAL by construction |

Hand check: 3-ring, k = (0.407, 3.667, 3.667); u₁ = 2 → π_pp = 0.2727/1.773
= 0.1538, fold 1.692 exactly as solved.

**Pathway totals predicted before any arm** (curator, this build):
α/β **315 = ctrl** (+56 vs bal01), γ **263 = bal01** (+34 vs ctrl), HDR
**302 = ctrl** (+49 vs bal01). Net vs ctrl on the three: **+34**.

## Q6. What can go wrong, and the test that precedes the measurement

- **Byte-identity claim.** Every case in a pathway with no pool-fed step input
  and no machine path must equal `bal01` bit for bit: 12 of the 18 pool
  pathways and all 74 without pools. The diff is run **before** scoring; any
  other movement is an implementation bug, not a result.
- **The seven that may move, with direction:** α/β +, γ 0, HDR + (all vs
  bal01); Insulin (pool1 dropped; pool2: GRB10 KO → INSR pool modified state
  up, GRB10 80x → down, INSR overexpression no longer cancelled); EGFR
  (pool1 dropped, reverts to ctrl); Mitotic G2 (CDK1 catalyst linearised:
  CDK1/CCNB1 perturbations weaker, may lose cases where the switch was
  carrying them); MET (Ub input, expected inert).
- **Known costs the rule keeps:** SOCS 80x cannot deplete the pool (Q1); a
  ring's second step de-repressed reads NORMAL at the state (γ SOCS1 KO passes
  only through downstream loops); autocatalysis linearised.
- **Post hoc, third design.** Reported exploratory; adopt only if both axes
  pass the pre-registered gates **and** the net outside α/β, γ and HDR is ≥ 0,
  with gains spanning ≥ 2 pathways and 5 perturbations and distinct readouts
  counted. Arms: `ctrl039`, `bal01`, `bal01d3` at one solver commit on a
  rebuild whose only change is the machine-path drop (count them: 3 paths,
  3 pools).
