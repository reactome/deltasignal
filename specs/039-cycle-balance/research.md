# specs/039 — Interconversion cycles as conserved pools

> **Status (2026-09-28): measured, NOT adopted.** Pre-registered, amended three
> times before any arm, measured on build pools039d: see "Result" and
> "Amendment 4 result" below. The rule ships default-off (`DS_CYCLE_MODE=off`).
> Open: problem3.md (where a pool ends).

**Status:** design, 2026-09-27. Awaiting Adam's decision on the fixed
assumption (φ₀) and on where detection lives, before pre-registration.

## Why

specs/034 §11–12:
- our loops are mostly curated interconversion cycles (A ⇄ B);
- the solver multiplies around them (gain 1, the knife-edge, the RAF collapse);
- MP-BioPath's curators cut such loops by hand and added "reverse enzyme ⊣
  modified form" inhibitions.

The catalog has **259 curated interconversion pairs in 55 pathways**.
- 93 have a catalyst in both directions, 150 in one, 16 in neither.
- Phospho 89, ubiquitin 59, GTP/GDP 32, binding/dissociation 31, methyl 28.

## Method: a blind double derivation (agreed with Adam)

The same problem statement (`problem.md`) went to two models independently:
- this session (Opus 5.5): `derivation-opus.md`, written before the other was
  read;
- a Fable 5.1 subagent: summarised below.

## The two derivations agree on the rule

A steady-state rate balance on a conserved pool, A + B = total:

    r = (u + ε) / (v + ε)       D = φ₀·r + (1 − φ₀)
    B_fold = s · r / D           A_fold = s / D
    F_fold = R_fold = J = s / ((1 − φ₀)/u + φ₀/v)      (the cycle flux)

- u is the forward reaction's drive, evaluated with the existing AND/divide
  semantics with A held at baseline.
- v is the same for the reverse reaction, with B held at baseline.
- s is the supply fold of the cycled protein: its pin, else the OR-mean of its
  external producers, else 1.
- φ₀ is the baseline fraction of the pool in the modified form.

**Multi-state cycles:** J = s / Σⱼ (φⱼ / uⱼ), and Aᵢ = J / uᵢ.

**Properties:**
- exactly baseline when unperturbed;
- no gain-1 loop and no all-zero root;
- bounded: B ≤ s/φ₀;
- the direction of every enzyme perturbation is correct for any φ₀ in (0, 1).

## They differ on φ₀, and the difference matters

| φ₀ | 80x through a 4-tier cascade (RAS → RAF → MEK → ERK) | source |
|---|---|---|
| 0.5 (Opus, uninformative) | 1.98, 1.33, 1.14, **1.07: ERK reads NORMAL** | the compression flagged in both derivations |
| **0.1** (Fable; modified forms a minority at rest: RAS ~5–10% GTP-loaded, basal phospho-occupancy ≤10–20%) | stays UP for about 20 tiers | knockouts are exact 0 at any depth, for any φ₀ |

**Proposed: φ₀ = 0.1, fixed a priori** on the biological justification.
- Orientation: the product of the catalysed, modifying direction is B.
- φ₀ = 0.5 would be a declared sensitivity arm, never used to choose.

Fable notes that the φ₀ → 0 limit (B = s·u/v, A = s) is exactly MP-BioPath's
hand cut plus its added inhibitions, with no parameter at all. That holds only
where the reverse reaction is curated (24% of the added inhibitions, specs/034
§12b).

## Added by the Fable derivation

- **Reaction nodes carry the flux J**, so co-outputs (ADP, Pi) read flux.
- **Drop the E_f ⊣ A depletion edge inside a cycle**; the rule already lowers A.
- **Autocatalysis** (B catalyses F): bistable. Take the nonzero root that is
  continuous from baseline.
- **Detection** by shared member leaves (forms of one reference entity).
  Degradation must not be mistaken for a reverse step.
- **Failure modes:**
  - φ₀ near 1 (constitutively modified sites);
  - unmodified-form readouts read NORMAL under an enzyme perturbation at
    φ₀ = 0.1 (1.11);
  - a set-member knockout gives u = (N−1)/N, so B reads DOWN only for N ≤ 6;
  - overlapping cycles need the Markov-chain form;
  - a pinned member fixes the other form by conservation.

## Fixtures (φ₀ = 0.1, s = 1)

| perturbation | A | B | F = R |
|---|---|---|---|
| none | 1 | 1 | 1 |
| E_f KO | 1.111 | 0 | 0 |
| E_f 80x | 0.1124 | 8.989 | 8.989 |
| E_b KO | 0 | 10 | 0 |
| E_b 80x | 1.1096 | 0.01387 | 1.1096 |
| gene 80x (KO) | 80 (0) | 80 (0) | 80 (0) |
| both enzymes 80x | 1 | 1 | 80 |
| two tiers, kinase₁ 80x | A₂ = 0.556 | B₂ = 4.997 | 4.997 |

## Open decisions (Adam)

1. φ₀ = 0.1 as the fixed assumption.
2. Detection in the generator, which ships `cycles.csv` (forward, reverse, A, B
   per virtual-reaction pair), rather than in the solver.

## Decisions (Adam, 2026-09-27)

1. **φ₀ = 0.1** for every cycle, fixed a priori. φ₀ = 0.5 runs only as a declared
   sensitivity arm and is never used to choose.
2. **Detection in the generator.** It detects each protein's pool once (its
   forms, by shared reference entity, and the curated reactions converting
   between them) and ships it with the network. DeltaSignal solves the pool on
   every run.
3. **General form: π = πP** (Adam's framing).
   - Each pool is a Markov chain: the states are the protein's forms, and each
     transition rate is its baseline rate × its drive.
   - The stationary distribution gives the split between forms.
   - The two-state closed form above is the special case.
   - For a branched pool, the baseline rates need one extra stated assumption:
     an equal split of flux between exits, or detailed balance. It is to be
     chosen and stated in the pre-registration.

The explanation for readers is in `docs/MODEL.md` §4.

## Pre-registration (2026-09-27, before any code or arm)

**Pool shapes** in Neo4j across the 92 pathways, as connected groups of forms
linked by curated forward/reverse pairs: **213 pools**, of which 178 have 2
forms, 35 are chains or branches of 3–9 forms, and none is a simple ring.

**Detection (generator, `pools.csv`):**
- A pool is found at **node (uuid) level**. For a curated pair (F: A → B,
  R: B → A in one pathway), the node-level loop A_u → F → B_u → R → A_u must
  exist: the same A node feeds F and receives R's output.
- Stid-level pairs whose occurrences the uuids separated are NOT pooled, so
  artefact merges are not re-created.
- Forms joined by such loops form one pool.
- The file lists, per pool, each form node (and whether it is the base form)
  and each transition (from form, to form, reaction node).

**Base form** (the least modified):
- the fewest modified residues summed over the form's leaves (Neo4j);
- then the fewest components;
- remaining ties go to the smaller stId, and are counted and reported.

**Baseline** (the one assumption, **φ₀ = 0.1**, Adam): detailed balance, with
each step away from the base form holding (φ₀ / (1 − φ₀)) = 1/9 of its
predecessor's share.
- 2 forms: 0.9 / 0.1.
- A 3-form chain: 0.890 / 0.099 / 0.011.
- The baseline rates follow from this; branched pools need no extra assumption.

**Solve (`DS_CYCLE_MODE=off` default | `balance`, with `DS_CYCLE_PHI=0.1`):**
- Each transition's drive is its reaction evaluated with the existing semantics
  and its source form held at baseline, divided by baseline.
- Rates are baseline rate × drive.
- π is solved from πQ = 0, Σπ = 1.
- The supply s is the OR-mean of the forms' producers outside the pool, else 1;
  a pinned form sets the pool.
- Form value = baseline × s × π / π₀. Each transition's reaction node carries
  its flux, relative to baseline.
- This is applied inside the component iteration every sweep, and is
  order-invariant.
- The catalyst ⊣ source-form depletion edge inside a pool is not applied (the
  balance already contains it).
- Autocatalytic pools take the root continuous from baseline.
- A pool the solve cannot handle falls back to iteration, and is counted.

**Arms:**
- **Build:** one build from the generator with `pools.csv`; its networks are
  otherwise identical to the canonical build.
- `ctrl039`: `DS_CYCLE_MODE=off`.
- `bal01`: `balance`, φ₀ = 0.1. **Primary.**
- `bal05`: `balance`, φ₀ = 0.5. A declared sensitivity arm, never used to choose.

**Adopt `bal01` only if all hold** (against `ctrl039`, same build):
- curator held-out net > +15, p < 0.05;
- experimental net ≥ +15, p < 0.05, excluding RAF (RAF reported separately);
- no pathway, tuning or held-out, loses more than 10;
- gains span at least 2 pathways and 5 perturbations.

**Predictions:**
- convergence improves (fewer non-converged solves);
- RAF's collapse stops (RAF reported);
- the knockout and overexpression of reversal enzymes (phosphatases, GAPs)
  gain cases.

## Amendment 1 (2026-09-27, after the adversarial review, before any arm)

The first build (`20260927-2257_c210f21_pools039`, 102 pools) was reviewed
before any arm ran. The review found that detection did not implement what this
pre-registration describes, plus defects in the solve. No arm was run on that
build, and it is not used.

**Detection, as found.** Only 24 of 102 pools were pure interconversions of one
protein.
- 61 were enzyme binding cycles (E + S → E:S → E + P).
- 17 had forms sharing no reference entity.
- 11 used proteasomal degradation as the "reverse" step.

The query asked only for "A → B and B → A". It never checked "forms of one
protein", which this pre-registration's own text requires.

**Detection, amended (generator `1358d03`):**
- Each transition reaction's only non-small-molecule input is the source form,
  and its only non-small-molecule output is the other form. Ubiquitin counts as
  a co-substrate. This excludes binding, release, degradation, and reactions
  taking two forms.
- The two forms share a reference entity of an EWAS leaf.
- A reaction node that would be a transition more than once is dropped, with
  its loops. Its node would be written by two fluxes, which gave a
  label-dependent value and no convergence (R-HSA-5654736, residual 0.79).
  They are counted.

**Base form, amended.**
- The first build decided a donor rule first, then components, then residues.
  The donor rule was not pre-registered, and the order was not the one
  pre-registered. 30 of 78 exchange pools were oriented backwards; all 24 pure
  pools were oriented correctly.
- Now: residues first (as pre-registered). Then the donor rule, added because
  RAS:GDP and RAS:GTP have no residue difference: the product of the direction
  that consumes a group donor (ATP, GTP, SAM, acetyl-CoA, NAD+, ubiquitin) is
  the modified form. Then components.
- Undecided pools fall back to the smaller stId, and are counted.

**Solve, amended:**
- **Supply** comes only from producers outside the pool's cyclic component.
  A producer inside it can be fed by the pool itself (11 pools).
- **Regulariser:** rates are k·(u + 1e-9), not k·u + 1e-12. With every drive at
  zero, the split now tends to π₀. Before, the constant chose it: a GEF+GAP
  double knockout read CDC42:GTP at 5.0x.
- **Depletion and inhibitor edges into a form from outside the pool** are
  applied as a fold on that form. Before, pool management silently dropped them
  (13 edges). The depleted share is removed, not redistributed.
  - The pool's own catalyst ⊣ source-form edges stay unapplied, as
    pre-registered.
- **Transition nodes** carry their flux as drive × the fold of the source form.
  This equals s·u·π_i/π₀_i when no modifier applies.
- **`DS_SCC_METHOD` `pool*` / `minimize`** never run the sweep that updates
  pools, so under them every pool is reported unmanaged. Before, they were
  miscounted as solved.

**Stated limitations (not fixed):**
- Influence scores do not see the pool rule.
- The autocatalysis rule is not implemented. No detected pool has a form
  catalysing its own pool's transition; this is checked on the build and
  reported.
- Several pools in one component are updated in a fixed order each sweep. The
  result is order-invariant only at a unique fixed point.

The arms, gates and predictions above are unchanged. They run on a rebuild
(`pools039b`) from generator `1358d03`. Its pool counts are recorded below
before any arm.

**Rebuild `20260927-2322_1358d03_pools039b` (recorded before any arm):**
- **25 pools**, all two-form (50 forms, 204 transition-reaction copies), in 12
  pathways. RHO GTPase cycle (R-HSA-9012999) has 12 of them; RAF/MAP (R-HSA-5673001) has 1.
- 0 orientation ties; 0 multi-use reactions dropped.
- At the stId level the amended query finds **55 form pairs in 19 pathways**.
  The node-level check keeps the 25 whose loop exists in the network.
- **Why so few:** the loose query found 203 distinct form pairs. 164 fail the
  one-protein-interconversion test, almost all because they are the enzyme
  step E ⇄ E:S (USP9X ⇄ Ub-SMAD4:USP9X, PPM1A ⇄ p-SMAD:PPM1A, …).
  - Reactome usually curates a modification cycle in several steps: bind,
    modify, release.
  - The substrate cycle SMAD4 → … → Ub-SMAD4 → … → SMAD4 runs through
    enzyme:substrate complexes, and a one-reaction A → B → A pool cannot see it.
  - These **multi-step modification cycles** are the larger class. This rule
    does not cover them, and treating them needs a design of its own: the
    enzyme is a second conserved species. It is not attempted here.
- **Consequence for the gates:** with 25 pools in 12 pathways, the
  held-out > +15 and experimental ≥ +15 gates may be out of reach. The gates are
  not lowered. A result under them is reported as it is.

## Amendment 2 (2026-09-28): multi-step cycles, pre-registered before any code or arm

**Why.** The Fable review of amendment 1 traced RAF/MAP. The GAP reaction is
curated as bind then hydrolyse (RAS:GTP + GAP → RAS:GTP:GAP → RAS:GDP + GAP), so
the single-step pool never sees the GAP. NF1 KO then reads RAS:GTP 1.000 (pERK
1.0) where `off` read UP by accident. That contradicts prediction 3 in the
pathway the pre-registration names. No arm was run on amendment 1.

**Method.** A second blind double derivation from `problem2.md`:
`derivation2-opus.md` was committed (7778402) before `derivation2-fable.md` was
read. The derivations agree on the structure. The rule below takes Fable's
choice on each point of difference:

| point | Opus | Fable (adopted) | why |
|---|---|---|---|
| baseline share of an enzyme complex | r = 1/9 | **0** (β₀ → 0) | Only the limit reduces to the pre-registered two-form rule (kinase 80x: 8.99; bound 1/φ₀). A share of 1/9 compresses to 3.67 and moves the bound to 12.2. |
| states | forms carrying exactly the core proteins | least-bound form per **modification signature** | Signature = residues on R's leaves ∪ small molecules in the complex whose direct children include R's leaf. It attaches GTP to RAS, not to the GAP; R97 nests these (checked). |
| intrinsic-step weight | per pair | **per exit** from a state: ε_int = 1e-3 for a non-enzyme exit beside an enzyme-driven one; self-catalysis counts as non-enzyme | RAS's intrinsic hydrolysis competes with the GAP *path*, which a per-pair weight cannot see. |
| enzyme free form | loop left iterated | release-step producers dropped | Tested as its own arm (below): it resembles specs/035, which lost. |

**Rule (generator: detection).**
- **R-steps.** For each reference entity R: a reaction whose exactly one input
  and exactly one output contain R, both at stoichiometry 1, and which differ.
- **Pool candidate.** A strongly connected set of the R-graph with ≥ 2 forms and
  ≥ 2 signatures. A set with one signature is a *carrier* loop (an enzyme, a
  scaffold), not a pool.
- **States.** Per signature, the forms with the fewest non-small slots. The rest
  are intermediates. The base state follows amendment 1 (residues, donor,
  components).
- **Transitions.** Simple directed paths state → intermediates* → another state,
  ≤ 6 steps. Longer paths are dropped and counted. A path is **enzyme-driven**
  if a step has a catalyst or a joining input containing a protein other
  than R.
- **Node level.** As in amendment 1: every step exists as a reaction node
  between form nodes, and each uuid path is a parallel transition.
- **Shared nodes.** A node claimed by two pools is removed from both, recomputed
  and counted.
- **Carriers.** The enzyme's free-form node, and the release-step nodes that
  produce it.

**Rule (solver).**
- **π₀.** Over states only, as pre-registered (ρ^m). Intermediates hold 0.
- **Baseline flow.** J₀ is the stationary flow of the random walk on states
  whose exits are weighted 1 (enzyme-driven) or ε_int. Then k_p = J₀_p / π₀_from,
  so π₀ Q₀ = 0 on any graph. For two forms this is the amendment 1 rule, and the
  fixtures are unchanged.
- **Drive.** u_p = Π over steps of (the step with its source at baseline) /
  baseline. Rate k_p·(u_p + 1e-9); π Q = 0.
- **Values.**
  - A state reads baseline × s × modifier × π/π₀.
  - A path's flux fold is u_p · (state fold of its source).
  - A step reaction node or intermediate reads Σ_p J₀_p·flux_p / Σ_p J₀_p over the
    paths through it.
- **Carrier (enzyme) free form.** It reads baseline × (the mean fold of its
  producers other than release steps, else 1) × modifier. Pinned carriers keep
  their pin.
- **Unchanged from amendment 1.** Supply, pins, the unapplied catalyst ⊣
  source-form edge, and fallbacks (counted).

**Stated approximations.**
- Enzyme in excess, not conserved: substrate 80x drives its complexes and flux
  80x.
- Sequestration of an enzyme between substrates is invisible.
- A product of drives over-multiplies when two steps of one path are perturbed.
- A stable complex that is not an enzyme intermediate (p-SMAD2,3:SMAD4, HRR's
  13-form repair machine) reads as flux: direction right, magnitude not.
- A set-member enzyme KO reads (N−1)/N.

**Expected fixtures** (Fable hand-trace on pools039b RAF; numbers asserted in
`test_cycles.jl`):

| perturbation | RAS:GTP |
|---|---|
| NF1 KO | 1.29 UP |
| all GAPs KO | 9.99 |
| SOS1 80x | 8.99 |
| KRAS 80x | 80 (all forms) |
| GEF KO | 0.001 |

Six-form ring:

| perturbation | S* |
|---|---|
| kinase KO | 0 |
| phosphatase KO | 10 |
| kinase 80x | 8.99 |

**Arms** (one new build; same pathway list):
- `ctrl039` (`off`);
- **`bal01`** (balance, φ₀ 0.1, carriers on). **Primary.**
- `bal01nc` (carriers off). Isolates the release-drop.
- `bal05` (φ₀ 0.5). Sensitivity, never used to choose.

**Gates and predictions:** those of the original pre-registration, unchanged.
Adopt `bal01` only if all its gates hold. If only `bal01nc` passes, adopt that
configuration instead, and record why.

**Coverage (R97, stId level, recorded before the build):** 78 pools in 37
pathways (32 multi-step, 46 single-step). Pools reach 7 of the 9 experimental
pathways. The node-level counts are recorded after the build, before any arm.

**File contract** (generator → solver):
- `pools.csv`: `pool_id,node_uuid,stable_id,role,is_base`, with role `state` or
  `intermediate`.
- `pool_transitions.csv`: one row per step, with
  `pool_id,path_id,step,source_uuid,target_uuid,reaction_uuid,reaction_stid,enzyme_driven`.
- `pool_carriers.csv`: `pool_id,carrier_uuid,release_reaction_uuid`.

The amendment 1 format stays readable, as one-step paths.

**Amendment 2, clarification (before any arm): parallel copies are per step.**
- Implementing "each uuid path is a parallel transition" literally produced
  paths as the Cartesian product of reaction copies per step:
  - RAF: 2,304 paths where the hand trace has 48;
  - HRR: 262,145 paths;
  - EGFR: 64,834 paths.
- **Rule.** A path is identified by its node sequence plus the curated
  reaction at each step. Each step lists its uuid copies of that one reaction.
  - A step's drive is the mean of its copies, since each carries an equal share
    at rest. With equal copy weights, the sum over the Cartesian product of the
    copy drives factorises into exactly this: the same rule, stated compactly.
  - A copy node reads its path's flux × (its drive ÷ the step's mean).
  - Distinct curated reactions between the same nodes stay distinct paths, so
    GEF and intrinsic exchange keep their per-exit weights.
- **Generator interpretations adopted** (Fable, generator commit message):
  - self-catalysis means a catalyst with no protein other than R;
  - direct-child signatures are read through set membership;
  - co-travelling proteins with identical pools are merged, not dropped;
  - states on no kept path are omitted and counted;
  - variant signatures are the parent's;
  - base-rule residues are summed over all leaves.

**Amendment 2 build `20260928-0037_ad5e06e_pools039c`** (generator `ad5e06e`;
recorded before any arm):
- Nodes 72,322 and edges 214,197: identical to pools039b, so only the pool
  tables differ.
- **41 pools in 19 pathways**: 103 states, 15 intermediates, 128 paths,
  582 step-copy rows, **7 multi-step pools**, 169 carrier rows.
- 0 orientation ties; 0 long-path or budget drops; 10 shared nodes dropped
  (R-HSA-170834, R-HSA-8848021).
- Autocatalysis: 177 source-form catalyst pairs, 0 product-form. So the
  pre-registered bistable root rule is not exercised.
- The generator's scratch run reported 20 pathways. Counting non-empty
  `pools.csv` in the build gives 19, and 19 is the figure recorded.

## Amendment 3 (2026-09-28): fixes from the pre-arm review of amendment 2

This is Fable's adversarial review, with an end-to-end RAF solve on the pool
tables. No balance arm had run. `ctrl039` ran at `173b36a` and is re-run at
the final commit, so every arm shares one solver commit.

**Solver fixes:**
1. **Supply reachability follows mass flow only** (activator edges).
   - Through a depletion edge (RAS:GTP:BRAP ⊣ mature RAS), the pool "reached"
     RAS's own maturation supply, so KRAS 80x and KO read NORMAL.
2. **A carrier that any pool writes is not a carrier.** Such cases are counted
   (`cycle_carrier_conflicts`).
   - In Insulin, pool2's carrier was pool1's intermediate. It was written twice
     and did not converge.
3. **An intermediate reads v / k_cat**, as the adopted derivation states: the
   path flux over the drive of the step that exits it.
   - The implementation had given it the whole path's flux. A blocked exit then
     read 0 where the complex accumulates: calcineurin KO put
     14-3-3:p-S99-BAD at 0.
   - This is an implementation error, not a rule change.
4. **A pinned intermediate or step copy makes its pool fall back to the
   iteration**, counted (`cycle_pools_pinned_fallback`), as pre-registered.
   - It had been silently overwritten: root_cycle pins RAS:GTP:GAP for
     SPRED1/2/3 KO.

**Generator fix (states must be core-only):**
- Curators annotate a receptor's residues inconsistently across complexes, so
  enzyme carrier loops became pools: FGFR4 ⇄ FGFR4:PLCG1 ⇄ FGFR4:p-4Y-PLCG1,
  and EGFR:ERBB2 with PLCG1/PTK6, the latter an experimental-axis pathway.
- **Rule:** the core is the set of proteins present in every form of the
  candidate. Only core-only forms can be states. The candidate is a pool only if
  its core-only forms carry ≥ 2 signatures. This is the core rule of
  `derivation2-opus.md` §1, combined with the adopted signature.

**Recorded in advance (not fixed):**
- **RAF's downstream collapse is untouched.** p-MEK, p-ERK and the p-MAPK set
  read about 0 for every perturbation in every mode. The MAPK cycle closes only
  through sets, which the stId-level detection does not see.
  - **Prediction 2 ("RAF's collapse stops") is expected to fail** for scored
    RAF readouts. It holds only for RAS:GTP and convergence.
- **`bal01` vs `bal01nc` differ only through UNPERTURBED enzymes.** Root pinning
  already pins the perturbed enzyme's carriers.
  - There the carrier rule is load-bearing. With carriers off, SOS1 80x reads
    RAS:GTP 1.00, because the release loop restores the gain-1 knife-edge.
- **`bal05` compresses by construction:** NF1 KO → 1.14 (NORMAL), SOS1 80x →
  1.98.
- **Expected fixture on the real RAF bundle after the fixes:** KRAS 80x → RAS
  forms UP. Checked on the rebuild, before the arms.

`test_cycles.jl` goes 174 → 188. The mutants for fixes 1 and 3 each go red.

**Amendment 3 build `20260928-0115_bdaed8e_pools039d`** (generator `bdaed8e`;
recorded before any balance arm):
- **40 pools in 18 pathways**: 99 states, 16 intermediates, 123 paths,
  567 step-copy rows, **7 multi-step pools**, 168 carrier rows.
- 1 carrier conflict excluded; 0 ties; 0 long-path or budget drops.
- The core rule dropped one pool: FGFR4, now a carrier loop.
- EGFR:ERBB2 (R-HSA-1227986) now has 2 core-only states (p-6Y / p-7Y
  ERBB2), with the PLCG1 and PTK6 complexes as intermediates. The other 39 pools
  are unchanged.
- **Kept, flagged:** R-HSA-177929 pool1's reverse step (183084, "CBL escapes
  CDC42-mediated inhibition") is a curation shortcut, not a dephosphorylation.
  It is not special-cased.
- **Uuids are re-minted per build**, so logic networks differ byte-wise from
  pools039c. All four arms, the control included, run on this build at one
  solver commit.

**Amendment 3, second end-to-end check (before any balance arm).** On pools039d
at aef2f34, fixes 1–6 were verified:
- Insulin converges.
- Calcineurin KO → 14-3-3:p-BAD 10 (UP).
- The SPRED pin falls back and is counted.
- FGFR4 is gone, and ERBB2 states move together under ERBB2 80x.

**One new defect: the pool suppressed its own supply.**
- The supply reaction "mature p21 RAS binds GDP" consumes mature RAS. The
  generator's derived depletion edge RAS:GTP:BRAP ⊣ mature RAS comes from a node
  the pool feeds.
- Every RAF number was therefore multiplied by s ≈ 0.48: SOS1 80x 4.38 (fixture
  8.99), all-GAP KO 4.61, KRAS 80x 11.8 (fixture 80), KRAS KO 0.63.
- **Rule:** under balance, a DERIVED depletion edge whose source the pool feeds
  (mass-flow reach) and whose target is in the activator ancestry of the pool's
  supply producers (inside the component) is not applied. It is counted
  (`cycle_supply_depletions_held`).
- Why: it is the pool's protein consuming its own precursor, which the
  conservation already holds. This is the same reasoning as the pre-registered
  "catalyst ⊣ source-form edge is not applied". Curated inhibitor edges are
  untouched.
- `test_cycles.jl` goes 188 → 194. The mutant goes red.

## Result (2026-09-28): not adopted

Build `20260928-0115_bdaed8e_pools039d` at solver `b119943`, all four arms
(`results/b119943/{ctrl039,bal01,bal01nc,bal05}`). Protocol `root_cycle` in
every arm. Balance arms: 932 curator / 212 experimental pool-solves (126 / 68
multi-step).

| arm vs ctrl039 | curator held-out | curator tuning | curator all | experimental |
|---|---|---|---|---|
| **bal01** (primary) | −5 (61/66, p 0.72) | −40 (p 8e-5) | −45 (p 0.003) | −8 (1/9, p 0.021) |
| bal01nc (carriers off) | −8 (p 0.53) | −63 | −71 | −13 (p 0.001) |
| bal05 (φ₀ 0.5) | −12 (p 0.31) | −49 | −61 | −7 (p 0.016) |

- **Gates:** bal01 fails the curator held-out and experimental gates, so it is
  **not adopted**. Nor is any other arm.
- **Convergence** (prediction 1) improved slightly: 1,671 → 1,679 of 1,725
  curator solves, 228 → 231 of 244 experimental.
- **The carrier rule is positive on its own:** bal01 vs bal01nc is +26 curator
  (29/3, p 3e-6) and +5 experimental (5/0).
- **Concentration** (bal01, all curator cases): 12 pathways moved.
  - Worst: Interferon α/β −56 (identical in all three arms), HDR −49.
  - Best: Interferon γ +34, ERBB2 +14, intrinsic apoptosis +13.
- The losses are traced below before any further change.

### The losses, traced (Fable, harness reproduces the arm values)

- **Interferon α/β, −56** (JAK1 80x and IFNAR2 80x read exactly 1.0; ctrl read
  100). **Rule as specified, meeting an artefact.**
  - The pool's first transition, "Activation of JAK kinases", is divide-inhibited
    by the SOCS-bound receptor 912681.
  - 912681 is built by an AND from the pool's own base state (909703).
  - So the inhibitor equals the base fold, and the flux base → pp is
    1.306 × (1/1.306) = 1.000, for any φ₀, supply or carriers.
  - ctrl was right only because its gain-1 loop railed to 100.
  - specs/022's damping does not fire: `containment.csv` lists 912681's
    components flat.
- **HDR, −63 / +14** (every 80x reads about 0; KOs read 100). **A detection
  defect plus the rule.**
  - The "pool" RPA ⇄ p-RPA runs a 5-step path through the resection machine.
  - Its exit step's co-inputs (p-RAD51, p-BRCA2:SEM1) depend on p-CHEK1, which
    is made from the pool's own intermediates.
  - With the top bounded, the only attractor is the zero root.
  - This is the "13-form repair machine" flagged in derivation2-fable §8.
- **Interferon γ, +34.** The same SOCS mechanism with the opposite sign, and a
  gain for the right reason in kind: conservation sets the total, and SOCS only
  redistributes within the pool.
- **Experimental −8** (1 won, 9 lost):
  - 6 of the losses are the HDR mechanism.
  - MYC ×2 (S phase) and CCNB1 (Prophase) read a base-state form at 1.11: the
    pre-stated "unmodified-form readouts read NORMAL at φ₀ = 0.1" failure mode.
  - The one gain is RAF NF1 KO (0.17 → 1.29, UP), the case this spec set out to
    fix.

## Amendment 4 (2026-09-28): POST HOC, designed after the arms above

**Not a rescue.** The pre-registered result stands: not adopted. This rule was
written after reading the traced losses, so its measurement carries the usual
forking-path caveat and is reported as exploratory.

**Rule.** In a pool's drives, a step input that the pool itself feeds by mass
flow is read at baseline, as the source state already is. "Feeds" means
reachable from the pool's states through activator edges, not through a pin.
Carriers and the source are excluded.
- This extends the pre-registered "catalyst ⊣ source-state edge is not
  applied" and amendment 3's "the pool does not suppress its own supply" to
  the pool's own transitions.
- The pool's protein acts on its own transitions only through its states, never
  a second time through complexes built from them.
- It covers Interferon α/β (the SOCS inhibitor built from the base state) and
  HDR (p-CHEK1 from the pool's own intermediates driving its exit).
- Expected effect on Interferon γ: pp 8 → 80, still UP.

**Arms and gates:** the same build, a new solver commit, `ctrl039` re-run at
that commit, `bal01a4` (balance, φ₀ 0.1, carriers on), and the same gates.
Because the rule is post hoc, it is adopted only if it passes on **both** axes.
The case-level concentration is reported with the pathways that motivated it
(Interferon α/β, HDR) excluded.

### Amendment 4 result (exploratory): not adopted

Arms at solver `393876f` on pools039d; ctrl039 is byte-identical to the
b119943 ctrl (0 of 23,511 curator and 0 of 845 experimental cases differ).

| bal01a4 vs | curator held-out | curator tuning | experimental |
|---|---|---|---|
| ctrl039 | **−41** (61/102, p 0.0016) | +11 (p 0.42) | −7 (3/10, p 0.09) |
| bal01 (pre-registered) | −36 | +51 | +1 |

- **HDR is fixed** (+51 against bal01): the self-fed exit was the mechanism.
- ~~Interferon α/β is unchanged at −56; the hold did not reach the SOCS
  inhibitor.~~ **Wrong, corrected 2026-09-28** (re-traced by the problem-3
  Fable derivation; confirmed on the case files).
  - The hold did reach the inhibitor: JAK1 80x and IFNAR2 80x were fixed.
  - SOCS1 KO and SOCS1 80x (28 + 28) were lost instead, because holding the
    SOCS-bound receptor at baseline silenced SOCS entirely.
  - The −56 is a different 56 cases.
- ~~Interferon γ lost its +34.~~ **Wrong, corrected.**
  - IFNG 80x and JAK1 80x stayed correct (+34 against ctrl).
  - SOCS1 KO and 80x were lost (−34) by the same silencing, so the net is 0.
- Also corrected: the bal01 γ gain was **the 100x cap, not conservation**.
  SOCS binds the p-JAK2 state (a 9% share at rest), which wants 870-fold and
  is capped at 100, so the inhibitor stops tracking. In α/β, the bound state
  holds 82% at rest, is not capped, and cancels exactly. The earlier "gain for
  the right reason in kind" is withdrawn.
- It fails on both axes. With the pre-registered arms and this one, no cycle
  rule variant is adopted. The code stays default-off (`DS_CYCLE_MODE=off`).

**What stands:**
- The rule fixes the case it was designed for, RAF NF1 KO (0.17 → 1.29 UP,
  experimental).
- It improves convergence.
- The carrier rule is positive against carriers-off.
- Its losses come from pools whose own protein regulates them through complexes
  outside the pool (SOCS, the HDR machine). The single-component pool
  abstraction does not yet capture that. The next design question is how far a
  pool's boundary extends, and it needs its own derivation, not another post-hoc
  arm.

## Amendment 5 (2026-09-28): the pool boundary. POST HOC; pre-registered before code

This is the third design after seeing results. It comes from a blind double
derivation of `problem3.md`: `derivation3-opus.md` was committed (bdf9ff2)
before `derivation3-fable.md` was read. It is judged on a stricter gate
(below).

**Agreement.** A pool's protein must not act on its own transitions a second
time through complexes built from it. Amendments 1, 3 and 4 were pieces of
this.

**Differences, and the choice (Fable's on each):**

| point | Opus | Fable (adopted) | why |
|---|---|---|---|
| forms that leave and do not return (A:SOCS) | sink states | stay outside the pool and read the live product | A sink needs a baseline leak-to-cycling ratio, a second φ. Opus's "no new parameter" claim was wrong. |
| a pool-fed input of a pool step | edge dropped | **re-evaluated** from its own inputs with the pool at baseline | Amendment 4's hold read A:SOCS at 1, silencing SOCS (SOCS1 KO/80x lost). Re-evaluation reads the SOCS fold. The double count is the regulator's R-dependence; the rest is genuine feedback. |
| machine paths (HDR) | reachability | **a multi-step path is an enzyme cycle only if every non-R entity a step consumes is output again by the same path** | The enzyme is regenerated. Stated at generation time. |

**Rule.**
- **Solver.** For a pool step's inputs that are mass-flow descendants of the
  pool (reach through activator edges, not through pins, excluding carriers):
  - re-evaluate them, and the ancestors inside the reach they depend on, from
    their own reactions;
  - do this with the pool's states, intermediates and copies at baseline;
  - use one pass, in a fixed order (reach BFS order over sorted ids).
  - The drive reads these re-evaluated values. This replaces amendment 4's hold.
- **Generator.** Drop multi-step paths that fail the regeneration test.
  Counted.

**Predictions** (Fable, on scratch copies of the three bundles, φ₀ 0.1):
- α/β JAK1 80x → UP (ISG20 100), SOCS1 KO pp 1.69, SOCS1 80x 0.03. Pathway
  total equals ctrl.
- γ IFNG 80x → UP, SOCS1 KO pp 1.05, SOCS1 80x 0.22. Pathway total equals
  bal01.
- HDR: the pool is dropped, total equals ctrl.
- Census: 3 of 9 multi-step paths fail regeneration (HDR, EGFR pool1, Insulin
  pool1). 6 of 40 pools have pool-fed step inputs (α/β, γ, HDR, Insulin pool2,
  Mitotic G2 CDK1, MET).
- **Predicted net +34 vs ctrl.**

**Checks before scoring:**
- Every pathway other than the seven named (α/β, γ, HDR, EGFR, Insulin,
  Mitotic G2, MET) must equal bal01 bit for bit. Otherwise the change leaked,
  and the arm is invalid.

**Gate:**
- curator held-out and experimental, as pre-registered;
- **and** net ≥ 0 on the cases outside the three motivating pathways (α/β, γ,
  HDR).

**Amendment 5, generator census (recorded before the build).** LNG
`6815dfd` (`feat/pools-regen`), counted on the pools039d networks:
- **37 pools; 4 multi-step** (BAD:14-3-3, RAS:GAP, SHC1:INSR, and PTK6 ⇄
  p-Y342-PTK6, which was newly visible once the shared-node rule no longer
  removed it); 92 states, 6 intermediates, 114 paths, 513 copy rows, 40
  carriers.
- **11 of 14 multi-step paths fail regeneration; 7 pools dropped.**
- This differs from the pre-registered census (3 of 9). The extra failures:
  - ERBB2 (R-HSA-1227986) ×4 and its copies in R-HSA-8848021 ×2. ERBB2's
    paths consume PLCG1 / PTK6 and release p-4Y-PLCG1 / p-Y342-PTK6: the partner
    leaves modified, so ERBB2 is the enzyme there, not a pooled substrate.
  - TGF-β (R-HSA-170834) ×2: the receptor complex and ZFYVE9 are not returned.
- The rule is applied as written. No "up to modification" clause is added,
  since that would be tuning to the census.
- **The byte-identity check's exempt set therefore widens.** Pathways allowed to
  differ from bal01: α/β, γ, HDR, EGFR, Insulin, Mitotic G2, MET, **ERBB2,
  PTK6 signalling (8848021), TGF-β (170834)**, and the pathways whose pool set
  changed in the new build. That last group is listed from the build, before
  scoring.

### Amendment 5 result (post hoc): curator gate passes, experimental gate fails

Build `20260928-0320_6815dfd_pools039e`. Arms: ctrl039 and bal01a5 at
`1897747` (solver = `151e940`), and bal01 at `b119943` on the same build for
the byte-identity check.

| bal01a5 vs | curator held-out | curator tuning | curator all | experimental |
|---|---|---|---|---|
| **ctrl039** | **+77** (87/10, p 1.8e-16) | −5 (p 0.18) | +72 (p 6e-13) | **−2** (1/3, p 0.62) |
| bal01 (same build) | +56 (56/0), α/β only | 0 | +56 | 0 |

- **Byte-identity check: passes.** Against bal01 on the same build, only
  Interferon α/β differs, so the solver change leaked nowhere. The rest of the
  change against ctrl comes from the generator's regeneration test. bal01 on this
  build is itself +21 held-out against ctrl (p 0.11); on pools039d it was −5.
- **Fable's pathway predictions: all three hold.** α/β equals ctrl, γ equals
  bal01 (+34), and HDR equals ctrl. Net against ctrl was predicted at +34 and is
  +72 (the extra is the regeneration drops elsewhere: PTK6 +19 …).
- **Concentration:**
  - 9 pathways moved: best γ +34, PTK6 +19, intrinsic apoptosis +13,
    chromatin +7; worst S phase −6, Prophase −1.
  - 28 perturbations.
  - Net outside the three motivating pathways: **+38**.
- **Experimental:** 1 won, RAF NF1 KO (0.17 → 1.29 UP). 3 lost, all the
  pre-stated "readout on the base state" failure mode at φ₀ = 0.1: MYC 80x ×2
  (S phase) reads 1.11, and CCNB1 80x (Prophase) reads 0.11.
- **Gates:**

  | gate | result |
  |---|---|
  | curator held-out > +15, p < 0.05 | **pass** |
  | experimental ≥ +15, p < 0.05 | **fail** (−2, n.s.) |
  | no pathway loses > 10 | pass |
  | ≥ 2 pathways and ≥ 5 perturbations | pass |
  | net outside the motivating pathways ≥ 0 | pass |

- **Verdict under the pre-registered rules: not adopted**, because the
  experimental gate fails. It is the first variant of this spec that is
  positive on the curator axis and neutral (not significantly negative) on the
  experimental one. Whether to adopt on those terms is Adam's decision. It is
  not decided here.

### The three amendment-5 experimental losses, traced (2026-09-28)

- **MYC 80x, S phase ×2 (pre-RC and p-FZR1/p-RB1 read 1.11; expected UP).**
  - The pool is CCNA:CDK2 ⇄ CCNA:p-Y15-CDK2 (WEE1 forward, CDC25A/B reverse).
    MYC raises CDC25A.
  - Y15 phosphorylation is **inhibitory**, so the unmodified base state is the
    *active* kinase. φ₀ = 0.1 gives it 90% of the pool at rest, so activation
    can raise it only to 1/(1 − φ₀) = 1.11.
  - φ₀'s justification is "activated forms are a minority at rest" (§Decisions).
    The implementation reads it as "modified forms are a minority". For an
    inhibitory modification the two are opposite.
  - This is Fable review finding 7 (amendment 2), now shown to cost cases.
- **CCNB1 80x, Prophase (p-lamin reads 0.11; expected UP).**
  - The lamin readout's ancestry reaches the lipin pools (LPIN1/2/3 ⇄
    p-S106-LPIN) eight edges up, through a set-pool catalyst.
  - CDK1 phosphorylates lipin, which inactivates it. Resting (active) lipin
    drops to 0.11, and the lamin reaction ANDs CDK1's direct action with the
    lipin-dependent chain.
  - The pool's direction (lipin inactivated) is plausible. The loss comes from
    AND composition downstream, not from the pool rule.
- Neither needs the pool boundary changed. The first is a question about what φ₀
  refers to, which is Adam's decision (φ₀ was his call).

## Amendment 6 (2026-09-28): φ₀ refers to the ACTIVE form. Pre-registered before code

**Decision (Adam, 2026-09-28):** φ₀ is the share of the **active** form at rest,
not of the modified form. This is the meaning its justification always had
("activated forms are a minority at rest").

**Rule (generator, orientation only; the solver is unchanged):**
- A pool's state is **active** if it acts downstream of the pool: at node level
  it has a `catalyst` edge, or a positive `regulator` edge, into a reaction
  that is **not** one of the pool's own steps.
- If exactly one state of a two-state pool is active, the other state is the
  base (the resting, inactive form), and φ₀ goes to the active one.
- For pools with more than two states: if exactly one state is active, the base
  is the non-active state furthest from it; ties fall back.
- Otherwise (no state active, or several) the current rule applies: residues,
  then donor, then components.
- Counted: pools oriented by activity, pools whose base **flips** relative to
  pools039e, and pools that fall back.

**Expected:**
- CCNA:CDK2 ⇄ CCNA:p-Y15-CDK2 flips (base = p-Y15). MYC 80x then drives active
  CDK2 up to 1/φ₀ → S-phase readouts UP (2 experimental cases).
- RAS keeps its orientation (base GDP).
- Lipin: active = unphosphorylated, which is currently the base, so it flips
  (base p-S106). CCNB1 80x → active lipin falls further; the lamin case stays
  lost (AND downstream).

**Census before the build:** every flipped pool listed by name and hand-checked
for plausibility, before any arm. If a flip is implausible (the "active" state
is not the biologically active one), that is recorded, not special-cased.

**Arms:** one build, ctrl039 and bal01a6 (balance, φ₀ 0.1, carriers on) at one
solver commit, plus bal01a5 re-run on the same build for byte-identity. Only
pools that flipped may change anything.

**Gates (to adopt `balance` as the default):**
- curator held-out > +15, p < 0.05, vs ctrl on the same build;
- experimental net ≥ 0, and not significantly negative;
- no pathway loses > 10;
- gains span ≥ 2 pathways and ≥ 5 perturbations.

The experimental gate is relaxed from the original "≥ +15, p < 0.05". The
reason is Adam's decision to hold amendment 5 "until the resting-form issue is
addressed". Adopting on that basis is **Adam's to confirm** before the default
changes.

**Amendment 6 census (LNG `8e527d7`, on pools039e; before any build or
arm).** As implemented, it flips **0 of 37** pools. 1 is oriented by activity
(PTK6, unchanged), and 36 fall back. On its own, the amendment is therefore
inert.
- **RAS does not flip,** as expected: RAS:GTP reaches RAF and PI3K as an
  *input*, not a catalyst.
- **CCNA:CDK2 does not flip, contrary to the expectation.**
  - Neither state catalyses anything. The acting kinase is CCNA:p-T160-CDK2,
    made from CCNA:CDK2 by CAK one reaction outside the pool.
  - Seeing it would need a rule extension ("the state the catalytic form is
    made from"). That is a new rule, not recorded as part of amendment 6.
- **Lipin: an implementation gap, and a clarification (not a rule change).**
  - Unphosphorylated LPIN1/2/3 (the current base) reach catalysis through
    `set_member` edges into the set-valued catalyst's pool node. Since
    specs/033, that is how a member of a set-valued catalyst is represented.
  - A `set_member` edge into a pool node that is itself a catalyst or positive
    regulator of a non-pool reaction therefore counts as "acts downstream".
  - Expected flip: lipin base → p-S106 (active = unphosphorylated). The lamin
    case is still expected lost: AND downstream.
- **Correction:** the generator report said "no lipin pool exists in this
  build". That is wrong: pools 2–4 of R-HSA-68875 are LPIN1/2/3 ⇄ p-S106.
