# specs/125: a learnable model — design note

**Status:** design note, 2026-11-02. Nothing here is registered as an arm, and nothing is adopted. Each stage below
becomes its own pre-registered spec before it is built. Asked for by Adam, 2026-11-02, after a review of where the
project stands. In his words: "my non linear equations were something where you could learn parameters ... I had plans
of learning the parameters (shapes of the hill like curves) from data. and then the shapes of the curves would be
somehow associated with rates of reaction ideally"; the 100× ceiling "has been mentioned to be arbitrary"; and "I always
felt uncomfortable of 1 meaning normal amount ... in a lot of cells a particular gene would normally be off and then
in other cells be turned on".

## 1. Where we are

**What holds up.** The generator's networks: after ~120 specs they represent what curators recorded, including variant
complexes, sets, catalyst reuse, recycled copies and specialised complexes. On MP-BioPath's experimental cases our
propagator scores 610 / 849 on MP-BioPath's networks and 643 / 849 on ours (specs/047, 111): network construction is
where the gains came from.

**What the external data say.** On data the rules were not built on, the pathway logic adds little beyond "is there a
route from the target to the readout": HPN-DREAM 0.693 against reachability 0.667 (specs/112); decryptM and korkut below
simple baselines (specs/114); on CPPA's training drugs the per-drug correlation is far below a profile that ignores the
drug's identity (specs/121, numbers private). Direction is usually right near the target; beyond it the model calls too
many changes (specs/115), and its magnitudes have no units.

**What we drifted from.** The August 2025 design (Adam, "DeltaSignal — A Pathway Perturbation & Dynamics Engine")
specified per-edge Hill curves with learnable shape (steepness h, n; half-point K; inhibitor strength β), a stability
rule (Section 9: "keep maximal slopes < 1 in feedback SCCs"), and a time model (Section 4). To make network work
measurable we fixed one curve for every edge and spent the effort on networks. No curve parameter has been learned;
the training rounds (specs/116–120) fitted correction tables on top of the model's output, not parameters inside it.

## 2. Today's equations

Internal scale 0–1; every node rests at 0.01 (UI 1), so the ceiling is 100× (specs/124 measures how much depends on it).

| element | today | fixed or tunable |
|---|---|---|
| AND (inputs of one reaction) | product of folds, smooth cap at the top of the scale (`hill_sat`, specs/010, 074) | fixed |
| OR (alternative producers) | mean of folds | fixed |
| inhibitor | output ÷ inhibitor fold (`divide`), de-repression ≤ 10× | fixed (10× arbitrary) |
| depletion | as inhibitor, suppression ≥ 0.1× | fixed |
| loops | iterated to a fixed point, with structural rules: cycle-balance pools, carrier loops, twins, self-inhibitor and self-activator weights, cofactor / drug / pathogen holds | rules are structure |
| perturbation | root pins: 0 (knockout), 80× (over-expression); drugs as knockout or activity cut (specs/122–123) | fixed doses |

Every edge passes a change at full strength, so a positive loop at rest has gain exactly 1 (the knife-edge behind
rails, specs/014, 054, 082), and nothing attenuates a signal with distance (the false change calls). Several loop rules
patch in conservation (a pool's total, a catalyst returned unchanged) or identity (two copies of one entity) that the
fold-change formalism does not have.

## 3. Proposed model

Three of Adam's concerns are one change: from "fold change in a generic cell" to "levels in a specific cell".

**State.** Each node n has a level x_n ∈ [0, 1], a fraction of the most it can reach. Each cell c has a resting level
r_{n,c}: near 0 for a gene that is off, from data where we have it (§5), from a class default otherwise. The ceiling
on a node is 1 / r_{n,c} — measured, not chosen — and "1 = normal" disappears: a fold change is an output, x_n
(perturbed) / x_n (unperturbed), for that cell.

**Reaction output** (the 2025 design, Section 2):

- activators combine as today (AND: product of normalised inputs; OR: mean of alternatives);
- each input enters through a Hill curve, y = s^h / (s^h + K^h), with K (half-point) and h (steepness);
- an inhibitor suppresses by 1 / (1 + (I / K_I)^m);
- parameters are shared by **edge class** (catalysis of a phosphorylation, binding, transcription, degradation,
  depletion, ...) and may deviate per edge only where data reach that edge, with a penalty that pulls deviations to 0.

**Resting state is a steady state.** The parameters must make the unperturbed cell's resting levels a solution. Given
r, this fixes each node's scale and anchors K near its input's resting level, as the 2025 design's "K anchored near
s_baseline". It is a constraint that ties K to measurements and reduces what training must pin down.

**Loops.** The structural rules (pools, carriers, twins, holds) stay: they decide *what each reaction reads*. The
learnable curves sit on top. Stability is the 2025 Section 9 rule, enforced: in every strongly connected component the
loop gain at rest is held below 1 (a penalty or a parametrisation), so every loop has one answer and no rail. A rail
then means the data asked for one, not that the arithmetic had nowhere else to go.

**Starting point.** Parameters start where the model reproduces today's calls as closely as the new forms allow (the
"untrained" model is the frozen one, approximately). Every learned change is measured against it.

**Time.** Each node: dx_n/dt = (F_n(x) − x_n) / τ_n. At rest the left side is 0 and the model is the steady-state model
above, with the same K, h and r. Time adds one number per node (τ, shared by class first): how fast it responds.

**Rates.** K reflects the balance of forward and reverse rates (a Km or binding constant, relative to the resting
level); h reflects cooperativity or a multistep process; τ reflects absolute speed. Steady-state data identify ratios;
only time courses identify speeds.

**Flux (later).** A reaction consumes its inputs and produces outputs (the 2025 Section 4). Conserving material makes
explicit what the pool and carrier rules approximate today; it is the last step, once the time model works.

## 4. Solver

The current solver has hard switches (rule branches, iteration caps, the Levenberg–Marquardt fallback and restarts),
so it gives no clean gradient. The learnable model needs a smooth forward map, a fixed point found by iteration (it
converges because loop gains are below 1), and gradients through the fixed point by implicit differentiation (one
linear solve per fit step; ForwardDiff is already a dependency). The qualitative propagator stays as the reference
model and as the explanation view for the Pathway Browser.

## 5. Data

Roles follow specs/051 (`bench/datasets.tsv`): FIT may be fitted, SELECT chooses, TEST is scored once.

| parameter | what it means | data we hold (role) | candidate additions (each needs Adam's approval) |
|---|---|---|---|
| resting level r | per node, per cell | DepMap expression (approved), CPPA untreated samples (training portion) | Gonçalves et al. 2022 *Cancer Cell* (949 lines, DIA proteomics, PRIDE PXD030304); Nusinow et al. 2020 *Cell* (CCLE proteomics); CCLE RPPA (same antibody family as CPPA) |
| curve shape K, h | how strongly, and where it levels off | decryptM (doses; FIT) | Kuster lab 2025 decryptM extension (bioRxiv 2025.11.18.689017; 133 inhibitors × 5 lines × 11 doses) — check first it is not already inside our decryptM pull; it cannot be a test either way (same lab and method) |
| interactions (AND, inhibition) | how two inputs combine | korkut2015 drug pairs (FIT) | — |
| breadth across drugs and cells | shared parameters | CPPA training portion (79 drugs × 238 lines) | — |
| speed τ | time constants | HPN-DREAM time courses (FIT) | Tognetti et al. 2021 *Cell Syst* (67 breast lines, 10 points over 60 min, 5 inhibitors, CyTOF); Gross et al. 2022 *Commun Biol* (LINCS MCF10A, RPPA 0–48 h; Synapse terms unverified) |
| drug action | knockout / activity / partial | CPPA training, decryptM | — |

Access and licence for every candidate were read from public pages only (no download, no terms accepted); details
marked unverified there stay unverified until checked.

**Coverage constrains everything.** Only part of each network is measured (HPN-DREAM: 20 of 48 readouts in a network),
so parameters are learned by class first; per-edge parameters only where data reach.

## 6. How it is tested

Training on everything we hold would leave nothing to test. Three layers:

1. **Choices inside training** (forms, sharing, regularisation): cross-validation over held-out drugs and cell lines
   within the fitting data, never the test sets.
2. **Transfer between datasets** (internal, not proof): fit on CPPA training + decryptM, predict korkut and HPN-DREAM
   unfitted, and rotate.
3. **Tests, frozen code and parameters, scored once:**
   - the **CPPA test portions** (52 drugs, 79 cell lines, locked since specs/121): the primary test;
   - **Release 98** curator expectations, when published;
   - a **reserve dataset**, chosen now by coverage alone and not read until the end. Candidate: **Lun et al. 2019
     *Mol Cell*** (doi:10.1016/j.molcel.2019.04.021): 649 kinases and phosphatases over-expressed one at a time in
     HEK293T, CyTOF with ~30 phospho readouts and the over-expressed protein measured in each cell (a dose per cell);
     open, CC BY 4.0 (Mendeley Data doi:10.17632/3kh7ypz232.1). It covers the model's weakest area (genetic
     perturbation; everything we hold is drugs) and no dataset we use overlaps it. Caveat: over-expression in HEK293T is
     artificial. A second reserve could be Gross 2022 or Feng et al. 2023 (MCF10A, EGF and inhibitor dose and time
     series; data route unverified);
   - a **prospective experiment** (§7), the strongest test.

Rules: code and parameters tagged before any test is read (as specs/111); each test scored once by a scoring written
in advance; a reserve is chosen by coverage only (which readouts are in our networks; no values), so choosing it does
not spend it.

## 7. A prospective test (draft; no one contacted)

Predictions written down before an experiment exists are the most convincing evidence a model can give.

- **Design.** One or two cell lines with good baseline data (in CPPA or DepMap); 20–40 perturbations, chosen from the
  model's predictions to be informative: where it disagrees with reachability and with the drug-blind profile, plus
  cases where it predicts no change. Mix knockdowns, over-expressions and drugs at two doses; untreated and vehicle
  controls; replicates to estimate the ceiling.
- **Readout.** An RPPA panel (the platform of CPPA, ~200–400 antibodies), or a phospho-flow or western panel for a
  smaller study.
- **Sealing.** Predictions, model tag and scoring script committed publicly with a hash (or an OSF pre-registration)
  before samples are run; the lab doing the work does not see them until the data are in.
- **Scoring, fixed in advance.** Direction calls (balanced accuracy) and per-perturbation correlation, against
  reachability, the drug-blind profile and the replicate ceiling.
- **Who.** A collaborating lab or a core facility running RPPA as a service. Choosing and contacting one is Adam's.

## 8. Stages (each its own registered spec)

| stage | what | judged by |
|---|---|---|
| 0 | specs/124: how much rests on the 100× ceiling | registered predictions P1–P5 |
| 1 | smooth Hill forms, class-shared parameters, contraction in loops; start reproducing today's model | agreement with the frozen calls on curator and experimental (no-harm), behaviour checks (specs/107: convergence, monotonicity, no rails) |
| 2 | resting levels per cell from data; fold change as an output | CPPA training cross-validation; genes that are off now produce no change |
| 3 | learn class parameters (and drug action) | cross-validation by drug and cell line; transfer between datasets |
| 4 | time constants τ | HPN-DREAM (and Tognetti, if approved) cross-validation on held-out time points and lines |
| 5 | one scoring of the frozen learned model | CPPA test portions, Release 98, the reserve; the prospective test when it exists |
| 6 | material conservation (flux) | as stage 1, then the stage 5 tests |

## 9. Risks

- **Identifiability.** Many parameters, sparse readouts: class sharing, priors from today's model, and the resting-state
  constraint are the defences; report parameter uncertainty, not only fits.
- **Overfitting to SELECT and FIT data.** Every choice so far has seen the curator and experimental sets; only §6 tests
  answer whether learning helps.
- **The test budget.** The CPPA test portions are spent by one scoring. If the learned model fails there, the next
  attempt needs the reserve or new data.
- **Compute.** Differentiating through 92 networks with large loops is heavier than today's solve; stage 1 measures it.
- **The rules interact.** Learnable curves under structural rules each measured against old defaults may behave
  unexpectedly; stage 1's no-harm and behaviour checks come before any fitting.
