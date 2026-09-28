# How the DeltaSignal model works

A plain explanation of the model as it stands, for explaining it to others.

- **Why each choice was made:** in the `specs/` directory, cited inline. Every
  choice below was measured before it was adopted.
- **The numbers:** `docs/RESULTS.md`.
- Each section says whether it is **implemented** or **designed, not yet
  measured**.

## 1. The question the model answers

Given a Reactome pathway and a perturbation (a gene knocked out, or
overexpressed), which downstream entities go UP, stay NORMAL, or go DOWN?

The model propagates the perturbation through a logic network generated from
Reactome, and reports a fold-change for every node.

## 2. The network (implemented)

The logic-network generator turns a Reactome pathway into a directed, signed
graph.

- **Nodes** are physical entities (proteins, complexes, small molecules) and
  reactions.
- **Positive edges** are inputs, outputs, catalysts and positive regulators.
- **Negative edges** are negative regulators and derived depletion edges.

**Entities get a node per place, not per Reactome id.**
- In a Reactome diagram, the same entity can appear in several places that are
  not connected to each other.
- The generator gives each separate occurrence its own node (a uuid), so
  occurrences the diagram does not connect are not merged.
- Merging them by id was a major source of artefact loops in the earlier
  MP-BioPath networks, which were then cut by hand. The uuids already resolve
  about 86% of those hand-cut loops (specs/034 §11).

**Sets.** A Reactome set means "any one of these".
- **Set inputs and outputs** are expanded into virtual reactions, one per
  member. That lets one member be perturbed and followed.
- **Set-valued catalysts and regulators** are one pool node, the set itself.
  Its members feed it, and it feeds each reaction once (specs/033, 038). A
  reaction is made unique by its inputs and outputs, not by its regulators.

## 3. Propagation (implemented)

Every node carries a **fold-change relative to normal**: 1 is normal, 0 is
absent, and the cap is 100.

- **A perturbation is a pin.** The gene's root input is fixed at 80
  (overexpressed) or 0 (knocked out). This is the protocol of the MP-BioPath
  publication.
- **A reaction** multiplies the fold-changes of everything it requires
  (substrates, catalysts, positive regulators; AND), and divides by each
  inhibitor's fold-change. Doubling a required input doubles the reaction;
  doubling an inhibitor halves it; removing a required input stops it.
- **An entity** made by several reactions takes their average (OR).
- **A set pool** multiplies its members' fold-changes (`DS_SET_POOL_MODE=product`).
  So every member's change passes through: one member doubled doubles the
  pool, and one knocked out stops it. Treating members as interchangeable
  alternatives (the maximum, or the average) was measured and is worse on both
  benchmarks: the ground truth expects a single member to matter (specs/038).
- **Cofactors** (ATP, ADP, water, …) take part in reactions but are held at
  normal, so a perturbation cannot travel through them (specs/007).

The network is split into strongly connected components (loops), and they are
solved upstream first. Inside a loop, the solver iterates until the values
stop changing.

## 4. Loops, and why they need their own maths (measured, NOT adopted: specs/039)

> **Status.** Everything in this section is implemented behind
> `DS_CYCLE_MODE=balance` and was measured on both benchmarks. The default
> remains `off`: loops are iterated as in §3.
> - **As first pre-registered, it was not adopted:** curator held-out −5,
>   experimental −8 against the rule off.
> - **With the pool boundary fixed (amendment 5):** curator held-out is +77
>   (p 1.8e-16), and experimental is −2 (not significant). That fails the
>   pre-registered experimental gate, so it is **held off by default** (Adam,
>   2026-09-28).
> - The amendment 5 fix: a regulator built from a pool's own form keeps its
>   own effect but loses its tracking of the pool, and a multi-step path must
>   regenerate its enzyme.
> - **The remaining loss:** readouts on a pool's resting form. At φ₀ = 0.1,
>   driving a pool toward its modified form moves the resting form only to
>   1.11.
>
> The numbers and the traced losses are in
> `specs/039-cycle-balance/research.md`.

**What the loops are.** Once the artefact loops are separated by uuids, most of
what remains are **curated interconversion cycles**, where one reaction turns a
protein from form A into form B and another turns it back:

| cycle | A | B | forward enzyme | backward enzyme |
|---|---|---|---|---|
| RAS activation | RAS-GDP | RAS-GTP | GEF | GAP |
| phosphorylation | protein | phospho-protein | kinase | phosphatase |
| ubiquitination | protein | Ub-protein | ligase | deubiquitinase |
| methylation | histone | methylated histone | methyltransferase | demethylase |

Reactome usually curates such a cycle in steps: the enzyme binds, modifies and
releases (S + kinase → S:kinase → p-S:kinase → p-S + kinase). The count is
recorded per build in specs/039.

**Why the current maths gets them wrong.** Propagation multiplies around a
loop: A feeds the reaction that makes B, and B feeds the reaction that makes A.
At normal levels every step multiplies by 1, so a loop is balanced on a knife
edge. A small push runs it up to the cap or down to zero, and which way can
depend on the order the solver visits nodes in. In the RAF/MAP kinase pathway
this collapses the whole cascade to zero under every upstream perturbation.

The root problem is that the model treats "how fast B is made" as "how much B
there is". A cycle has no sense of scale without the one fact the maths is
missing: **the total amount of the protein is conserved; the cycle only moves
it between forms.**

**The fix: each protein's forms are one conserved pool, solved at steady
state.** Treat the forms of one protein as the states of a Markov chain:

- **States:** the protein's curated forms (A, B, and more for multi-site
  proteins).
- **Transitions:** the curated reactions between forms. Each transition's rate
  is its baseline rate × its **drive**: the enzyme's fold-change × the other
  inputs' fold-changes ÷ inhibitors. These come from the network, as for any
  reaction.
- **Stationary distribution:** π, with π = πP (or πQ = 0 in continuous time),
  gives the fraction of the pool in each form, once the forward and backward
  flows balance.
- **Fold-change of each form** = (perturbed π ÷ normal π) × (the fold-change in
  the protein's total supply, e.g. its gene overexpressed).

For a two-form cycle this reduces to a closed form, derived independently by
two models (specs/039):

    r = forward drive / backward drive
    modified form B  = supply × r / (φ₀·r + 1 − φ₀)
    unmodified form A = supply     / (φ₀·r + 1 − φ₀)

φ₀ is the fraction of the protein in the modified form at rest. **φ₀ = 0.1**
(decided by Adam, 2026-09-27): activated or modified forms are a minority in a
resting cell (RAS is about 5–10% GTP-bound; basal phosphorylation is usually
≤10–20%). It is fixed in advance, not fitted.

**Cycles curated in steps.** The states of the chain are the protein's
distinct modification states (S and p-S; RAS:GDP and RAS:GTP). The enzyme
complexes on the way (S:kinase, RAS:GTP:GAP) are intermediates:
- A transition is the whole path from one state to the other. Its drive is the
  product of its steps' drives, so the kinase acts through the binding step, in
  proportion to its fold.
- Intermediates hold no share of the protein at rest. They are transient, and
  their level is set by the flux through them, which is what they read.
- With that, a multi-step cycle gives exactly the numbers of the two-form rule
  below. A share of 1/9 for each complex would instead cap a kinase-driven rise
  at about 3.7-fold.
- An uncatalysed step (intrinsic GTP hydrolysis) beside an enzyme-driven one
  carries a 1/1000 share at rest, because it is orders of magnitude slower.
- The enzyme's free form is released by the cycle it drives. Counting that
  release as a new supply of enzyme would close a loop of its own, so the free
  form reads only its other producers (`DS_CYCLE_CARRIERS`, tested separately).

**What this gets right:**
- **Direction:** kinase or GEF up → modified form UP; phosphatase or GAP knocked
  out → modified form UP; kinase knocked out → modified form to zero. This holds
  for any φ₀.
- **Gene overexpression:** both forms rise together, because the pool is
  bigger.
- **No knife edge and no collapse to zero:** each pool is bounded (the modified
  form can rise at most 1/φ₀ = 10-fold), and normal is exactly normal.
- **Cascades still propagate:** with φ₀ = 0.1, an 80-fold signal stays UP
  through about 20 stacked cycles. With φ₀ = 0.5 it would be compressed to
  NORMAL within four (RAS → RAF → MEK → ERK), which is why 0.5 is only a
  sensitivity check.
- **It is what the MP-BioPath curators did by hand.** As φ₀ → 0 the rule
  reduces to cutting the loop and adding "backward enzyme ⊣ modified form",
  which is the edit they made in 219 places. The rule derives this from the
  curated reactions wherever the reverse reaction is curated (about a quarter
  of their additions; specs/034 §12b).

**What it does not cover:**
- Conservation holds within a protein's pool, not across the network: a signal
  can still amplify from pool to pool, through enzymes.
- Between pools, propagation stays as in §3. The forms of one pool are the
  enzymes driving the next pool's transitions (RAS-GTP drives RAF's, RAF drives
  MEK's).
- The outer solve still iterates, but every pool is now bounded.
- **Branched pools** (a protein modifiable at two sites) need one extra stated
  assumption about baseline rates: an equal split of flow between exits, or
  detailed balance.

**Known limitation of the implementation:** influence scores (explainability)
do not see the pool rule; under `balance` they are computed from the ordinary
reaction model.

**Who does what:**
- **The generator** detects each pool once, when it builds the network from
  Reactome: which nodes are forms of one protein (a shared reference entity),
  and which reactions convert between them. It ships the list with the network.
- **DeltaSignal** solves π = πP for each pool on every run, because the drives
  depend on the perturbation. The generator is the side that can read "these
  are the same protein" reliably; the solver would have to guess it from ids.

## 5. What is compared, and how

- **Ground truth:** curator predictions (about 24,000 cases) and published
  experimental results (849 cases), from the MP-BioPath publication.
- **Split:** the paper's ten pathways are used for tuning; the other 70 are
  held out and reported.
- **Discipline:** every change is pre-registered with its adoption criteria
  before it is measured, and measured against a freshly regenerated control on
  both axes.
- **Current numbers:** `docs/RESULTS.md`.
