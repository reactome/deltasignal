# specs/040: a node's own downstream multiplied back into it

**Status:** pre-registered 2026-09-28, before any arm. Problem: `problem.md`.
Method: a blind double derivation. `derivation-opus.md` was committed
(25ca753) before `derivation-fable.md` was read.

## The derivations agree

The RAF collapse is **two defects**, needing two rules:

- **Rule A:** self-fed (recycled) inputs are multiplied back into the step
  whose downstream produces them. This covers amplifiers 1 and 2.
- **Rule B:** a self-contained inhibitor that shares *leaves*, not a node, with
  its step's input. This is specs/022's class, missed by its whole-node test.
  It covers amplifier 3, and also damps amplifier 4 (drug-bound copies of the
  step's own input) without the refuted "inert drugs" rule.

## Where they differ, and the choice (Fable's on each)

| point | Opus | Fable (adopted) | why |
|---|---|---|---|
| detection of A | static: every producer of i reachable from the step's target | **per solve:** i is self-fed iff every activator path from a signal-carrying entry to i passes through Out(i), the targets i feeds | Fable ran both static definitions on RAF and both fail (ATP makes a reaction copy an entry; or no OR-producer enters). |
| what A reads | open-loop re-evaluation | baseline on the edges into Out(i), via specs/018's supply plumbing | These are equivalent under the per-solve definition: a self-fed node has no signal source except through its own outputs. Baseline is simpler. |
| B's reach | shared leaf with the input | a shared non-cofactor leaf whose carrier reaches both the input and the inhibitor | The PEBP1 complex is built from a sibling variant two generations from the input. |

## Rules (as pre-registered)

**A: `DS_SELF_FED_MODE=off` (default) | `entry`.**
- In an iterated strongly connected component, the **signal-carrying
  entries** are:
  - the nodes pinned off baseline;
  - unpinned targets with an activator from outside the component that is off
    baseline.
- A pin at baseline (an inert cofactor, for example) is a constant and blocks
  the search.
- Nodes owned by a specs/039 pool are excluded.
- A node u is **self-fed** iff every activator path from an entry to u passes
  through Out(u). Set pools are transparent to the search.
- u's edges into Out(u) read baseline.
- Recomputed per solve; deterministic and label-free.

**B: `DS_SELF_INHIBITOR_LEAVES=0` (default) | `1`.**
- specs/022's pair test is widened: an inhibitor is self-contained if it shares
  a non-cofactor leaf with the step's input, and the leaf's carrier reaches both
  the input and the inhibitor.
- The same formula applies (w = 0.1, weaken-only).

Off is byte-identical for both. Baseline is exact to 2e-12 on all 92 bundles
(Fable prototype).

## Prototype on RAF (canonical build, balance on; Fable harness)

With **A + B** (activated dimer / MAP2K-binding / dissociation / p-MAPK):
- KRAS 80x: 100 / 69 / 45 / 100, UP.
- BRAF 80x: 100 / 63 / 46 / 100.
- NF1 KO: 1.47 / 1.43 / 3.5 / 100, UP.
- KRAS KO: 0.35 / 0.38 / 0.03 / 1e-6, DOWN.
- HRAS 80x → 100, HRAS KO → 0.004; NF1 80x DOWN.

A alone fixes only the dimer: the drug copies divide the MEK step about 400x.
B alone rails every KO to 100.

**Open defect:** 3 of 8 test cases end at a residual of about 1, meaning
something still oscillates. **It must be traced before the arms.**

## Census (Fable, canonical build)

- **A, over all single root pins:**
  - 67 of 92 pathways, 175 components, 1,803 nodes (165 of them set-pool
    members).
  - It touches 432 of 11,328 pins.
  - Node overlap with specs/035's held set: 340 of 535, 256 of them in Class I
    MHC. The case overlap is small.
  - PIP3 is flagged under 1 of 322 root pins (035 held it always).
  - RAF: 98 nodes, 14 of 325 pins.
- **B:** flagged pairs go from 538 to 1,406, and 41 pathways gain some.
  - RAF goes 0 → 41.
  - **Chromatin goes 0 → 308, via shared histone leaves (promiscuous).**

## Arms and gates

- **Arms:** ctrl, A, B and A+B on one build, at one solver commit. Convergence
  is logged per case.
- **Gates**, for each of A+B, A and B against ctrl:
  - curator held-out > +15, p < 0.05;
  - experimental ≥ 0, not significantly negative;
  - no pathway loses > 10;
  - gains span ≥ 2 pathways and ≥ 5 perturbations;
  - net ≥ 0 outside RAF.
- **Required checks:**
  - relabel churn no worse than ctrl;
  - pathways with zero flags bit-identical to ctrl;
  - Class I MHC and Chromatin named, with their nets;
  - B's per-pathway concentration.
