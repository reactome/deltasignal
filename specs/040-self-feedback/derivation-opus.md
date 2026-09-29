# specs/040: derived independently (Opus), before reading the Fable derivation

## 1. The three amplifiers are two defects, not one

Working the RAF trace through, amplifiers 1 and 2 are feedback from a
reaction's own downstream, and amplifier 3 is not.

- **Amplifiers 1 and 2: recycled inputs.** The p-MEK dimers in the catalyst pool,
  and the 22 scaffold leaves of the bundle, are produced **only** by reactions
  downstream of the step they feed ("Dissociation of RAS:RAF complex"). They
  are catalytic recycling: the species leaves the step's product and returns to
  its input.
- **Amplifier 3 is a shared-input double count.** The PEBP1 inhibitor is built
  from RAS:GTP and RAF1, which are *upstream* inputs of the step it inhibits.
  That is specs/022's class, "an inhibitor containing its own reaction's input",
  missed only because 022 matches whole nodes, not shared leaves. It needs a
  detection clarification of an adopted rule, not a new rule.

So I do not think one rule honestly covers all three. The candidate wording ("a
node's own downstream must not multiply back into it") fits 1 and 2 and would
be stretched to fit 3.

## 2. Rule A: recycled inputs are read open-loop

**Definition.**
- Take a reaction r with target t, inside a strongly connected component.
- An input i of r (an AND activator, or a member of a set-pool node feeding r)
  is **recycled** if **every** producer of i is reachable from t by activator
  (mass-flow) edges inside the component, not through a pinned node.
- An input with any producer outside t's downstream is not recycled and is
  untouched.

**What it reads.** In r's evaluation, i reads its **open-loop value**: i
re-evaluated from its own producers with t held at baseline. This is one pass,
in a fixed order (BFS from t over sorted ids), exactly the machinery of
specs/039 amendment 5.
- i keeps every influence from outside the loop, such as a pinned scaffold gene,
  or MEK itself perturbed upstream of the dimer's production.
- It loses only the part that is t's own output coming back.

**Genuine feedback (question 2).** MEK → RAF is real positive feedback. In a
fold model a positive loop with gain at or above 1 has no bounded answer: it
rails to 100, or to 0 once another multiplied leaf dips, which is what happened.
The open-loop reading is the first-order response, i.e. the feedback's
contribution to direction, without its runaway amplification.
- The sign is kept for every perturbation entering the loop from outside.
- The magnitude is under-stated wherever real feedback amplifies.
- It is bounded by construction, because a node no longer multiplies itself.

**Relation to earlier specs:**
- specs/018 read derived-edge closures at their *entry* value, which is roughly
  baseline inside a component: close to this rule, but with no re-evaluation.
- specs/039 amendments 3–5 are the pool-local version: a pool's own output is
  not supply, and pool-fed inputs are re-evaluated.
- Rule A is that principle applied to any reaction, restricted to inputs with no
  outside producer.

**The main risk, and a prediction.** Amplifier 2's scaffold leaves are released
**unchanged** by the dissociation. That is exactly specs/035's
`DS_CONSERVED_MODE=inert` class, which held such inputs at baseline and lost
**experimental −122**.
- Rule A differs only in re-evaluating instead of holding. For a scaffold with
  no external perturbation, the two give nearly the same value.
- **So Rule A will reproduce much of 035's loss wherever it touches 035's set**,
  unless the census shows the overlap is small.
- p-MEK dimers (amplifier 1) are transformed inside the loop, so they are not in
  035's set.
- **Pre-arm check:** census A's recycled inputs against 035's `conserved_held`
  set on the same build. If most of A's catalog-wide footprint is 035's set,
  expect 035's loss, and narrow to **A1: set-pool members only**. That covers
  amplifier 1, which is new, and leaves the bundle leaves.

## 3. Rule B: a leaf-sharing self-contained inhibitor (a specs/022 clarification)

- An inhibitor of r is self-contained if its leaf set shares a leaf with the
  leaf set of one of r's inputs, and that shared input reaches the inhibitor by
  mass flow. This replaces 022's whole-node containment test.
- The adopted rule then applies unchanged: the overlapping log-fold is kept at
  weight w = 0.1.
- For RAF: the PEBP1 complex's RAS:GTP and RAF1 leaves are shared with the
  step's RAS:GTP:RAF dimer input, so the inhibitor reads roughly the input's
  fold^0.1, not its full fold.

## 4. RAF, worked (qualitatively, from the traced values)

- **A alone:** the catalyst pool's p-MEK members read about 1, so the pool
  stops being fold³. The bundle leaves read about 1, so the step stops being
  0.47²². RAS:GTP:RAF dimer 80x → activated dimer UP.
  - The next step is still divided by the PEBP1 inhibitor tracking the input.
    The trace's "hold both" row gave 2.47x, which is UP, but weak.
- **A + B:** PEBP1 damped to fold^0.1. The trace's "hold PEBP1" row gives the
  step 100, but the MEK/MAPK tier is still blocked by drug-bound complexes
  (0.001).
- **Without amplifier 4 (out of scope) the p-MAPK readouts may stay DOWN.**
  Honest expectation: the activated-dimer readout 5672718 and RAS-level cases
  flip; the p-MAPK dimer readouts may not. That is perhaps 7–10 of the ~19
  experimental cases, not all of them.

## 5. Census before any arm (Q5)

For each rule, on the canonical build:
- the number of recycled inputs (A) and leaf-sharing inhibitors (B), by pathway;
- the overlap with specs/035's conserved-held set (A), and with specs/022's
  current flags (B, which should be a superset);
- which scored readouts lie downstream of an affected reaction, by pathway, with
  their counts.

## 6. Arms and gates

- **Arms:** ctrl; A1 (set-pool members only); A (all recycled inputs); B; A1+B.
- **Gates**, as for specs/039:
  - curator held-out > +15, p < 0.05;
  - experimental ≥ 0, not significantly negative;
  - no pathway loses > 10;
  - gains span ≥ 2 pathways.
- RAF is reported separately, and the gates are also computed excluding RAF.

## 7. What can go wrong

- **A removes load-bearing feedback:** positive loops that are genuine
  amplifiers (MAPK cascades, the insulin/IRS loop). Expect magnitude loss, not
  sign loss.
- **Label independence:** the BFS order and "every producer downstream" are
  graph properties. Only the order of re-evaluation is fixed by sorted ids.
- **Cost:** re-evaluation per affected reaction per sweep is bounded by the
  census counts.
