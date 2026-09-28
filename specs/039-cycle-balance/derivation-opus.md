# Cycle rule, derived independently (Opus), before reading the Fable derivation

## Derivation

Take a two-state cycle A ⇄ B with total T = A + B, forward flux k_f·f·A and
backward flux k_b·b·B.
- f is the product of the forward reaction's non-cycle positive inputs: the E_f
  fold and the co-input folds, divided by the forward inhibitor folds (the same
  AND/divide semantics as a reaction node, excluding A itself).
- b is the same quantity for the reverse reaction, excluding B.

Steady state gives k_f f A = k_b b B, so B/A = r0·ρ with ρ = f/b and r0 = B0/A0
the baseline ratio.

With the modified fraction at baseline φ0 = r0 / (1 + r0), and folds relative to
baseline (t = T/T0):

    B_fold = t · ρ / (φ0·ρ + (1 − φ0))        A_fold = t / (φ0·ρ + (1 − φ0))

Check: at ρ = 1, both folds equal t.

**The supply t** is the fold of the cycled protein's external supply: the mean
of the producers of A and B that lie OUTSIDE the cycle, or 1 if there are none
(the cycled protein at baseline). For a root form of the gene, t is its pin.

**φ0 is the one assumption.** Take **φ0 = 1/2**, the uninformative equal split.
- ρ = 2 (kinase doubled): B = 1.33, A = 0.67.
- ρ → ∞ (phosphatase KO): B → 2, A → 0.
- ρ → 0 (kinase KO): B → 0, A → 2.

With a **small φ0** (the modified form a minority, often true of activated
forms): B ≈ t·ρ (capped at 1/φ0) and A ≈ t. That is closer to today's
multiplicative behaviour for B, but it removes the feedback through A.

## Sensitivity of the calls

- **Direction** is independent of φ0 for any ρ ≠ 1: B moves with ρ, and A moves
  against it. So expectations 1 and 2 hold for every φ0 in (0, 1).
- **Magnitude** depends on φ0. With φ0 = 1/2, B is capped at a 2-fold rise.
  - The fold cutoff is 1.15, so B = 1.33 at ρ = 2 is UP.
  - An upstream 80-fold activator gives B ≈ 2: UP, but not large, which weakens
    propagation downstream through AND products.
  - With φ0 = 0.1, B can reach 10-fold.
- **Expectation 3** (gene overexpressed → both UP) holds through t.
- **Expectation 4 is the risk.** B's rise is bounded by 1/φ0, so a strong
  upstream signal is compressed at every cycle it passes through. This is the
  real trade-off against today's unbounded multiplication.

## Composition with the network

- The cycle's two entity nodes (A, B) are computed by the rule; they are no
  longer iterated.
- Reaction nodes F and R read the resulting A and B for downstream display.
- Downstream consumers read B (or A) as ordinary entity values.
- ρ is read from the current values of E_f, E_b, the co-inputs and the
  regulators. In a larger cyclic component those can themselves depend on B, so
  the rule replaces only the A → F → B → R → A closure, not the whole component.
  The remaining feedback (B → … → E_f) is iterated as today, but it no longer
  passes through a gain-1 AND chain.

## Cases

- **One direction catalysed.** The uncatalysed direction has b = 1 (a constant
  basal rate). Then ρ = f, so the enzyme on the curated side controls the split.
- **Multi-state A → B → C → A.** The steady state of a linear chain with a
  return: fractions proportional to the products of the upstream rate ratios
  (the standard Markov-chain stationary distribution on the cycle), with a
  uniform baseline split of 1/n.
- **Binding / dissociation (A + L ⇄ AL).** This is not a pure interconversion:
  the ligand L is a second conserved species. Treat it as an interconversion of
  A with ρ including L's fold (L as a co-input of F), and cap at the minority
  species. This is an approximation, to be flagged in the results.

## Failure modes

- **Compression.** With φ0 = 1/2, strong upstream signals are compressed (see
  expectation 4 above).
- **Wrong cycle detection.** Two reactions whose A and B share a stId but are
  not the same molecule pool, such as different compartments merged by id.
  This needs the uuid identity, not the stId.
- **Co-input double counting.** A co-input that is ALSO consumed elsewhere is
  still read by fold.
- **Supply ambiguity.** When A and B have external producers that are
  themselves modified forms, t should be their mean. It is an approximation.

## Fixtures (φ0 = 1/2, t = 1 unless stated)

| case | ρ | B | A |
|---|---|---|---|
| nothing perturbed | 1 | 1 | 1 |
| E_f at 2x | 2 | 1.333 | 0.667 |
| E_f knocked out | 0 | 0 | 2 |
| E_b knocked out | → ∞ | 2 | 0 |
| E_b at 2x | 0.5 | 0.667 | 1.333 |
| gene at 2x, enzymes at baseline | 1 | 2 (t = 2) | 2 |
| E_f at 2x and gene at 2x | 2 | 2.667 | 1.333 |
| one direction catalysed, E_f at 3x, b = 1 | 3 | 1.5 | 0.5 |
