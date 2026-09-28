# specs/039 — Interconversion cycles as conserved pools (design; not yet pre-registered)

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
