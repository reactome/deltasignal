# specs/128: response curves on edges (specs/125 stage 1)

**Status:** registered 2026-11-02, before any code. First stage of the learnable model (specs/125 §8). Nothing is
fitted here and nothing is adopted: this stage builds the curves, shows that today's model is one setting of them, and
measures how far a few fixed settings move the calls, so the fitting stage (3) starts from known ground.

## Why

Today every edge passes a change at full strength: a 2× input gives a 2× output, at any distance and any level
(specs/125 §2). That is what makes a positive loop's gain at rest exactly 1 (rails, specs/014, 054, 082), and what lets a
change reach far from its target undiminished (the false change calls, specs/115). The 2025 design gave each edge a
response curve with a learnable shape. This stage puts the curves in, fixed by class, with today's behaviour as one
point in their range.

## The curves

Each reaction slot reads its input's fold f (its level over its resting level), after every structural rule has
decided what it reads (specs/125 §3: the rules decide *what* a reaction reads, the curve *how strongly*). The curve
replaces f by g(f):

- **Activating slot** (substrate, catalyst, positive regulator):
  g(f) = f^h · (1 + K^h) / (f^h + K^h).
  g(1) = 1 (rest stays rest), g(0) = 0 (a knockout stays a knockout). Slope at rest, in log terms (elasticity):
  e = h·K^h / (1 + K^h). The most the edge can pass up is 1 + K^h. K is where the input rests relative to its half-point
  (K = 9: the input rests at a tenth of saturation). As K → ∞, g(f) = f^h; with h = 1 it is today's identity.
- **Inhibiting slot** (negative regulator, under `DS_INHIBITION_MODE=divide`):
  factor = (1 + k^m) / (1 + k^m · f^m), which is the 2025 design's 1 / (1 + (I/K_I)^m) normalised to 1 at rest. The most it
  can de-repress is 1 + k^m. Today's `divide` with its 10× cap is close to k = 9, m = 1 (exactly 10 at f = 0, 0.526
  instead of 0.5 at f = 2).
- **Depletion slot:** the same form as an inhibitor; the depletion floor (`DS_DEPLETION_H_MIN`, 0.1×) is kept as a
  clamp, because specs/011 measured it.
- **Never curved:** output edges and the identity links (variant pool and split, specialization, set member, recycled
  and RNA supply, carrier twins, recycle). These say two nodes are the same material, not that one responds to the
  other.

A slot that is both a reaction's input and its catalyst (a complex modifying itself, specs/109–110) takes the
substrate curve.

**Classes** (parameters shared within a class, per specs/125 §3): `input`, `catalyst`, `regulator_pos`,
`regulator_neg`, `depletion`. Assembly and composition slots, where present, take the `input` curve.

**Configuration:** `DS_EDGE_CURVES=off` (default; byte-identical to today) or the path of a TSV with one row per class
(`class, K, h` for activating classes; `class, k, m` for inhibiting ones). A class not listed, or listed with `K = inf`
and `h = 1` (`k = divide` for inhibiting classes), keeps today's behaviour exactly. A malformed table, an unknown
class, or a non-positive parameter is a startup error. The solve reports `edge_curve_rule` (off, or the table's
sha256 and contents) and the number of slots curved per class.

## Arms (canonical frozen build `20261026-1349_8b6527e`; curator, experimental and korkut axes)

| arm | table | what it asks |
|---|---|---|
| `ec_off` | off | control |
| `ec_id` | every class listed at its identity setting | the machinery is exact (D0) |
| `ec_inh` | `regulator_neg` and `depletion` k = 9, m = 1; activating classes identity | today's inhibition, smoothed |
| `ec_sub90` | `input` K = 9, h = 1 (e = 0.9, at most 10×); others identity | substrates below saturation |
| `ec_all90` | `input`, `catalyst`, `regulator_pos` K = 9, h = 1; inhibiting k = 9, m = 1 | mild contraction everywhere |
| `ec_all75` | as `ec_all90` with K = 3 (e = 0.75, at most 4×) | strong contraction |

Reported for each arm against `ec_off`: curator held-out and tuning (RAF set aside and kept), experimental, korkut
(FIT; per the specs/062 protocol), McNemar p, pathways moved, converged solves, the share of readouts at a rail (fold
≥ 0.99 × ceiling or ≤ 0.01), and mild-dose rail-dependence (specs/108) for `ec_all75`.

## Predictions (written before any arm runs)

- **P1.** `ec_id` reproduces `ec_off` exactly on every case (validity, call, value as printed, convergence).
- **P2.** `ec_inh` changes fewer than 1% of held-out curator and experimental calls.
- **P3.** The contraction arms lose held-out curator calls, as blanket loop elasticity did (specs/094: −40): `ec_all75`
  is net negative at p < 0.05. Correct moderate calls fall back below the 1.15 / 0.85 cutoffs once each hop passes only
  part of a change.
- **P4.** Contraction removes rails: the share of readouts at a rail falls by more than half under `ec_all75`.
- **P5.** On korkut (drug knockdowns, where specs/114–115 found too many change calls far from the target), `ec_all75`
  makes fewer change calls than `ec_off`, and its accuracy is not lower.

## Reading (fixed now)

Nothing here is adopted: the default stays `off`. P1 failing blocks stage 3 until it is fixed. P3–P5 tell stage 3 where
its starting values and priors should sit: if contraction costs the curated calls, fitted curves must earn their
contraction from data, and the identity setting is the prior. Fitting (stage 3) uses only FIT data (specs/051); the
curator held-out split stays SELECT.
