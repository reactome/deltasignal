# specs/083: an activator built from its own reaction's input counts that input once

**Status:** pre-registered 2026-10-20, before the implementation and the arm. Analysis
from specs/082 and this directory's counts.

Adam, 2026-10-20: "done. continue", after specs/082 (PTEN overexpression in PIP3).

## The defect (traced)

PTEN overexpression in PIP3 (R-HSA-1257604) reads p-AKT 83.5× UP, the same as PTEN
knockout. The experiment and the curators both say DOWN (9 experimental cases). specs/082
found the AKT → MKRN1 ⊣ PTEN ⊣ PIP3 loop, which rails. What sends it the wrong way is
upstream:

1. **PTEN transcription** (R-HSA-8944104) has the PTEN gene as its `input` and five positive
   regulators that are TF:PTEN-gene complexes (TP53, EGR1, PPARG, ATF2, MAF1). Each of them
   contains the gene. Under the AND product the gene's fold is counted six times: PTEN gene
   at 1.01× gives mRNA 1.0605× (1.01⁶), and 80× gives 100× (capped).
2. That mRNA copy builds the 12 miRNA:PTEN-mRNA silencing complexes, which inhibit PTEN
   translation. Their 6% rise against a 1% shared input is treated as independent of the
   mRNA, so specs/022's self-inhibitor discount does not reach it. Twelve full-strength
   inhibitors halve translation (0.4993× at 1.01×), so more PTEN gene means less PTEN
   protein, and the MKRN1 loop rails that direction to 100×.

Specs/022 makes an *inhibitor* built from its own reaction's input count that input at weight
w. Nothing does the same for a positive *regulator or catalyst* built from its own input.
Biologically the rate of transcription follows how much TF:gene complex there is, so the
gene should be counted once.

## Structural count, before the implementation (`count_self_activators.py`, canonical `20261017-1246_5e663d9`)

Reactions with at least one positive regulator or catalyst whose entity contains one of the
reaction's `input` entities:
- **251 distinct reactions in 36 pathways**: 2,188 reaction nodes, 2,739 activator slots;
- TP53 54, RUNX1 31, R-HSA-446728 14, Mitotic G1 14, R-HSA-452723 13, RUNX2 13;
- largest per reaction: R-HSA-9909623 (CD28 family) 13, PTEN transcription 5, R-HSA-9760965 5.

Most have one such activator, which squares the input's fold. specs/012 measured
de-duplicating the *same entity* as catalyst and substrate at −15 held-out. That is
deletion of an identical node; this rule is containment-based damping, as specs/022 is.

## The rule (`DS_SELF_ACTIVATOR_WEIGHT`, default `off`; values `off` or w in [0, 1])

The mirror of specs/022:
- **Pairs.** For each positive `regulator` or `catalyst` edge A → reaction t: the
  `input`-type activators a of t whose entity A contains (containment table), and from
  which A is computed (A reachable from a). Strict containment only; no leaf clause.
- **Damping.** Let n be the number of such activators of t. Each is read in log-fold as
  l_eff = (l_A − e) + (w/n)·e, where e is its overlap with the shared inputs' log-fold
  (same sign as both, magnitude min(|l_A|, |l_s|); 0 if they move apart). Weaken-only, as
  specs/022. With w = 0 the shared part is removed; the activator keeps its independent
  change.
- **Reporting.** The solve reports `self_activator_rule` and `self_activators`.

## Arm

| arm | build | solver | against |
|---|---|---|---|
| `sa` | the canonical rebuilt for specs/081 (recycle pools on) | `--server DS_SELF_ACTIVATOR_WEIGHT=0.1` | that build's current-defaults run (`dflt`) |

Read with `holdout_report.py --set-aside RAF_MAP_kinase_cascade`, on both axes, and checked case by case.

**Gate (as specs/077):** curator held-out net ≥ 0, experimental net ≥ −5, and no pathway
below −10 (RAF set aside).

## Predictions

1. PTEN overexpression in PIP3: at least 5 of the 9 experimental cases read DOWN.
2. PTEN knockout in PIP3 stays UP (9 experimental, currently right).
3. No prediction for the held-out net. It is a 36-pathway change, and specs/012's related
   de-duplication lost.

## Implementation and a disclosed check (before the arm)

**Implementation:**
- `self_activator_setup` (`src/solvers/steady_state.jl`);
- the activator damping in `compute_reaction_output_vec`;
- `IndexedReaction.activator_shared`;
- `ReactionEvalConfig.self_activator_weight`;
- reporting (solve, API, bench; the bench refuses an "on" run that flags nothing).

Tests: `test_self_inhibition.jl`, 98 (+16 for this rule):
- G counted three times: 8× at 2×;
- once at w = 0;
- 2^1.1 at w = 0.1, and 0.5^1.1 below baseline;
- the TF's own change keeps full strength;
- not computed from the input: untouched;
- no containment table: inert;
- invalid values are errors.

Other suites are unchanged (config 217, invariants 121, determinism 126, cycles 516,
loop elasticity 175 + 1).

**Check on PIP3 (`pten_keys.jl`), not the arm.** With w = 0.1, PTEN transcription follows
the gene (1.01× gives 1.01×, 2× gives 2×; before, 1.06× and 59.7×). PTEN overexpression
still reads p-AKT 50.5× / 83.5× because of two other defects:
- **PTEN mRNA exists as two copies.** The starting copy, pinned by the benchmark, feeds
  translation. The transcribed copy builds the 12 miRNA silencing complexes. The silencers
  are not computed from the copy translation reads, so specs/022 does not treat them as
  self-inhibitors, and they inhibit at full strength (translation 0.896× at 1.01×).
- **The AKT–MKRN1–PTEN loop is railed** (specs/082).

So prediction 1 is now expected to fail. The prediction stands as registered.
