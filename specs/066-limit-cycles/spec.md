# specs/066: limit cycles — find the fixed point when relaxation cannot

**Status:** pre-registered 2026-10-16, before the arm. From the week review
(docs/reviews/2026-10-15-week-review.md, section 5): 23 curator held-out cases
change class with the iteration cap, all from five perturbations
(RNF8 ×15 and RNF168 ×3 in DSB Repair / DSB Response, FANCA and FANCC in
Fanconi Anemia, SH3GL1 and ARHGEF7 in EGFR).

## Diagnosis (traced 2026-10-16, canonical `20261016-0006_da08370`, solver `a0f6671`)

Probes are in `~/deltasignal-catalogs/probe_osc/`.

1. **These are limit cycles, not slow convergence.** FANCA overexpression,
   solved at every cap from 1 to 1,200: the Fanconi loop has period exactly
   30 iterations (25 iterations near 100×, 5 near 0), so the scored class is
   whichever phase the cap lands on.
2. **Damping does not cure it.** `DS_SCC_DAMPING` 0.5 / 0.25 / 0.1 / 0.05 /
   0.02 (caps up to 10,000), and damping in log space (x ← x^(1−λ)·F(x)^λ,
   λ 0.2 to 0.02, scratch build): no component converges except EGFR under
   SH3GL1 knockout, and the class still moves with λ.
3. **Why: a Goodwin ring.** Traced per iteration in Fanconi (18-node
   component): FANCA (pinned) and nine free FA core subunits AND-assemble
   into FA Core Complex:ICL-DNA (R-HSA-6785124); five steps on, p-FA core
   complex (R-HSA-6788398) has a `depletion` edge back to each of the nine
   free subunits (and UBE2T). A ring of about six first-order lags closed by
   negative feedback is unstable once the loop gain exceeds sec(π/6)^6 ≈ 2.4;
   here the AND multiplies nine depletions by the same complex, a log gain
   near 9. The fixed point exists (F is continuous on a compact box) but is
   unstable even in the continuous-time limit, so no relaxation step reaches
   it.
4. **Where the edges come from.** All three pathways have the same generator
   shape: the ubiquitin pass of `_emit_substrate_depletion_edges` takes as
   substrate every input leaf whose gene also appears in the output, and the
   catalyst of these reactions is the input complex itself (modified in
   place). So a complex depletes the free forms of its own partners:
   FA core subunits and UBE2T (Fanconi, "Monoubiquitination of
   FANCD2:FANCI", R-HSA-6785361, whose input and catalyst are both
   R-HSA-6788398); KDM4A/B, RNF168, PIAS4 ("RNF8 and RNF168 ubiquitinate
   KDM4A,B" / "... H2AFX", DSB); EGF, HBEGF, BTC, EPGN, EPN1, EPS15L1 (CBL
   ubiquitination reactions, EGFR). This is the specs/057 class; its
   generator fix (`LNG_UB_SUBSTRATE_CHANGED`) was measured and not adopted
   (held-out −39, −15 with the Ub-X test). This spec does not revisit it.
5. **A root finder reaches the fixed point.** Scratch build: when the
   damped iteration ends above tolerance, re-solve the component by
   Levenberg–Marquardt on the pure residual |F(x) − x|² (the existing
   `minimize_component!`, mu = 1, gamma = 0), warm-started from the last
   iterate:
   - Fanconi (18 nodes): residual 1.7e-10, converged; readout 6785124 reads
     2.0× at caps 500 and 2,000 (cap 500 alone read 0.0066).
   - EGFR (373 nodes): with 400 LM iterations, residual ~1e-11, converged;
     every SH3GL1 and ARHGEF7 readout reads DOWN.
   - DSB Repair (454 nodes): residual 3e-4 after 400 LM iterations, not
     converged.

## Rule

`DS_NONCONV_FALLBACK=off` (default) | `lm`. Under `lm`, a cyclic component
whose damped fixed point ends with residual ≥ tolerance is re-solved by
Levenberg–Marquardt on |F(x) − x|² (mu = 1, gamma = 0, no prior), started
from the last iterate, with `DS_NONCONV_LM_ITERS` (default 400) iterations.
The LM state is kept only if its residual is lower than the iterate's;
otherwise the iterate is restored. Components that converged are never
touched, so every converged solve is bit-identical to `off`. The solve
reports `nonconv_fallback` (rule), `nonconv_components`, `nonconv_rescued`
(residual now below tolerance) and `nonconv_improved`.

Not part of the rule: no change to damping, the cap, or the network.

## Arm and measurement

Arm `nclm`: `run_arm.sh nclm --server DS_NONCONV_FALLBACK=lm` on canonical
`20261016-0006_da08370`, against `results/0406a4b` (the arm's solver differs
from `0406a4b` only by this rule, default off). Report:
- held-out, tuning and experimental net with McNemar p, with and without RAF
  (`holdout_report.py --compare B --set-aside RAF_MAP_kinase_cascade`);
- solves converged, before and after;
- every held-out class change, labelled by whether its component was
  rescued; the 23 cap-decided cases separately;
- the per-pathway net and the distinct readouts that move.

## Predictions (from the probe)

- FANCA and FANCC at readout 6785124: DOWN → UP (2 held-out fixed).
- SH3GL1 knockout at readout 167679: UP → DOWN (1 held-out broken).
- ARHGEF7: every readout stays DOWN (no change).
- RNF8 / RNF168 (DSB): mostly unchanged; their component is not rescued.
- Elsewhere (≈1,400 non-converged solves, 907 of them in tuning, mostly
  TP53): no prediction. Converged solves cannot move.

## Decision rule

The rule makes the scored value the model's fixed point instead of the
cap's choice, so it is judged as a faithfulness fix, with the specs/048
gate:
- **Adopt** if curator held-out net ≥ −15 and experimental net ≥ −15, both
  read with RAF set aside, and the number of converged solves rises.
- **Do not adopt** otherwise, and record which pathways carry the loss.
- Either way, report both readings (with and without RAF), and trace every
  held-out class change in a component the fallback did not rescue (those are
  still the cap's choice, so a gain or loss there is not credited to the
  rule).

## Amendment 1 (2026-10-16, before the arm): two implementation details

- **Components holding cycle-balance pool nodes are skipped** and counted
  (`nonconv_skipped_pools`). Their pool nodes are set by the pool's own
  π = πP solve, not by F, and the least-squares step would overwrite them
  with reaction outputs. None of the five traced components holds a pool
  node, so the predictions above stand.
- `DS_NONCONV_LM_ITERS` (default 400) sets the LM budget; the existing
  `DS_LM_ITERS` (60, `DS_SCC_METHOD=minimize`) is not reused, because 60 left
  EGFR's 373-node component at residual ~5e-5.
- Test: `test/test_limit_cycles.jl` (31 assertions). A nine-subunit Goodwin
  ring cycles under `off` (the readout differs between caps 500 and 501) and,
  under `lm`, lands on the analytic fixed point d = P^(1/10) at both caps;
  an acyclic network and a gain-1 ring are bit-identical under both modes.
