# 022 — Self-contained inhibition: the input sets the direction, the inhibitor only damps it

## Motivation

specs/012 found that an inhibitor which *contains* its own reaction's input
(a sequestering complex such as WIF1:WNT on *WNT binds to FZD and LRP5/6*)
double-counts that input under `divide` inhibition:
- **one** such inhibitor cancels the input exactly (x · 1/x = 1);
- **two** invert it, giving x · min(x⁻², 10) (specs/021, WNT5A trace).

Deleting those edges was measured at −61 held-out, but that arm was run on a
hash-seed-dependent flag set and is not reproducible (specs/012 addendum).

Adam, 2026-09-25: *"when we have a direct inhibition from an input I feel like
the inhibition should be a lot less than the input. so if the input is up the
output of the reaction should just be a little less up and if it is down
regulated it should be down regulated but just a little less."* And: *"later on
we could learn from data the actual weight it should be."*

## The rule (`DS_SELF_INHIBITOR_WEIGHT=w`, off by default)

An inhibitor I of a reaction is **self-contained** when Reactome's containment
table says I contains the stable id of one of that reaction's activators (the
shared inputs S). Only `divide`-mode inhibitor edges are affected; depletion
edges are not.

I's fold f_I is split into two parts:

- **the part explained by the shared inputs**, f_S = ∏ fold(s ∈ S). This part
  is kept only at a small power w / n, where n is the number of self-contained
  inhibitors on the reaction, so that the reaction's total damping exponent is
  w however many there are;
- **the independent part**, f_I / f_S. This keeps full strength, so a WIF1
  overexpression still inhibits exactly as today.

    effective fold of I = (f_I / f_S) · f_S^(w/n)          (f_S → 0: independent part = 1)

With every self-contained inhibitor tracking its input exactly, the reaction
reads x^(1−w): 2x in gives 2^(0.9) = 1.87x out, and 0.5x in gives 0.54x out.
At baseline the rule changes nothing. w is a single structural parameter. It is
**fixed at 0.1 before the run**, not tuned on the evaluation set, and is the
first candidate for learning from data later.

## Pre-registration (committed before any arm runs)

**Arms.** Both run on build `20260925-1039_d4f4f64`, code defaults otherwise.
- **Control:** `DS_SELF_INHIBITOR_WEIGHT` unset. It must reproduce the
  production scoring (`ae84de9`) byte-for-byte, or the arm is not interpreted.
- **Arm:** w = 0.1.

**Decides:**
- curator **held-out** net cases (fixed − broken) and macro-F1, with an exact
  McNemar p;
- how many pathways moved out of how many were scored, and how many distinct
  readouts the discordant cases span.

**Reported alongside:** tuning split, experimental axis, per-pathway net, and
the WNT5A readout's dose ladder.

**Reading rule, fixed now:**
- **Adopt as the default** only if held-out net > 0 with p < 0.05, the gain is
  not concentrated in one pathway or one readout, and the experimental axis is
  not worse by more than the noise floor (15 cases).
- **Record as neutral** if held-out |net| is within the regeneration noise
  floor (15).
- **Record as negative** if it is below −15.

**Prediction:** a small effect. Only 26 distinct reactions (174 variant nodes)
carry two or more self-contained inhibitors, and 628 pairs have at least one.
Most of them are in PIP3, which is a tuning pathway.

## Result (build `20260925-1039_d4f4f64`, solver `ffb3aa1`, pinned worktree, two containers)

**The control is byte-identical to production** (`ae84de9`) on both axes. The
probe confirmed the arm damped 43 inhibitor slots in Signaling by WNT, and the
control damped 0. Every case is present in both arms.

| curator | control | w = 0.1 | net | fixed / broke | McNemar p | pathways moved | genes |
|---|---|---|---|---|---|---|---|
| **held-out** | 0.8687 / mF1 0.8299 | 0.8701 / mF1 0.8327 | **+27** | 62 / 35 | 0.008 | **3 / 71** | 9 |
| tuning | 0.7694 / 0.7612 | 0.7741 / 0.7662 | +24 | 25 / 1 | < 1e-4 | 3 / 11 | 6 |
| all | 0.8477 / 0.8145 | 0.8498 / 0.8178 | +51 | 87 / 36 | < 1e-4 | 6 / 82 | |

Experimental axis: +7 (8 fixed, 1 broke, p = 0.039). Of the 8 fixes, 8 are
Cell Cycle Checkpoints, a tuning pathway.

Per pathway (curator):

| pathway | split | net | perturbations |
|---|---|---|---|
| Transcriptional regulation by RUNX1 | held-out | **+27** (30 / 3) | **MIR675 only**, KD and OE |
| Cell Cycle Checkpoints | tuning | +21 (21 / 0) | ATM, ATR, CHEK2 |
| Pre-NOTCH Expression and Processing | held-out | +11 (11 / 0) | CCND1, E2F1, JUN, NOTCH1 |
| Signaling by WNT | tuning | +2 | WNT5A OE |
| Signaling by ERBB2 | tuning | +1 | ERBB2 OE |
| Transcriptional regulation by RUNX2 | held-out | −11 (21 / 32) | CBFB, ESR1, PPM1D, RUNX2 |

**Verdict under the pre-registered rule: NOT adopted as the default.**
- Held-out net is positive and significant, but it is **concentrated**: one
  gene in one pathway (MIR675 in RUNX1, +27) equals the entire held-out net.
- Without RUNX1 the held-out split is 32 fixed against 32 broken, which is
  exactly zero.
- The reading rule required the gain not to be concentrated in one pathway, so
  it fails that condition.
- The experimental gain is also one pathway (a tuning one).
- Nothing got worse beyond the noise floor on either axis. The rule is safe,
  and it is a real correction of an operator that double-counts, but on this
  catalog it is not a demonstrated accuracy improvement.

**The WNT5A reversal that motivated it is NOT fixed.** Under the arm, the
readout still goes 100, 100, 1.9e-5, 1.6e-5. On 12 of the 16 variants of *WNT
binds to FZD and LRP5/6*, the second inhibitor (WNT3A:sFRP) does not contain
the variant's ligand. It is correctly left at full strength, but it falls
anyway, because the benchmark's set pin knocks WNT3A down too. On a WNT8A
variant the output is then x · x^(−0.1) · x^(−1), which reads slightly UP at
KD 0.5; at KD 0.05 the ceiling binds and it reads DOWN. The loop rails both.
This confirms the corrected specs/021 trace: that case is caused by the pinning
rule, which this change does not touch.

**Not yet traced:** why MIR675 moves 30 RUNX1 readouts. One plausible reading
is a miRNA:mRNA complex that inhibits translation of an mRNA it contains, which
is exactly the self-contained structure. Until one case is traced, the +27 is
not claimed as that mechanism.

**What stays:** the flag, off by default, with its tests, as a learnable
structural weight. It is not adopted: the evidence is one gene. The next
evidence would come from learning w on data, as Adam proposed, rather than
from sweeping it on this benchmark.

## Re-measured under root pinning (specs/023): ADOPTED

With the protocol of record, w = 0.1 on top of `DS_ASSEMBLY_LIMITING=0` gives:
- curator held-out **+90** (127 / 37, p 1e-12), 10 of 11 moved pathways up,
  25 perturbations;
- **+57** without the top pathway (RUNX1);
- tuning +37; experimental +9.

It meets every condition of the pre-registered reading rule, so it is the
default. The broad-pin result above (+27, one gene) was real but masked. The
pins set the complexes whose inhibitors this rule corrects. Full table in
specs/023.
