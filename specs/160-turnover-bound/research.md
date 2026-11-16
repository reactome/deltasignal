# specs/160 research: the measurement (2026-11-16)

Offline on the recorded 98 korkut TSV, scorer `bench/analysis/160/score_160.py`
(fixture-checked; table and scorer committed in 7f95d7d BEFORE the record was
read). 10,469 scored rows; 6,904 on total readouts; 27 of 44 total readouts
covered by the half-life table.

| | accuracy | balanced | Δ balanced [95% boot, 78 drug sets] | calls moved |
|---|---|---|---|---|
| recorded | 0.8117 | 0.4469 | — | — |
| per-protein bound | 0.8145 | 0.4479 | **+0.0010 [+0.0007, +0.0013]** | 29: 29 fixed, 0 broken |
| tier-1 control (46 h, all) | 0.8204 | 0.4476 | +0.0008 [−0.0030, +0.0033] | 97: 94 fixed, 3 broken |

RAF set aside: per-protein +0.0011 [+0.0008, +0.0014].

## Verdict: the registered gate passes, and the effect is small

The interval excludes zero, and the per-protein delta exceeds the control's,
which are the two registered conditions. The honest reading:

1. **The motivating class is untouched.** The 232 railed total-readout rows are
   right 7 times before and 7 after. The two readouts long-lived enough for a
   complete drop to flip (BAK, PKCα; floors 0.90 and 0.88) carry no railed
   false alarms. The railed read-through false alarms sit on proteins whose
   floors are below 0.85 even with per-protein turnover. Turnover kinetics
   does not explain them: these proteins could physically fall that far in
   24 h, yet the measurement says they did not. The cause is upstream (the
   drug's actual effect on the protein's production), not the readout.
2. **What moves is moderate predicted drops.** A predicted 0.6× on a protein
   with a 0.8 floor reads 0.92× at 24 h, so the call goes to "no change". All
   29 per-protein moves were fixes, and none broke.
3. **A pre-registration expectation was wrong, recorded as such.** spec.md
   expected the tier-1 control to be call-neutral. It is neutral on railed
   rows (floor 0.70 < 0.85), but it moves moderate drops too (97 calls, 94
   fixed). Its balanced-accuracy interval includes zero, and it also applies
   to the 17 uncovered readouts, so the two arms differ in coverage as well as
   in turnover information. The comparison the gate registered still holds;
   the expectation about the control did not.

## Serving

spec.md US3 registered that a pass implements the rule behind
`DS_TURNOVER_BOUND` (default off), with adoption as Adam's call. The effect
is +0.001 balanced accuracy over 29 of 10,469 calls, and the railed class is
untouched. So the flag is NOT built speculatively: building it only enables
an adoption, and that adoption is Adam's decision, put to him with these
numbers first. This deviation from the registered step is recorded here.
