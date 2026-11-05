# specs/138: korkut's zero calls, and which drugs reached their targets

**Status:** analysis, 2026-11-05; no rule registered. Amendment 1: korkut reported by engagement. R1: the engaged
cases traced; no rule indicated (S6K → IRS1 feedback not curated). Follows specs/137 (the change calls that matter sit at zero, not near
the cutoff; Adam: "continue"). Control: korkut `builds/20261026-1349_8b6527e/results/e4809b8/k134_abund_a2` (the default).

## The zero calls

853 of the 1,788 korkut change calls are predicted at ≤ 0.01 of baseline (a target knocked out and read straight
through), and 608 of them are wrong; 587 of the wrong ones are measured NOT. They concentrate in three drugs' targets:
STAT3 (175), the JAKs (128 alone, more in pairs) and SRC (107). By readout: STAT3 pY705 150, SRC pY416 108, STAT5 pY694
60, then STAT3's transcription targets in *Interleukin-4 and Interleukin-13 signaling* (MYC, BCL-XL, BCL2, cyclin D1,
p21, COX2, about 23 each).

## The ground truth: three of those drugs did not move their own target

Each drug's direct marker on the array, named by mechanism in `bench/analysis/138/engagement.py` (the target's own
activation mark or its direct substrate). Engaged: a marker moves |log2| ≥ log2 1.5 the expected way at either dose.
I had already seen the values for AK, 901, P6, SR and ST when the markers were named; the others were read after.

| drug (target) | marker, log2 at each dose | engaged |
|---|---|---|
| PD0325901 (MEK) | MAPK pT202 −0.85 / −1.09 | yes |
| AKTi 1/2 (AKT) | AKT pS473 −2.84 / −2.53; pT308 −2.06 / −2.00 | yes |
| PLX4720 (BRAF V600E) | MEK pS217 −0.90 / −0.78; MAPK pT202 −0.65 / −0.63 | yes |
| Temsirolimus (mTORC1) | p70S6K pT389 −1.39 / −1.60; S6 pS235 −2.96 / −4.12 | yes |
| ZSTK474 (PI3K) | AKT pS473 −2.16 / −2.29 | yes |
| Nutlin (MDM2) | p53 +0.20 / +0.16; p21 +0.26 / +0.05 | no |
| P6 (JAK) | STAT3 pY705 −0.17 / −0.29; STAT5 pY694 −0.13 / −0.12 | no |
| RY (CDK4) | RB pS807 −0.56 (one dose) | no (just under) |
| SR (SRC) | SRC pY416 −0.29 / −0.24 | no |
| Stattic (STAT3) | STAT3 pY705 −0.22 / −0.07 | no |
| HNHA (HDAC), RO-31-7549 (PKC) | no direct substrate on the array | no marker |

In these cells (SKMEL-133) at these doses, JAK–STAT and SRC carry little signal that the drugs remove, or the drugs were
weak; Nutlin raising neither p53 nor p21 fits a cell without functional p53. The single drugs Stattic, P6 and SR
change 0 of 36, 0 of 37 and 1 of 49 scored readouts. The model removes each target entirely, so it calls about a quarter
of their cases changed.

**503 of the 608 wrong zero calls, and 804 of all 1,408 wrong korkut change calls (57%), come from drug sets whose
removed targets are all Stattic's, P6's or SR's.** 12 right zero calls do.

## Scored by engagement (reported, no gate)

Drug sets split by whether every drug, some, or none engaged; model, always-NOT and reachability on the same cases.
Without each drug set's own markers (they are measured changed by definition when engaged):

| drug sets | cases | model balanced / accuracy | always-NOT accuracy | reachability balanced / accuracy |
|---|---|---|---|---|
| all | 9,476 | 0.403 / 0.854 | 0.938 | 0.396 / 0.825 |
| every drug engaged (15 sets) | 1,571 | 0.439 / 0.899 | 0.912 | 0.442 / 0.863 |
| some engaged (35) | 4,423 | 0.425 / 0.859 | 0.930 | 0.413 / 0.826 |
| none engaged (28) | 3,482 | 0.331 / 0.829 | 0.961 | 0.323 / 0.806 |

With the markers kept, every-drug-engaged reads 0.492 / 0.868 against always-NOT 0.845, the only korkut subset where the
model beats always-NOT, but that win is the markers themselves.

## What this says

- Most korkut errors are not the network: they are drugs that did not act in these cells, which the model, removing
  the target, cannot know. A Reactome pathway does not say whether it is running in a given cell, or how hard a drug
  hits. That needs cell context (expression, a pathway-activity readout, measured potency), not a structural rule.
- Where the drugs did act, the model still does not beat always-NOT on the readouts beyond the drug's own marker, and
  only ties reachability. That is the remaining question for korkut, and a smaller one (1,571 cases).
- Engagement is read from the same measurements that are scored, so splitting by it is a description, not a test. A
  scoring protocol that reports by engagement would apply to every method alike, but it is a protocol choice: Adam's.
- specs/120 (fitted inhibition depth) and specs/122–123 (activity vs knockout per drug) asked how a drug acts; this asks
  whether it acted at all here. Neither of those can turn a non-engaged drug into "no change".

## Amendment 1 (2026-11-05): korkut is reported by engagement, beside the raw score

Adam ("ok", to: report both, raw stays the headline). `bench/analysis/138/engagement_report.py CASES` gives the raw score
and the three splits for the model, always-NOT and reachability. Fixed with it:

- the marker list in `engagement.py` is frozen as of this commit;
- each drug set's own markers are left out of every split row;
- the split is never used to choose or gate a rule; gates stay on the raw cases as registered.

## Amendment 2 (2026-11-05): a drug without a marker is unknown, not "not engaged"

**Found by code review** (2026-11-05). `split_of` filed every drug that is not "yes" as not engaged, so HNHA and
RO-31-7549 (no direct substrate on the array, engagement unknown) put 13 of the 28 "none engaged" and 10 of the 35
"some engaged" drug sets in the wrong rows; no text said so. The report now has a fourth row: a set is filed by its
drugs whose engagement is known only when the unknown ones cannot change the answer (an engaged and a non-engaged drug
together are "some engaged" whatever HNHA or RO do); otherwise it is "engagement unknown". The marker list is unchanged.

Re-run on the same cases as the table above (`k134_abund_a2`), each drug set's own markers left out:

| split | cases | model (balanced / accuracy) | always-NOT accuracy | reachability (balanced / accuracy) |
|---|---|---|---|---|
| every drug engaged (15 sets) | 1,571 | 0.439 / 0.899 | 0.912 | 0.442 / 0.863 |
| some engaged (25) | 3,178 | 0.416 / 0.847 | 0.926 | 0.405 / 0.812 |
| none engaged (15) | 1,796 | 0.485 / 0.818 | 0.970 | 0.479 / 0.790 |
| engagement unknown (23) | 2,931 | 0.391 / 0.861 | 0.947 | 0.380 / 0.838 |

The earlier "none engaged" balanced 0.331 came mostly from the HDAC and PKC sets. What stands: the model is below
always-NOT in every row, and the every-drug-engaged row (unchanged) does not beat always-NOT beyond the drugs' own
markers. What changes: on sets whose drugs are known not to have acted, the model's balanced accuracy is not lower than
elsewhere; its errors there are change calls on readouts that do not move (always-NOT is 0.970). Five drugs are classed
not engaged (Nutlin, P6, RY, SR, Stattic); three of them (P6, SR, Stattic) carry the error concentration above.
