# specs/123: a drug acts by its mechanism (knockout or activity, per drug)

**Status:** registered 2026-10-30, before the mechanism labels exist, before any code and before any arm. Adam,
2026-10-30, on the specs/122 result: "shouldn't we use it when it works and not other times", then "I agree with your
plan" to choosing per drug by documented mechanism rather than per dataset.

**Why.** specs/122 found that neither drug rule fits every drug. Modelling a drug as blocking its target's activity
(`activity`) gains on CPPA's training drugs and loses on korkut and HPN-DREAM; a decomposition of the CPPA gain (private
note) puts almost all of it, and of the loss, in the drugs' own target readouts. Which way a target's own phosphorylation
goes is known pharmacology:
- an inhibitor that stops its target being activated (the PH-domain allosteric AKT inhibitors MK-2206 and AKTi 1/2,
  which keep AKT off the membrane) or removes it (a degrader) lowers the target's activated forms, as a knockout pin does;
- an inhibitor of the target's catalysis (ATP-competitive kinase inhibitors; allosteric MEK inhibitors that leave RAF's
  phosphorylation of MEK) leaves those forms, and often raises them by relieving feedback, as `activity` models.
The same allosteric AKT inhibitor loses under `activity` in CPPA and in korkut, so the datasets do not disagree; their
drug mixes differ. Choosing the rule per dataset would be fitted to the scores. Choosing it per drug from its documented
mechanism, fixed before scoring, is not.

## Change

- **Labels** (one per drug, on its primary target, from the literature): `knockout` (K1 degraders; K2 inhibitors that
  block the target's own activation step: membrane recruitment, its upstream phosphorylation, activating dimerisation or
  SH2 engagement; K3 antibodies or antagonists blocking a receptor's ligand binding or dimerisation), `activity`
  (catalytic-site, ATP-competitive or allosteric inhibitors of what the active target does, effector-binding PPI
  inhibitors, bromodomain inhibitors, rapalogs), or `unknown` (an anonymised compound such as CPPA's "AKTINHIBITOR01",
  or an undocumented mechanism; identities are never guessed from a class name).
- **Who labels, and blind to what.** A separate agent labels every drug of korkut, decryptM, HPN-DREAM and CPPA from
  public literature, with a source for each, reading only a list of drug names and target genes: no prediction, score
  or measurement. Disclosed: I (the author of this spec) have seen per-drug CPPA training correlations under both pure
  rules (specs/122) and wrote the rule above with that knowledge; the per-drug labels are not mine. Labels are committed
  (public datasets in `bench/drug_mechanism.tsv`; CPPA privately) before any arm, and are not edited after scoring.
- **`DS_DRUG_ACTION=mechanism`** (bench-side): each drug's target genes are perturbed by its label (`knockout`: root
  pins at 0 as today; `activity`: the specs/122 edge cut); `unknown` uses `DS_DRUG_MECHANISM_UNKNOWN` (default
  `knockout`, today's rule). A drug set (korkut) perturbs each member by its own label in one solve. A drug missing
  from the table is an error.

## Measures (as specs/122, against the `knockout` control already scored there)

- **Primary (M1):** CPPA training per-drug correlation (fold arm, specs/121's measure): mechanism − knockout, 95%
  bootstrap over drugs (seed 123), interval above 0.
- **No-harm (M2):** korkut and decryptM balanced accuracy (by held-out drug, as specs/117) and HPN-DREAM mean AUROC over
  its 30 official contexts: no interval for mechanism − knockout entirely below 0.
- **Reported:** mechanism − activity on every axis; the same with `DS_DRUG_MECHANISM_UNKNOWN=activity`; per-label drug
  counts; M1 restricted to drugs with a known label; convergence.
- **Traced, not gated:** the five CPPA training drugs whose `activity` loss specs/122's mechanism reading did not explain
  (palbociclib and AMG900 on downstream readouts; dasatinib, fulvestrant and topotecan on their own targets).

M1 is partly in-sample (the training drugs informed the rule, not the labels). **The test of record is the CPPA test
portions (unseen drugs, unseen cells, both), read once at the end of specs/121.** CPPA numbers stay in the private
analysis; this repository records verdicts and public-dataset aggregates.

**Adoption** (Adam's decision) needs M1 and M2. If adopted, `mechanism` becomes the drug benches' default and the
specs/121 phase-3 learning starts from it.

## Labels frozen (2026-10-30, before any arm)

The labelling agent read only the drug list and public literature. Public datasets: `bench/drug_mechanism.tsv`
(korkut 8 activity / 2 knockout / 2 unknown; decryptM 14 / 4 / 0; HPN-DREAM's AZD8055 activity). CPPA: a private table
(sha256 `3f97b5b7e462…`), 114 drugs, 45 unknown (mostly anonymised compounds such as the eight "MEKINHIBITOR" entries).
Within the MEK inhibitors the literature splits: trametinib blocks RAF's phosphorylation of MEK (knockout); PD0325901,
selumetinib, refametinib and cobimetinib do not (activity). Five labels are low-confidence (chloroquine, chloropyramine,
mitoxantrone, podophyllotoxin, AUY); they stand as assigned.
