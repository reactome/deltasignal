# specs/108 results

Arm `mild108`: canonical build `20261024-2211_7f4b28a`, solver 00f0283 (defaults), bench
`DS_PERTURB_UI_DOWN=0.5 DS_PERTURB_UI_UP=2` (the log's `Protocol:` line confirms both). Compared case
by case with `dflt102` (`rail_audit.py`); a correct case is **not credible** when its mild-dose
readout is ≤ 0.01× or ≥ 99×.

## Audit

| axis | split | right | rail-dependent (not credible) | wrong | rail-dependent |
|---|---|---|---|---|---|
| curator | held-out | 16,778 | 1,343 (8.0%) | 1,703 | 134 (7.9%) |
| curator | tuning | 4,409 | 275 (6.2%) | 621 | 12 (1.9%) |
| experimental | tuning (all) | 641 | 61 (9.5%) | 204 | 7 (3.4%) |

Curator rail-dependent correct cases sit in 34 pathways (29 held-out). The largest: DNA Double-Strand
Break Repair 513, Transcriptional regulation of pluripotent stem cells 290, HDR 133, Interferon α/β
104, Mitotic G2-G2/M 64, DNA Double Strand Break Response 60, Base Excision Repair 43, PTK6 40.
Experimental: PIP3 28, RAF 20, HDR 13.

**Predictions:** Q1 holds (8.0% ≥ 2%). Q2 fails (34 pathways, predicted fewer than 15; half the
cases are in two). Q3 fails (below).

## The ten specs/107 rule 3b cases

| case | defaults | mild dose | rail-dependent | criteria 1 / 2 / 3 / 4 | justified |
|---|---|---|---|---|---|
| 9 E3 ligase knockouts at Antigen peptide bound class I MHC [plasma membrane] | 0.667× | 0.667× | no | yes / yes / **no** / yes | **no** |
| HLA-A knockout at the same complex [Golgi membrane] | 0.0028× | 0.0019× | yes | yes / yes / yes / yes | **yes** |

The Golgi-membrane readout of the same ligase knockouts is railed (0.0015× at both doses), but those
cases are scored right under both the defaults and rule 3b, so they do not enter. The
plasma-membrane readout the nine losses come from is not railed by the registered test, so under
specs/108 they are ordinary losses, not justified disagreements. Rule 3b's held-out −10 therefore
stands as −9 unjustified, 1 justified.

**Observed, not registered as an artefact:** the ligase knockouts read the same at 0.5× and at 0
(to six digits), and an overexpression at 2× reads the same as at 80× (34×, 100×). A response that
ignores dose is suspect, but specs/108 did not list it as an artefact class, so it does not count
here. Adding it would need its own registration.

## What this means for the scores

About one in thirteen correct curator cases, and one in ten correct experimental cases, is reached
through a loop that drives a mild perturbation to the limit. Scores stay raw (specs/096); these
counts are reported beside them. A rule that removes rails (specs/107 rule 3 and 3b did, by design)
will tend to lose some of these cases, and the audit shows which.
