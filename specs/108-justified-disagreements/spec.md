# specs/108: justified disagreements

Adam, 2026-10-25: a justified disagreement is one where there is a biological reason the loop is
acting as we would expect the biology to act. This spec makes that a test that can be applied
mechanically and written down before a case is looked at, so it cannot become a way to excuse
any loss after the fact.

## Definition

A case where a rule's prediction differs from the ground truth (curator or experimental) is a
**justified disagreement** only if all four hold:

1. **Traced.** The change is traced case by case to one network element (a carrier, an edge, a
   pool), by adding or removing that element alone and re-solving. A plausible pathway-level story
   does not count.
2. **Independent expectation.** The rule's behaviour at that element matches a biological
   expectation stated without the case's label: a physical constraint (a recycled chaperone's
   total is conserved), Reactome's own curation (read in Neo4j), or a cited source.
   "The curators must be wrong" is not one.
3. **The old answer was an artefact.** The answer the rule changed came through something the
   project already counts as a fault: a rail (B3 below), the iteration cap, label dependence, or a
   double count (specs/012, 022). "Right for the wrong reason."
4. **The new answer still fits the biology.** It agrees in direction with the expectation in (2),
   and does not merely remove the artefact.

**Symmetry.** The same test applies to gains: a case a rule fixes, or the defaults get right,
through an artefact in (3) is flagged **not credible**. Both lists are reported beside every score.

**Scores stay raw** (specs/096). Justified disagreements and not-credible cases are listed, never
subtracted. Adoption decisions may weigh them; adoption is Adam's.

## Rail, operationally

A case is **rail-dependent** when the same perturbation at a mild dose drives its readout to a rail:
re-run with a knockout pinned at 0.5× instead of 0 (`DS_PERTURB_UI_DOWN=0.5`) and an overexpression
at 2× instead of 80× (`DS_PERTURB_UI_UP=2`), and the readout reads ≤ 0.01× or ≥ 99× (the specs/107
B3 thresholds). A halving that empties a readout, or a doubling that saturates it, is the loop
amplifying, not the biology.

## Audit (pre-registered, run once)

Arm `mild108` on the canonical build `20261024-2211_7f4b28a`, solver defaults, bench overrides
`DS_PERTURB_UI_DOWN=0.5 DS_PERTURB_UI_UP=2`. For every valid case of `dflt102` (curator and
experimental), the mild arm's readout value decides rail-dependence. Reported:
- how many cases the defaults get right that are rail-dependent (not credible), per axis, split
  and pathway;
- the same count for wrong cases (rails that cost us);
- the ten specs/107 amendment 1 cases (Class I MHC knockouts) under the four criteria.

**Predictions:**
- Q1: at least 2% of the defaults' correct held-out curator cases are rail-dependent.
- Q2: rail-dependent correct cases concentrate in fewer than 15 pathways.
- Q3: the ten Class I MHC knockouts are rail-dependent under the defaults.

## Applying it to the specs/107 rule 3b cases (recorded before the audit)

- (1) Traced: yes, BAP31 oligomer alone moves all ten (`specs/107-loop-rules/mhc_iso.jl`).
- (2) Independent: Reactome curates BAP31 oligomer joining the peptide loading complex
  (R-HSA-983142), riding peptide loading (R-HSA-8951499) and dissociation (R-HSA-983161), and
  released unchanged by transport of the MHC heterotrimer to the ER exit site (R-HSA-983138):
  a closed cycle, so its total is conserved and it cannot amplify a ligase knockout.
- (3) Artefact: decided by the audit (Q3), not asserted here.
- (4) Still fits: under rule 3b the knockouts lower both readouts (0.75–0.92×), the curators'
  direction.
