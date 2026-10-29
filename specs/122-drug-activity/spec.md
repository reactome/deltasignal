# specs/122: a drug blocks its target's activity, not its abundance

**Status:** registered 2026-10-28, before any code or arm. Adam, 2026-10-28, on the specs/121 trace: "we need to
change delta signal then. we need to know why this is happening", then "ok" to this plan. The freeze (specs/111) is
lifted for this change by that decision.

**Why.** Every drug bench (korkut, decryptM, HPN-DREAM, CPPA) pins a drug's target gene to 0, as a genetic knockout.
An inhibitor does not remove its target; it stops what the target does. Under a knockout pin the target's own
phosphorylated forms vanish, so the model predicts them falling hard; a kinase inhibitor typically leaves them, or
raises them when the target's feedback onto its own activators is lost (MEK inhibitors raise p-MEK). A trace of a MEK
inhibitor on the CPPA training portion (specs/121; its values stay outside this repository under the CPPA terms) found
exactly this, and found that the targets' own readouts cost a large share of the correlation. Genetic perturbations
(the curator and experimental benchmarks) are knockouts and are not touched.

## The rule (`DS_DRUG_ACTION=activity`; `knockout`, the default, is today's behaviour)

For a drug with target genes G, in each network:
- **Nothing is pinned** for the drug. The target's forms are computed as any node is.
- **Its activity is cut.** A `catalyst` or `regulator` edge (positive or negative) from node S into reaction R is
  redirected to a new source node held at 0 when (a) S's member leaves include a protein of G, and (b) the
  catalyst activity or regulation of R that S stands for has an active unit containing a protein of G, or has no active
  unit (Reactome 97: `catalystActivity`/`regulatedBy` → `activeUnit`). So blocking MEK does not block RAF in a RAF:MEK
  complex whose active unit is RAF. Redirected, not dropped: a dropped catalyst edge would let the reaction run as if
  it needed none. A negative regulator held at 0 de-represses its target, as the drug removing it would.
- **Reachability** (the comparator) is unchanged: routes from the nodes holding G.
Stimuli (HPN-DREAM ligands, raised) are unchanged.

## Measurement (FIT data and the CPPA training portion only; CPPA test portions untouched)

Arms `knockout` (control) and `activity` on the top-level build `20261027-1346_8b6527e_top114`, every other setting the
default: korkut2015, decryptM (kinobeads map), CPPA training portion; HPN-DREAM on its own build.
- **Primary (A1):** CPPA training per-drug correlation (fold arm, specs/121's measure): activity − knockout, 95%
  bootstrap over drugs (seed 122), interval above 0.
- **No-harm (A2):** korkut and decryptM balanced accuracy, and HPN-DREAM mean AUROC: no interval for activity −
  knockout entirely below 0 (korkut and decryptM by held-out drug as specs/117; HPN-DREAM over its 30 contexts).
- Reported: the same at CPPA ≤ 12 h; the target-own readouts under the MEK inhibitors traced again; convergence.
CPPA numbers are reported in the private analysis; this repository records the verdicts and the aggregates for the
public datasets only.
**Adoption** (Adam's decision) needs A1 and A2. If adopted, `activity` becomes the drug benches' default and every
specs/121 arm is scored under it at the end; the CPPA test portions stay unread until then.
