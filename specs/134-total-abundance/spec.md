# specs/134: a total-protein antibody reads how much protein there is, not which forms it is in

**Status:** adopted 2026-11-04 (Adam) as amendment 1: `DS_TOTAL_READOUT=abundance` is the korkut and CPPA bench
default. The first rule (R1–R4) failed the korkut gate (balanced +0.0024 [−0.0067, +0.0104]); amendment 1 passes it
(+0.0172 [+0.0007, +0.0324], R6) and the CPPA training check holds (R7). `forms` reproduces earlier numbers.
Amendment 2 (a set or pool of the protein's copies carries it) **adopted 2026-11-05** (Adam: "yes"): korkut 144 fixed /
0 broken, exactly the registered cases (R8); CPPA training holds (R9). Amendment 3 (`recycle` is not a made-from edge;
every version keeps its own value) meets its criteria (R10–R11) and is **adopted 2026-11-05** (Adam: "ok"); `abundance_a2` reproduces earlier numbers.
Registered 2026-11-04, before any code or count for the rule. Follows specs/133's conclusion that the
remaining korkut errors belong to single pathways (Adam: "okay go ahead").

## What was seen (control, before the rule)

On the frozen model's korkut cases (specs/133 control `c0`, which reproduces the frozen build), *PIP3 activates AKT
signaling* is right on 521 of 1,083 cases, and 515 of its errors are DOWN where the cells show no change. Most are
total-protein antibodies (AKT, Caspase-9, GSK3, p21, p27, TSC2, SRC, ERα) under drugs that do not touch the protein.
Across all pathways:

| korkut readout kind | scored | measured not | accuracy | model calls a change | of them wrong |
|---|---|---|---|---|---|
| total protein | 7,497 | 7,157 (95.5%) | 0.626 | 2,661 | 2,578 |
| phospho site | 3,565 | 2,892 (81.1%) | 0.680 | 1,261 | 894 |

Traced (PIP3, HDAC1/2/3 knockout, AKT1 readout): the starting copy of AKT1 reads 1.0×, its phosphorylated and bound
forms 0.004–0.5×, and the bench's total readout is their mean, 0.45× (DOWN). The bench reads a total antibody as the
mean over every node that is, or contains, any EWAS of the gene (specs/062), including the gene's DNA and mRNA
(the gene-name match takes every reference entity).

## The rule

A total antibody measures abundance. Binding a partner or adding a phosphate does not change how much of a protein
there is; making it (translation) or destroying it (degradation) does. In the network:

- **P** (part 1): the readout's nodes are those that are, or contain, a *protein* EWAS of its genes (reference entity a
  `ReferenceGeneProduct`, isoforms included). Gene and mRNA nodes are no longer read by a total antibody.
- **Abundance copies** (part 2): a node of P is read only if nothing makes it *from* P: no reaction producing it (an
  `output` edge into it) takes an `input` in P, and no identity edge into it (`specialization`, `variant_pool`,
  `variant_split`) comes from P. Starting copies qualify; a protein translated from its mRNA qualifies; a complex
  formed by binding the protein, or a form made by modifying it, does not. Depletion edges (degradation) are not
  sources, so they still lower a copy they act on.
- The readout is the mean over the unpinned abundance copies, as before over all nodes. None unpinned: `trivial`, as
  before. P non-empty but no abundance copy (a loop with no starting copy): status `no_copy`, not scored.

Phospho readouts are untouched. Nothing in the solver changes; the curator, experimental, decryptM, HPN-DREAM and P100
benches never read a total antibody and are untouched by construction.

Bench setting `DS_TOTAL_READOUT` in `korkut_bench.py` and `cppa_bench.py`: `forms` (default, today's) | `protein`
(part 1 alone) | `abundance` (parts 1 and 2).

## Arms

On the frozen 92-pathway build `20261026-1349_8b6527e`, korkut axis, `scripts/run_arm.sh`:
`k134_forms` (control; must reproduce the frozen korkut predictions case for case), `k134_protein`, `k134_abund`.

On the 28-pathway build `20261027-1346_8b6527e_top114`, CPPA predictions (specs/121 protocol), scored on the
**training portion only** with `score_121.py`: `forms` and `abundance`. Numbers stay private under the CPPA terms;
the public record says only whether each prediction held. The CPPA test portions stay sealed.

## Gate

korkut, on cases scored in both arms, `abundance` minus `forms`: three-class balanced accuracy, 95% interval from
resampling drug sets (10,000, seed 134); the interval must be above 0. Accuracy and macro-F1 reported beside it, and
each also on totals and phospho separately. CPPA training: the registered primary (per-drug correlation, specs/121)
must not fall by more than its own interval half-width, and the secondary balanced accuracy is reported.

**Honest limit.** The korkut errors chose this rule, so a korkut pass cannot confirm it; the CPPA training portion was
not looked at for it and is the check here. A real test is one scoring on the CPPA test portions or Lun 2019, and
opening either is Adam's decision.

## Predictions

- **Q1.** Phospho readouts: every prediction identical in all three arms.
- **Q2.** Under `abundance`, total-protein change calls fall from about 2,660 to under 500, and accuracy on totals rises
  from 0.626 to above 0.90.
- **Q3.** `protein` alone moves little (under 5% of total cases change class): the DNA/mRNA nodes are a small share of
  each readout.
- **Q4.** The gate passes on korkut. Most of the gain is the "not" class; correct DOWN calls on totals (about 76) are
  partly lost, so the gain in balanced accuracy is smaller than the gain in accuracy.
- **Q5.** Total change calls left under `abundance` come from synthesis or degradation (TP53 targets, p21/p27, MDM2
  under nutlin), and fewer than half of them are right.
- **Q6.** CPPA training: secondary balanced accuracy rises; the primary correlation moves by less than 0.01.

Nothing changes a default; adoption is Adam's.

## Amendment 1 (2026-11-04): the rule as written, implemented

**Found by tracing, after R1–R4.** Total AKT under an AKT inhibitor in *Signaling by VEGF* reads 0 under `abundance`,
though every copy it reads has no producing reaction: they are starting complexes of AKT with PIP3 and PDPK1, which the
network builds from AKT's starting copy through `assembly` edges (boundary decomposition, specs/030). The code tested
only `output` edges from reactions and the three identity edges, so a complex assembled from the protein, or a
protein released from a complex holding it (`dissociation`), counted as an abundance copy. The rule above says such a
complex does not ("a complex formed by binding the protein ... does not"). R2's conclusion that every remaining total
change call comes from a synthesis copy is wrong for the same reason: re-counted with the copies classed by their
in-edges, 1,278 of the 1,881 come from readouts whose copies all have no in-edge in the reaction graph (R5).

**The corrected rule** (no other change): a node of P is made from P if a reaction producing it takes an input in P,
or if any material or identity edge into it comes from P: `assembly`, `dissociation`, `recycle`, `released_twin`,
`set_member`, `specialization`, `variant_pool`, `variant_split`. Catalyst, regulator and depletion edges are not
sources, as before. Same setting (`DS_TOTAL_READOUT=abundance`), same statuses.

**Arms:** korkut `k134_abund_a1` against the R1 control `k134_forms` (unchanged code path, so its predictions stand);
CPPA training `abundance` (amendment 1) against `8276c28/cppa_frozen`, numbers private. **Gate:** unchanged.

**Honest limit.** The correction is post hoc on korkut cases, though it implements the registered text rather than a
new idea; korkut still cannot confirm it, and CPPA training has now been looked at once (R4).

**Predictions** (before the arms run):
- **Q7.** Phospho identical again.
- **Q8.** Total change calls fall well below R1's 1,881; those left are mostly (a) synthesis copies and (b) a drug's own
  target, whose copies the knockout pin still reaches (an inhibitor modelled as removal).
- **Q9.** Accuracy rises above R1's 0.700; the gate still fails (balanced accuracy within ±0.01 of forms, interval
  across 0), because the correct DOWN calls on totals were already lost in R1.
- **Q10.** CPPA training: no harm on the primary; secondary balanced accuracy still above forms.

## Decision (Adam, 2026-11-04)

Adopt amendment 1's rule as the korkut and CPPA bench default ("adopt as it is fixing something and the accuracy is
significantly improved and the balanced accuracy is going in the right direction"), conditional on the CPPA training
check (Q10) showing no harm on the primary endpoint. `DS_TOTAL_READOUT=forms` stays, to reproduce earlier numbers.
Headline korkut and CPPA numbers before and after are not comparable (cases leave scoring) and are reported both ways.

## Amendment 2 (2026-11-04): a set or pool of the protein's copies carries it

**Found by tracing, after R7** (specs/135's census of the 695 total change calls left in R6). In *RAF/MAP kinase
cascade* (R-HSA-5673001), total SRC, fibronectin, BRAF and RAF1 each keep one made copy, the activated
RAF:scaffold:MAP2K:MAPK complex (R-HSA-5672720), made by R-HSA-5672972 from the RAF/MAPK scaffolds set and the
activated RAF complexes. Neo4j has p-Y530-SRC, FN1 and both RAFs inside it, through the set's focal-adhesion member and
the RAF dimers. The network gives each set a pool node, fed by `variant_pool` / `set_member` edges from its member
copies; the pool node is not a readout node, so amendment 1's test ("an input in P") never sees that the input holds
the protein. The complex is a bound form and the rule says it is not an abundance copy. 168 change calls come from it,
none right.

**The corrected rule** (no other change): a node *carries* P if it is in P, or is the target of a `set_member`,
`variant_pool`, `variant_split`, `specialization` or `released_twin` edge from a node that carries P (a set, pool or
identity copy of the protein's copies is the protein). A node of P is made from P if a reaction producing it takes an
input that carries P, or if an edge of amendment 1's material or identity set comes from a node that carries P. Only
identity edges spread "carries"; a complex or a released partner never does, so a protein released from a complex is
not treated as its partner.

**Arms:** korkut `k134_abund_a2` against `k134_abund_a1` (R6); CPPA training `abundance` (amendment 2) against the R7
run (`7cb140c/cppa134_abundance`), numbers private.

**How it is judged.** It implements the registered text, as amendment 1 did, so it is judged on faithfulness: every
case that moves must be one whose readout loses a copy built through a set or pool of the protein (checked case by case
by script), and CPPA training must show no harm on the primary. The korkut gate is reported, not applied. Adoption is
Adam's.

**Predictions** (before the arms run):
- **Q11.** Phospho identical.
- **Q12.** The 168 change calls on SRC, fibronectin, BRAF and RAF1 in RAF/MAP kinase cascade read NOT. Fewer than 50
  change calls move anywhere else.
- **Q13.** Accuracy rises by about 0.015; balanced accuracy rises with RAF kept and moves by less than 0.003 with RAF
  set aside (the cases are all RAF).
- **Q14.** CPPA training: no harm on the primary.

## Amendment 3 (2026-11-05): a recycle edge does not make the starting copy

**Found by code review** (2026-11-05, of PRs #216–#223), not by tracing a score. Amendments 1–2 list `recycle` among the
edges that make their target from their source. A recycle edge (specs/081, specs/086) runs from a released or returned
copy INTO the protein's starting copy; the solver strips it before any rule (`src/core/recycle.jl`) and the bench pins
the starting copy as a root. In a release pool the released copy carries the protein (`released_twin`, amendment 2),
so its recycle edge marked the true starting copy as made, against the rule's own "starting copies qualify".

**The corrected rule** (no other change): `recycle` is not a made-from edge.

**Every version keeps its own value** (constitution IV; amendment 2 had changed the adopted rule in place):
`DS_TOTAL_READOUT=abundance_a1` | `abundance_a2` | `abundance_a3`; `abundance` (the default) names the adopted
version, `abundance_a2` until Adam decides. Logs now print the version, not the alias.

**Counted before registration** (structural, no measurement read; `bench/analysis/134/copies_diff.py`): korkut,
frozen build, 1 of 272 (pathway, total readout) pairs differs: b-Catenin in Signaling by WNT (R-HSA-195721) gains its
starting copy (2 → 3 copies). CPPA, top114 build, 3 of 1,293 pairs: CD171 in L1CAM interactions (R-HSA-1266738),
beta-Catenin in Signaling by WNT (R-HSA-162582), FoxO3a in PI3K/AKT signalling (R-HSA-74160), each gaining one copy. No
readout enters or leaves scoring. p130/p107 (the specs/081 recycle pools) have no korkut or CPPA total readout.

**Arms.** korkut: `run_arm.sh --axis korkut --bench DS_TOTAL_READOUT=abundance_a3` (`k134_abund_a3`) against
`k134_abund_a2` (R8). CPPA training: `scripts/run_cppa.sh` (the runner behind R4/R7/R9, now committed) with
`DS_TOTAL_READOUT=abundance_a3` on the three pathways above, spliced into the R9 run for every other pathway (the
readout rule does not change any solve, only which solved nodes a readout averages); numbers private.

**How it is judged.** Faithfulness, as amendments 1–2: every case that moves must be in the pairs above, and CPPA
training must show no harm on the primary. The korkut gate is reported, not applied. Adoption is Adam's.

**Predictions** (before the arms run):
- **Q15.** korkut: phospho identical; only b-Catenin cases in R-HSA-195721 can move; any that change class move toward
  NOT (the starting copy is unpinned under every korkut drug and rests near baseline, diluting the mean).
- **Q16.** CPPA training: only the three pairs above move; no harm on the primary.

**Left open (recorded, not changed here).** A copy whose only in-edge is depletion and that has no out-edge is an
abundance copy but is never pinned (the bench ignores a depletion in-edge only for nodes with an out-edge), so a drug's
own target can be read from that unpinned copy instead of leaving scoring. Not counted yet.

## Decision on amendment 3 (Adam, 2026-11-05)

Adopted ("ok", to the recommendation to adopt on faithfulness with the CPPA secondary's beta-catenin effect recorded).
`abundance` now names `abundance_a3`; `abundance_a2` and `abundance_a1` reproduce the numbers of R8–R9 and R6–R7.
