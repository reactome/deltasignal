# specs/134: a total-protein antibody reads how much protein there is, not which forms it is in

**Status:** measured 2026-11-04 (research.md R1–R4). The korkut gate fails (balanced +0.0024 [−0.0067, +0.0104];
accuracy +0.049); the CPPA training check (Q6) holds. Not adopted under the registration; adoption is Adam's.
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
