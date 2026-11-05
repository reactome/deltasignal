# specs/135: what is left of the total-protein errors, and why no rule is chosen

**Status:** analysis, 2026-11-04. No rule registered or adopted (Adam: "ok", after the findings below). Follows
specs/134 amendment 2 (R8), whose korkut arm `k134_abund_a2` (code `e4809b8`, frozen build `20261026-1349_8b6527e`)
is the input here. korkut2015 is FIT (specs/062), so everything below is development, not a test.

## What is left

`bench/analysis/census_135.py` classes each scored total readout by its abundance copies: `starting` (no in-edge but
depletion), `made` (produced by a reaction) or `other` (an in-edge from something else).

| copies | cases | change calls | right |
|---|---|---|---|
| starting only | 5,821 | 0 | 0 |
| starting and other | 114 | 0 | 0 |
| made only | 742 | 409 | 10 |
| made and starting | 185 | 105 | 3 |
| made, other and starting | 42 | 13 | 0 |

Every one of the 527 total change calls left (13 right) is on a readout with a made copy, and in Neo4j every producing
reaction of those copies is a gene-expression event: a `BlackBoxEvent` from a gene (`ReferenceDNASequence`), or in one
case translation from an mRNA (PCBP4 and CDKN1A). The largest are *Expression of STAT3-upregulated cytosolic / nuclear
proteins* and *Expression of BCL2, BCL2L1* (Interleukin-4 and 13 signaling), PCNA and CCNE1 expression (Mitotic G1), and
RUNX3- and TP53-driven transcription. Of the 409 single-copy calls, 196 read exactly 0: the knockout removes the only
positive regulator of the expression event, and the regulator is read as required. BCL2's event is regulated by the
p-Y705-STAT3 dimer through a `PositiveGeneExpressionRegulation`, which the generator writes as an AND-joined `regulator`
edge, so a STAT3 or JAK knockout gives no BCL2 at all. korkut measures no change in those cases.

## Two candidate rules, and why neither is registered

Both were considered as the synthesis-stock spec this directory was opened for. Their effect was checked on the 409
single-copy calls, where the readout is the made copy (`census_135.py`, what-if section) before any code was written.

- **A stock.** The korkut lysates are taken 24 h after drug (Korkut et al. 2015, eLife 4:e04640, Methods). At the
  median mammalian protein half-life of 46 h (Schwanhäusser et al. 2011, Nature 473:337), a fraction
  f = 1 − 2^(−24/46) = 0.30 of a protein is new at 24 h, so a made copy reading v would read (1 − f) + f·v. A
  synthesis copy at 0 reads 0.70, still below the bench's 0.85 DOWN cutoff. It moves 39 of 409 calls to NOT.
- **A basal expression level.** Reactome separates `PositiveGeneExpressionRegulation` (299 regulations on 247 reactions
  in the 92 benchmark pathways) from `PositiveRegulation` and `Requirement`, and the two never share a reaction. Read
  as one input beside a basal level at baseline (the existing OR, a mean), a knocked-out regulator gives 0.5 instead of
  0. 0.5 is still DOWN; it also moves 39 of 409 (an upper bound: it assumes every change comes through the regulator).
- **Both.** 61 move; 197 land at 0.848, within 0.01 of the 0.85 cutoff. Whether they read DOWN or NOT would be decided
  by the third digit of the half-life or the hours of treatment, so a measured gain or loss would be the constants'
  coincidence, not evidence for the rule. Not registered.

## What it points to

The bench calls a change at 15% (UI cutoffs 0.85 / 1.15), while korkut's measured calls
need 50% (|log2| ≥ log2 1.5, specs/062). Under matched cutoffs, a knocked-out expression regulator with a 24 h stock
would read NOT with a margin. Changing korkut's call cutoffs after looking at its errors would be post hoc, and the
cutoffs are a protocol choice across datasets: left to Adam, not pursued here.

Also recorded, not pursued: `NegativeGeneExpressionRegulation` (101 regulations) has the mirror question for a
repressor's knockout; per-protein half-lives exist (Mathieson et al. 2018, Nat Commun 9:689) but are a new dataset.
