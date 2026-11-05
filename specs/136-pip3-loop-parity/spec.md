# specs/136: the PIP3 switch — one lipid pool, and PIP3 counted once per reaction

**Status:** registered 2026-11-04, before any code or arm (Adam: "ok"). Follows specs/135 and specs/133's pointer to
"PIP3's AKT readouts under drugs".

## What was seen (control, before any rule)

Control: korkut arm `k134_abund_a2` (code `e4809b8`, frozen build `20261026-1349_8b6527e`). In *PIP3 activates AKT
signaling* (R-HSA-1257604), 1,059 korkut cases are scored; 216 are wrong, 153 of them DOWN calls measured unchanged.
They are the AKT-axis readouts (AKT pS473 / pT308, GSK3 pS9 / pS21, TSC2 pT1462, SRC pY416) under the HDAC, SRC, mTOR
and PI3K drug sets, and TSC2 pT1462 under nearly every drug set.

**Traced (SRC knockout, AKT pS473 reads 0.21).** SRC is in one of about 30 PI3K activator complexes
(ESRs:STRN:ESTG:MyrG-pY419 SRC:PI3K alpha). 12 of about 400 PIP3-producing reaction copies go to 0, so PIP3's
producers read 0.97 on average. PIP3 reads 0.398 because its depletion by PTEN (R-HSA-199420) reads 2.42: the curated
loop PIP3 → AKT → p-AKT phosphorylates MKRN1 (R-HSA-8948757) → p-S109-MKRN1 depletes PTEN → PTEN depletes PIP3 ran to
its low state. Along the AKT activation chain PIP3's fold enters several times: in *mTORC2 phosphorylates AKT at S473*
(R-HSA-198640) through both the input (AKT:PIP3) and the catalyst (0.40 → 0.17); in *AKT binds PDPK1* through a second
PIP3 complex (0.17 → 0.07); in *PDPK1 phosphorylates AKT at T308* (R-HSA-198270) through one complex that is both input
and catalyst (0.07 → 0.005, the specs/110 squaring).

**Loop parity (Adam: "check if it is an odd or even number of negative interactions in a loop").** Every simple cycle
through PIP3 of at most 16 edges in its strongly connected component (425 nodes), by the parity of negative edges
(`depletion` or a negative regulator):

| cycles | negative edges | what they are |
|---|---|---|
| 39 | 2 (even) | the MKRN1 loop: PIP3 → AKT → p-MKRN1 ⊣ PTEN ⊣ PIP3. Positive feedback, as curated |
| 372 | 0 (even) | PIP3 → *PTEN dephosphorylates PIP3* (R-HSA-199456) → PIP2 → *PI3K phosphorylates PIP2 to PIP3* (R-HSA-2316434) → PIP3: a substrate–product interconversion read as a gain-1 positive loop |
| 5,580 | 1 (odd) | MKRN1 ⊣ PTEN, then PTEN's own conversion making PIP2 that PI3K makes into PIP3: PTEN read as raising PIP3 |
| 6 | 3 (odd) | the MKRN1 loop through THEM4/TRIB3's inhibition of AKT phosphorylation |

PTEN reaches PIP3 with both signs: lowering it through the curated depletion, raising it through its product. The sign
the biology has (PTEN lowers PIP3) is carried only by the depletion. The interconversion is the PIP2 ⇄ PIP3 lipid pool
that specs/107 rule 1 (`DS_METABOLITE_POOLS=on`) solves as a conserved quantity, the same kind of loop that rails in
PI Metabolism (specs/133 R5 group 6).

## Arms

All on the frozen build, `run_arm.sh`, axes curator (held-out, tuning), experimental and korkut:

- **P** (`k136_P`): `--server DS_METABOLITE_POOLS=on`, the existing rule 1, unchanged.
- **L** (`k136_L`): `--server DS_SELF_ACTIVATOR_LEAVES=1`, new. specs/083 counts a positive catalyst or regulator once
  when it CONTAINS its reaction's input and is computed from it. L extends the pair test as specs/040 rule B extended
  specs/022 for inhibitors: the activator and the input share a non-cofactor leaf, and a carrier of that leaf reaches
  both (built from the same species). Same weight (`DS_SELF_ACTIVATOR_WEIGHT`, 0.1), same w/n overlap; it only weakens.
  Default off; reports `self_activator_leaves` and `self_activator_leaf_pairs`. An activator that is the input itself
  (specs/110) is not touched.
- **PL** (`k136_PL`): both.

Control for every axis: the same code with both off (`k136_ctl`, which must reproduce `k134_abund_a2` on korkut and the
frozen curator and experimental predictions case for case).

## Gate

korkut, arm minus control on cases both score: three-class balanced accuracy, 95% interval from resampling drug sets
(`score_134.py`); the interval must be above 0, read with RAF set aside, both figures reported. Guard: curator held-out
must not lose significantly (McNemar p < 0.05 against), as specs/107. Tuning and experimental are reported, not gates.
An arm that passes is then run once on the CPPA training portion (no harm on the primary, as specs/134). Adoption is
Adam's.

**Honest limit.** korkut chose the trace, and korkut is FIT; a korkut pass is development evidence only. Curator
held-out is SELECT. A confirmation would be data none of this has touched.

**Reported beside the gate:** the parity census above for P (with the pool's conversion edges removed from the
iterated graph) and the SRC-knockout AKT pS473 value in PIP3 under each arm.

## Predictions (before any code)

- **Q1.** P removes the 372 even and 5,580 odd cycles through the conversion; the 39 + 6 MKRN1 cycles remain.
- **Q2.** P alone does not take the MKRN1 switch off its rail: SRC-knockout AKT pS473 in PIP3 stays below 0.85, and
  fewer than 50 korkut calls move. No curator class change (as measured in specs/107).
- **Q3.** L flags the S473 step (catalyst and input share PIP3) and fewer than 1,000 slots across the catalog. Alone it
  does not unrail the switch either, because the T308 step's squaring (input = catalyst) remains.
- **Q4.** PL moves more PIP3 calls than either; whether it unrails the switch is the open question. Probably the gate
  fails for all three arms (the switch is curated positive feedback, like RAF); a pass, if any, comes from PL.
- **Q5.** Curator held-out: P 0, L within ±20, PL within ±25.
