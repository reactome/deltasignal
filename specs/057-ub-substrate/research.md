# specs/057 research

## B1. The variant build (before any arm)

Build `20261014-1303_73d5e8c_ubs`: LNG `73d5e8c` (branch
`ub-substrate-changed`, LNG PR #113, = `01206ca` + the flag), `--env
LNG_UB_SUBSTRATE_CHANGED=1`, 92 of 92 pathways, status complete. `current`
not moved.

`edge_diff.py` against the canonical build, edge-type counts:

| edge type | change |
|---|---|
| depletion | −307 (in 32 pathways) |
| variant_pool | −356 |
| assembly | −56 (all in Signaling by EGFR, R-HSA-177929) |
| everything else | 0 |

`ub_census.py` on the variant build finds 0 depletion edges onto a
passed-through partner (301 on canonical; the census reads leaves to depth 3,
the generator's own query a little differently, hence 307 vs 301).

- **variant_pool −356**: the pool nodes that gathered a multi-variant
  catalyst's depletion (specs/046) exist only to feed depletions; with the
  depletion gone they are not built.
- **assembly −56, not predicted in the registration.** The generator emits
  depletion edges before it decomposes root complexes, and the hierarchy's
  "not downstream of the root" guard (`_existing_upstream`, LNG #97) follows
  every edge, depletion included. In EGFR, the produced copy of
  `R-HSA-182960` variants counted as downstream of root `R-HSA-182935` only
  through a removed edge (the 182935 pool -| free CBL, R-HSA-112199). Without
  it, the hierarchy joins that produced copy (as specs/030 intends) instead of
  building fresh copies from a free-CBL leaf: 2 fresh copies × 4 variants ×
  … fewer assembly edges. A knock-on of the removed edges through an
  existing rule, not new logic.

## R1. Arm `ubs` against the canonical results: the gate fails

Arm `ubs`: `scripts/run_arm.sh ubs --catalog 20261014-1303_73d5e8c_ubs` at
DS `ad72b42` (solver identical to `64d23cb`), results
`builds/20261014-1303_73d5e8c_ubs/results/ad72b42/ubs`. Control: canonical
`builds/20261012-1401_01206ca/results/64d23cb`. Scored with
`holdout_report.py --cases <control> --compare <ubs> --pair-on-valid
--set-aside RAF_MAP_kinase_cascade`. No case is valid in one arm only.

| axis | net | fixed / broke | p | pathways moved |
|---|---|---|---|---|
| curator held-out | **−39** | 49 / 88 | 0.0011 | 9 / 70 |
| curator tuning | +7 | 10 / 3 | 0.09 | 3 / 11 |
| experimental (RAF set aside identical) | 0 | 1 / 1 (MDM2, one perturbation) | 1 | 1 / 10 |

Headline, in release: held-out 89.33% → 89.12%. Convergence improves:
curator 1,652 → 1,697 of 1,725 solves, experimental 215 → 228 of 244.

**Gate: fails** (held-out net < 0 at p 0.0011; TGF-beta −25, DSB −14 and
EGFR −13 each lose more than 10). Not adopted; the flag stays off.

### Where it moved (held-out unless marked)

| perturbation | net | mechanism |
|---|---|---|
| MET KO, MET OE, HGF KO (MET) | **+30** | as predicted: GRB2 no longer depleted; plain MET now a root (its LRIG1 depletion edge is gone), so it is pinned |
| MIB1 KO (NOTCH1) | +7 | not traced |
| PLCG1, SRC (EGFR) | +6 | not traced |
| FANCA/C/D2 (Fanconi) | +4 | not traced |
| MDM4, AKT1/2, BRCA1 (TP53, tuning) | +5 / −1 | not traced |
| CBL KO+OE (MET) | **−20** | was right for the wrong reason (below) |
| EPS15, SH3GL1 (EGFR) | **−18** | same class |
| SMAD2, SMAD3 (TGF-beta) | **−24** | pin change (below) |
| RNF168 OE (DSB, DSB Response) | **−16** | non-converged in both arms; a coin (below) |
| others | −5 | not traced |

**CBL, EPS15, SH3GL1: negative regulators right only through the wrong
edges.** The curators expect a CBL-complex knockout to raise receptor
signaling (CBL ubiquitination sends the receptor to degradation). In the
canonical network a CBL knockout zeroed the p-Y-CBL pool, which de-repressed
free GRB2 (the removed edge) to 10×, and the readouts rose: right sign,
wrong mechanism. Under the rule CBL depletes only the changed substrate
(ubiquitinated MET/EGFR), and that leaf has no free node in the network, so
CBL has no effect at all and every readout reads 1.0. The negative regulation
Reactome curates (the active receptor complex ubiquitinated, then degraded)
is not represented by either rule.

**SMAD2/3: a defect in the rule, found after the arm.** In canonical every
SMAD2 form had an in-edge (the SMAD2:SMURF2 complex depleting free SMAD2
[cytosol], R-HSA-3006353), so `root_cycle` pinned four recycled forms,
including produced p-S465,S467-SMAD2, which carries the signal. The edge comes
from "SMURF2 ubiquitinates SMAD2" (R-HSA-173542), whose output **Ub-SMAD2
[cytosol] (R-HSA-2176443) is a Complex of SMAD2 [cytosol] and Ub**, not a
modified SMAD2. SMAD2 is the real substrate, but its stId reappears in the
output, so the rule called it passed through and removed the edge (verified in
Neo4j). Plain SMAD2 then became a true root and the only pin (2 copies); in
this pathway it feeds SMURF2 binding and one other reaction and reaches no
readout, so every SMAD2/3 case reads 1.0.

`bench/analysis/057/ub_complex_substrates.cypher`: in 24 of the 301 removed
edges the target sits, in the output, in a complex whose direct components
include Ub (Reactome's "Ub-X" complex representation): TGF-beta 9, R-HSA-1227986
(ERBB2) 5, R-HSA-1980143 (NOTCH1) 3, 7 more. The rule as registered is wrong
for these; a corrected test would count a leaf as the substrate when its stId
changes **or** its direct container in the output also holds Ub.

**RNF168 overexpression: a knife-edge.** The DSB solve under RNF168 at 80×
does not converge in either arm (`dump.jl`, conv false). Canonical rails the
readouts up (100×, counted right), the variant collapses them (~1e-6). The
removed DSB depletion edges (36) moved which side the loop falls to; neither
is signal (compare specs/054).

### What this means

- The rule fixes the case it was built for (+30, verified in Neo4j). It
  loses for two reasons. (1) 38 held-out cases (CBL, EPS15, SH3GL1) were right
  only through edges that deplete an unchanged partner: CBL-family negative
  regulators had no other effect. (2) The registered test misreads Reactome's
  "Ub-X" complex products: 24 of the removed edges hit a real substrate,
  which costs the SMAD2/3 cases (−24). The second is a defect in this rule,
  not in the idea.
- So a faithful network loses the curators' "a ubiquitin ligase is a negative
  regulator" cases unless degradation of the ubiquitinated substrate is
  represented. The candidate, not registered: a ubiquitination depletes the
  reaction's own input complex that holds the changed substrate (the active
  receptor), not a free partner. The generator docstring records that
  targeting the complex was −0.5pp experimental on the May networks; it has
  not been measured since.
- A corrected test (fixes reason 2) would still lose reason 1 (−38) and the
  RNF168 coin, so it would very likely still fail the gate by itself; with a
  degradation rule for reason 1 it might not. Both are new registrations.
- Nothing is adopted. `LNG_UB_SUBSTRATE_CHANGED` stays registered, off.
