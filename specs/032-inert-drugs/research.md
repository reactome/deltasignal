# specs/032 — Drug entities as inert participants

**Status:** pre-registered 2026-09-26, before any arm has run.

## Why

The experimental axis is where we trail MP-BioPath's published result:
- ours 67.14% (570/849);
- theirs 75.74%;
- paired against their supplementary table, 67.4% against 73.8%, −54 cases.

**One pathway carries 32 of the 54:** RAF/MAP kinase cascade (R-HSA-5673001),
where we score 28.6% and they score 93.9%. On MP-BioPath's own networks
DeltaSignal scores RAF at 100%, so the cause is our network, not the solver.

**Traced (canonical build `20260926-1221_590301c`, HRAS up):**
- The baseline is healthy: 1,571 of 1,571 nodes at 1.0, converged in 5
  iterations.
- Under HRAS up the solve does not converge (535 iterations), and 127 of the
  128 nodes of the second cyclic component read 0.
- That component's only negative edges are 4 regulator edges from
  **drug-bound complexes**: "activated RAF:scaffold:MAP2K:MAPK complex:dual
  (single) mechanism MAP2K inhibitors". They inhibit "RAF phosphorylates
  MAP2K dimer" and "MAP2Ks phosphorylate MAPKs".
- The drug is a root held at baseline, so raising RAS or RAF raises the
  drug-bound complexes, which shut the kinase steps off.
- Every RAS/RAF perturbation reads DOWN, including the NF1 knockdown, which
  should read UP.
- Reactome 97 curates these drug actions inside the pathway. MP-BioPath's 2019
  networks predate them.

**Scope in the catalog:**
- 202 drug entities, or entities that contain a drug: 134 ChemicalDrug, 29
  ProteinDrug, 39 complexes.
- They occur in 24 of 92 pathways: RAF 37, DNA Repair 28, Fanconi 24, Mitotic
  G1 24, MET 23, PDGF 19, KIT 18, IFN α/β 14, ERBB2 14, and others.
- No benchmark pin and no readout is a drug entity (checked on both axes).

## The rule

A benchmark case describes a cell **without the drug**. So an entity that is,
or contains, a Reactome `Drug` is held at baseline:
- it stays a participant (an AND input of fold 1.0; an inhibitor at baseline
  divides by 1);
- it can never carry or amplify a perturbation.

This is the cofactor rule (specs/007) applied to drugs. The curated edges stay
in the network, which keeps faithfulness to Neo4j. Pinning at 0 is rejected:
under `divide` inhibition, a regulator at 0 de-represses its target tenfold,
which models drug *withdrawal*, not drug absence. An explicit observation
still wins.

**Implementation:**
- **Generator:** ships `drugs.csv` beside `cofactors.csv`, listing the
  pathway's drug entities from the release's own data. It is an additive file;
  the network is byte-identical.
- **Solver:** `DS_DRUG_MODE` = `propagate` (default, byte-identical) or
  `inert`. The solve reports how many nodes were held, and the benchmark
  refuses an `inert` arm whose solves hold none.

## Pre-registration

**Arms:**
- `drugctrl`: code defaults.
- `druginert`: `DS_DRUG_MODE=inert`.
- Both on the same build, the first one that ships `drugs.csv`, through
  `scripts/run_arm.sh`.

**Adopt only if all three hold:**
- experimental net ≥ +10 with p < 0.05;
- curator held-out net ≥ −15 (the regeneration noise floor);
- no pathway loses more than 10 curator held-out cases.

The rule has no parameter, and nothing is tuned on either axis.

**Predictions, written before the arm:**
- RAF experimental rises by at least 20 of its 49 cases;
- RAF curator improves;
- the other drug pathways move by less than 10 each;
- pathways without drugs are byte-identical.

**Also reported:**
- curator tuning;
- per-pathway net on both axes;
- concentration: pathways moved, and distinct perturbations;
- the convergence count for RAF under each arm;
- the paired gap against MP-BioPath's supplementary table.

## Result — NOT ADOPTED (2026-09-26)

**Setup:**
- Build `20260926-1359_f2ba9cc_drugs`; both arms at solver 093282f, through
  `scripts/run_arm.sh`.
- `drugctrl`: code defaults, which report drug rule `propagate` with 0 held.
- `druginert`: `DS_DRUG_MODE=inert`; it held 6,940 drug node-solves on curator
  and 3,016 on experimental.
- Pins are identical between the arms (1,269 nodes over 864 perturbations), and
  so is convergence (1,659 of 1,725 on curator, 212 of 244 on experimental).

| split | net | fixed / broke | p | pathways |
|---|---|---|---|---|
| curator held-out | **+30** | 30 / 0 | 1.9e-09 | 4 (+4 / −0): MET +10, ROCKs +8, KIT +6, VEGF +6 |
| curator tuning | −16 | 5 / 21 | 0.0025 | RAF −17, ERBB2 +1 |
| **experimental** | **−5** | 3 / 8 | 0.23 | 1 (RAF) |

**Gates:**
- experimental must be at least +10 with p < 0.05. It is −5: **FAIL.**
- held-out must be no worse than −15. It is +30: pass.
- no pathway may lose more than 10 held-out cases: pass.

**The predictions failed.** RAF experimental did not rise by at least 20: it
fell. RAF curator did not improve: it is −17. **Not adopted.** `DS_DRUG_MODE`
stays, default `propagate`, as the record.

### Why the mechanism in this spec was incomplete

HRAS up in RAF (R-HSA-5673001), solved on the drugs build:
- `propagate`: 537 iterations, not converged, 287 nodes at 0.
- `inert` (40 drug nodes held): still not converged, 235 at 0.
- In the 128-node component, 127 nodes are at 0 under `propagate` and 124 under
  `inert`.

**The drug-bound complexes' inhibitory edges were real but not the cause.**
The component is a **positive loop closed through catalyst edges**:
- "Phosphorylation of RAF" (80 variant copies) takes the phosphorylated MAP2K
  dimers (p-S218,S222 MAP2K1 dimer; p-S222,S226 MAP2K2 dimer; the p-2S
  MAP2K1:MAP2K2 heterodimer) as **catalysts**.
- Inside this component those dimers are produced only by "Dissociation of
  RAS:RAF complex".
- That reaction consumes the activated RAF:scaffold:p-2S MAP2K:p-2T MAPK
  complex, which is made downstream of RAF phosphorylation.
- "MAP2Ks and MAPKs bind to the activated RAF complex" consumes F-actin, CNKSR2
  and Ca2+, which in this component come back only from the same dissociation.
- The external inputs are healthy: p21 RAS:GTP:RAF complex is at 100x, and ATP
  and the kinases are at 1x.
- Pushed hard, the loop does not settle, and it ends in the all-zero state.

This is the loop knife-edge (specs/013, 014) on a component built by
**catalyst and dissociation recycling**. Holding the drugs removes 3 of the
127 zeros.

**Open, not traced:**
- whether the catalyst edges (p-MAP2K → Phosphorylation of RAF) match what
  Reactome curates for that reaction, or were derived;
- why the held-out +30 (MET, ROCKs, KIT, VEGF) moved; the drug rule is sound
  there, but it was not gated separately.
