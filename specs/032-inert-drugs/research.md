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
