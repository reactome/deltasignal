# specs/100 research: arms rnd (control) and rnt

Variant build `20261023-2341_135af3c_rna2` (canonical content plus `rnas.csv`), solver
`61cdb8a` (main with specs/088's default). `rnd` (defaults) reproduces canonical's defaults
(`results/90f7df4/dflt088`) on every case and value, both axes. `rnt` adds
`--server DS_RNA_TWIN=on` (172 curator and 68 experimental RNA supply edge-solves). Read with
`holdout_report.py --set-aside RAF_MAP_kinase_cascade --pair-on-valid`.

| split | net (fixed / broken) | p |
|---|---|---|
| curator held-out | 0 (0 / 0) | 1 |
| curator tuning | **+28 (28 / 0)** | < 1e-4 |
| experimental (RAF set aside) | **+4 (4 / 0)** | 0.125 |

Every changed case is TP53 in PIP3 activates AKT signaling: knockout 0× → UP (14 curator, 2
experimental), overexpression 83.5–100× → DOWN (14 curator, 2 experimental). Converged solves
unchanged (1,724 / 1,725, 244 / 244). **The gate passes; nothing is broken.**

## Predictions

1. TP53 in PIP3, at least 20 of 28 curator and 3 of 4 experimental fixed: **holds, 28 / 28 and 4 / 4**.
2. Only the 4 pathways with an RNA pair can move: **holds** (only PIP3 moves a class).
3. PTEN overexpression stays right: **holds** (no PTEN case changes).
4. Held-out: 0.

## Reading

Recommended for adoption. The rule follows from what an mRNA is (one pool, however many node
copies the generator made of it), touches only RNA species (135 in 19 pathways), and fixes a whole
inverted class with nothing broken. Experimental would go from 635 to 639 of 849.
