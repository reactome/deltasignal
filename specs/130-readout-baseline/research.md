# specs/130 research

CPPA training portion only; scored with `bench/analysis/score_130.py` at this branch (solver-free: frozen predictions
`20261027-1346_8b6527e_top114/results/8276c28/cppa_frozen`). Under the CPPA terms every measured number is kept in the
private notes (`~/deltasignal-catalogs/analysis/121/research_private.md`, "specs/130"); this file gives the verdicts.

One fix after registration, before any result was read: `control_baselines` created an empty entry for a non-numeric
value and divided by zero; values are now parsed before an entry is made. The rule is unchanged.

| prediction | verdict |
|---|---|
| **P1.** the identity reproduces specs/121's frozen training score | **holds** (equal to every printed digit) |
| **P2.** ≥ 90% of training samples have controls for their (set, cell line) | **holds** (over 99%) |
| **P3.** CV gain positive and < 0.03 | **holds** as stated, but its 95% interval includes 0 |
| **P4.** secondary balanced accuracy rises | **fails**: it falls under the chosen θ and under every θ |

**Verdict: the readout baseline does not earn its place** (the registered gate needs the CV gain's interval above 0).
It is not carried to the CPPA test scoring. The lowest background quantile (θ = 0.05) was chosen in four folds of five:
the transform helps the correlation most where it barely compresses, and the compression it applies costs correct
change calls (the balanced accuracy falls), the same trade specs/128–129 found for attenuation along edges.

## Reading

Two different sources of cell context have now failed the same CPPA gate: expression gating of every node (specs/121
phase 2) and the readout's measured resting level (here). With specs/129, every way tried so far of making the model
call fewer or smaller changes trades correct change calls for correct no-change calls at about even rates. What the
model lacks is not a damping of its calls but information about **which** readouts respond, and neither Reactome's
wiring nor resting levels supply it on their own.
