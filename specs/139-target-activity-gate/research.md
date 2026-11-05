# specs/139 research

CPPA training portion only; scored with `bench/analysis/139/score_139.py` at `07b90ea` (solver-free), markers from
`bench/analysis/139/markers.py`. Predictions: specs/134 amendment 2 (`top114 e4809b8/cppa134_abundance`) and, beside
it, the specs/121 frozen run (`8276c28/cppa_frozen`). Under the CPPA terms every measured number and every per-class
finding is kept in the private notes (`~/deltasignal-catalogs/analysis/121/research_private.md`, "specs/139"); this
file gives the verdicts.

| prediction | verdict |
|---|---|
| **P1.** identity, markers left out, within 0.01 of specs/121's figure | **fails**: leaving the drugs' own markers out lowers it by more than 0.01. With the markers in, the identity reproduces specs/130's figure exactly (frozen run), so the scorer is right |
| **P2.** at least half of marker-drug samples eligible | **holds** |
| **P3.** more measured changes where the target is active at rest | **holds**, weakly |
| **P4.** the gate fails | **holds** |

**Verdict: the gate does not earn its place.** Neither prediction set passes the registered judge: cross-validation
chooses no gate (or, in one fold of the frozen run, a gate whose gain does not beat the random-gate null), and the
general phospho-level comparator chooses no gate anywhere. Not carried to the CPPA test scoring.

## Reading

Traced in the private notes: the model's change calls are right about equally often in cells where the target is
active at rest and in cells where it is not, so the gate removes right and wrong calls alike. Three sources of cell
context have now failed the CPPA gate: expression of every node (specs/121 phase 2), the readout's resting level
(specs/130), and the target's resting activity (here). The model's CPPA errors are not concentrated in particular
cells; a cell-level switch does not reach them.
