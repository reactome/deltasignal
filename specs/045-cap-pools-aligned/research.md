# specs/045: over-cap reactions as set pools, aligned with their expanded neighbours

## Why (2026-10-04)

**Adam's point.** Expanding sets into per-member copies is the same problem
solved once per member, so it should agree with a set-node network such as
MP-BioPath's (RAF experimental: 100% with our solver, against 29% for ours).
Where it does not, there is a bug.

**The bug located.** Expansion is exact only while every copy uses one member
per set.
- "MAP2Ks and MAPKs bind to the activated RAF complex" (R-HSA-5672972) needs
  **12,000** combinations: 80 × 25 × 2 × 3.
- Downstream steps multiply that again. Uncapped (`LNG_MAX_VARIANTS=10^6`),
  RAF alone ran a 30 GB machine out of memory.
- Past the cap (512) the generator bundles every alternative into ONE copy
  that requires all of them: 28 AND inputs. That is the opposite of
  expansion, and it is where specs/034 §8 traced RAF's collapse.

The tractable form of "solve once per member and recombine" is one OR node per
set. `LNG_CAP_POOLS=1` (specs/036) does that for over-cap reactions only, which
turns the step's 28 inputs into 4. It has only ever been measured combined
with `DS_DEDUP_ACTIVATORS`, which is itself damaging (specs/012 addendum
2026-10-04: held-out −56). So cap pools alone are unmeasured.

## Alignment (Adam: "make sure it aligns with the section split apart by sets")

RAF regenerated with `LNG_CAP_POOLS=1` at generator 1491276. The check is
`~/deltasignal-catalogs/analysis/045/alignment_check.py`:
- one pooled reaction;
- 13 member slots;
- 1 output, not a dead end.

Two member slots (MAPK1, MAPK3) read copies with no producer while other copies
are produced. **Traced; not misalignment:**
- the "producers" are (a) "Cytosolic DUSPs dephosphorylate MAPKs", whose
  curated output is a different set, "MAPK monomers and dimers"
  (R-HSA-5675361), dissolved into members by finding F10; and (b) the MAPK1/3
  dimer readout sinks of the boundary layer;
- nothing in RAF produces the consumed set, MAPKs (R-HSA-169291). It enters
  from outside;
- the diagram draws it as one glyph (2203), the input of all three consumer
  reactions.

Curation note: the MAPK cycle is not closed in Reactome. DUSPs release a set
different from the one the binding step reads.

The check is to be refined to glyph identity (specs/043's criterion) and run
catalog-wide on the cap-pools build.

## Pre-registration (before the arm)

- **Build:** `--variant cappools`, generator main 1491276, `--env LNG_CAP_POOLS=1`.
- **Control:** `f7ctrl` (20261003-1058_1491276_f7ctrl), the same generator commit.
- **Scoring:** solver main (`DS_SCC_ORDER=flow`), code defaults.

**Predictions:**
1. RAF: the binding step loses its all-required bundle.
   - Experimental RAF moves up. It is not expected to reach MP-BioPath's
     100%, because specs/036 traced further amplifiers (drug-bound
     inhibitors, input+catalyst squaring, one untraced).
   - Named: RAF experimental net > 0.
2. Held-out: within ±15. Over-cap reactions are a minority of the catalog.
3. Experimental: ≥ 0 overall.
4. Report both axes, held-out and tuning, the pathways and genes moved,
   McNemar p, convergence, and the catalog-wide alignment check.

**Decision rule:** adopt `LNG_CAP_POOLS=1` as the generator default if held-out
is ≥ −15, experimental ≥ 0, and no misalignment in the refined check is left
unexplained. If RAF does not move, trace it against MP-BioPath's network to the
first divergent node.

## Result (2026-10-04, solver 5221abf with `DS_SCC_ORDER=flow`, `cappools` against `f7ctrl_flow`)

| | Net | Fixed / broken | p |
|---|---|---|---|
| Curator held-out | **0** | 0 / 0 | – |
| Curator tuning | +21 | 21 / 0 | RAF only |
| Experimental | **+14** | 14 / 0 | 0.00012 (RAF only, 6 perturbations) |

- Convergence is identical (1,668 / 1,725 and 222 / 244).
- ARAF, BRAF and RAF1 are now right both as knockouts (DOWN) and as
  overexpressions (UP).
- Experimental is 68.90% → about 70.6% of all 849 cases. That is more than
  the 8 cases the label-free order cost.
- Every gain is in RAF, so this fixes RAF's bundle; it is not catalog-wide
  evidence. The other over-cap reactions moved no prediction.
- The glyph alignment check (`glyph_alignment.py`) finds the same result on
  `f7ctrl` and on `cappools`:
  - 77 produced/consumed copy pairs drawn as different glyphs (faithful);
  - 1 same-glyph pair in Class I MHC (R-ALL-29926), present before;
  - 9 pathways have no diagram.

  The check reads only each pathway's own (or nearest ancestor's) diagram, so
  it is a lower bound.

**Decision: adopt** (held-out ≥ −15, experimental ≥ 0, no new misalignment).
Generator PR `feat/cap-pools-default`: the default is on, and code-review F6
(the dependence of `CAPPED_IDS` on the cache) is fixed. The default build
equals the measured build, and a cached rebuild equals a fresh one (RAF,
Class I MHC).

**RAF is not finished.** MP-BioPath's network scores 100% on RAF experimental
with our solver; ours now scores higher than before, but not that. The next
step is a node-level trace of RAF against MP-BioPath's network, to the first
divergence.
