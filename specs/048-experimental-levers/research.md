# specs/048 research: arm results

All arms ran on solver commit `4cb152d` with `scripts/run_arm.sh`, protocol
`DS_PIN_SCOPE=root_cycle`. Scoring is `bench/analysis/arm_compare.py`, paired on
(pathway, gene, direction, readout), because a rebuild re-mints every uuid.

| arm | build | server override | results |
|---|---|---|---|
| canonical | `20261009-2307_132d4b8` | — | `results/0f0a11e` |
| `pathctl` | `20261010-0259_a43e20d_path` | `DS_PATHOGEN_MODE=propagate` | `results/4cb152d/pathctl` |
| `path` | `20261010-0259_a43e20d_path` | `DS_PATHOGEN_MODE=inert` | `results/4cb152d/path` |
| `bind` | `20261010-0341_a43e20d_bind` (`LNG_BIND_STOICH=1`) | — | `results/4cb152d/bind` |

Valid cases are identical in every arm: curator 23,511 (18,481 held-out, 5,030
tuning), experimental 845. Convergence is identical in every arm: 1,654 of
1,725 curator solves, 215 of 244 experimental solves.

## Control: `pathctl` reproduces canonical

- Every one of 24,100 curator and 849 experimental cases has the same
  prediction, validity and category as canonical.
- 9 curator values differ by more than 1e-6. Eight are DSB Repair
  overexpressions reading 75–89× that move by at most 0.0016 (relative 2e-5).
  The ninth is one EGFR value at 0.2965 that moves by 2e-6. None is near a
  threshold.
- These shifts are not traced. The build's networks have the same node and
  edge counts as canonical. Predictions match, which is what the amendment
  required, so `path` is read against `pathctl`.

## Lever A, arm `path`: NOT ADOPTED (fails the pathway floor by one case)

| axis | split | `pathctl` → `path` | net (fixed/broke) | McNemar p |
|---|---|---|---|---|
| experimental | all (tuning) | 69.11% → 69.94% | **+7** (7/0) | 0.016 |
| curator | held-out | 88.99% → 88.93% | **−11** (1/12) | 0.0034 |
| curator | tuning | 76.14% → 76.42% | +14 (14/0) | 0.00012 |
| curator | all | 86.24% → 86.26% | +3 (15/12) | 0.7 |

- The servers held 128 node-solves (experimental) and 14,642 node-solves
  (curator).
- Two pathways moved, on 10 perturbations:
  - PIP3, all PDPK1 knockout: +7 experimental, +14 curator.
  - DDX58/IFIH1: −11 curator held-out, across 9 perturbations.

**Predictions:**
- **A-P1 met.** All 7 PDPK1-knockout DOWN readouts go from 100× to 0×.
- **A-P2 met.** The 3 PDPK1-overexpression UP readouts stay UP; they read
  100× where before they read 64.5×.
- **A-P3 met.** Experimental net is +7, all in PIP3, and no other
  experimental case changes.
- **A-P4 met.** Held-out net is −11, inside ±15.

**Decision rule:**
- Experimental net > 0: yes.
- Held-out ≥ −15: yes.
- No pathway loses more than 10: **no.** DDX58/IFIH1 loses 11 (1 fixed, 12
  broken).

By the pre-registered rule, lever A is not adopted. `DS_PATHOGEN_MODE` stays
in the solver at `propagate` (the default), as a measured record, like
`DS_DRUG_MODE`.

### Traced: the DDX58/IFIH1 loss

- **What is held.** `R-HSA-168928/pathogens.csv` has 633 rows over 45 base
  entities. 33 of the 45 are human complexes that contain a viral component.
  Most of them are the RIG-I/MDA5 signalling platform:
  - `viral dsRNA:IFIH1, viral dsRNA:K63polyUb-DDX58:MAVS:...` with TRAF3,
    TBK1/IKBKE, IRF3/IRF7, RIPK1:FADD:CASP8/10, and IKK;
  - R-HSA-918199, 933470, 933479, 9705147, 8983916 and 933478 among them,
    with up to 48 variant rows each.
- **Why it breaks.** Every complex on the route from a host signalling
  protein to its readout contains the viral dsRNA ligand. Holding those
  complexes cuts every host perturbation that runs through the signalosome.
  - CYLD, IKBKB, MAP3K1 and CASP8, knocked out or overexpressed, all read
    exactly 1.0 at readouts 909690, 177673, 877351 and 933478, where they
    were correct before.
  - The one fix is a CREBBP knockout at 877351. It goes from 3.5× to 0.5×,
    because the `CREBBP:NS1` antagonist complex is now held.
- **The two curations differ.**
  - In PIP3, the pathogen is an antagonist: SARS-CoV-2 N:M sequesters a host
    kinase. A cell without the virus has no such complex, and holding it is
    faithful.
  - In DDX58/IFIH1, the viral RNA is the stimulus the pathway exists to
    sense, and the curator ground truth assumes it is present. A-P4 named
    this risk and did not predict its direction.
- **A narrower rule would separate the two.** It would hold a complex only
  when its pathogen part is a protein bound to a host protein, not a
  nucleic-acid ligand. That rule is post hoc. Its experimental upside is
  capped at the same +7, and it would need its own pre-registration and arm.
  It is not run here.

## Arm `bind` (B3, `LNG_BIND_STOICH=1`): ADOPT on faithfulness

| axis | split | canonical → `bind` | net |
|---|---|---|---|
| curator | held-out | 88.99% → 88.99% | 0 (0/0) |
| curator | tuning | 76.14% → 76.14% | 0 (0/0) |
| experimental | all | 69.11% → 69.11% | 0 (0/0) |

**Structural checks:**
1. **Met.** The build is complete: 213,171 nodes and 516,528 edges, against
   213,061 and 516,244 for canonical.
2. **Met.** Coverage is identical on both axes.
3. **Met.** In every copy of R-HSA-5675376, p-MAPK1 and p-MAPK3 go to their
   monomers and the phosphorylated dimers to their dimers.
4. **Counted before scoring.** Binding ties fall from 64 to 58 catalog-wide.
   - Edges change in 10 pathways: +84, +16, +40, +24, +72, +6, +12, +4 (RAF),
     +16 and +10.
   - Stable-id sets are identical in all of them.

**Predictions:**
- **bind-P1** registered no direction, and nothing moved.
- **bind-P2 missed.** RAF values move, but no prediction does:
  - 28 experimental values change, all in RAF.
  - The largest change is BRAF knockout at readout 1268261: 70.38× → 70.22×.
  - The 12 knockout→UP cases still read 36–70× through
    `R-HSA-169291::pool`.
  - Correcting the binding closes the curated MAPK ⇄ p-MAPK cycle, but the
    knockout→UP class does not depend on it. Two pools share 246 nodes, and
    that is still unsolved.
- 226 curator values move, none across a threshold. Most are IFN-γ values
  near zero that change in the sixth decimal.

**Decision rule** (held-out ≥ −15 and experimental ≥ −15): met at 0 and 0.
B3 is adopted as a faithfulness fix. That means `LNG_BIND_STOICH` defaults
on in the generator (an LNG PR, which Adam merges). Canonical numbers do not
change.

## What this leaves

- The RAF knockout→UP class (12 of 21 RAF failures) needs pools that share
  intermediates. That is the lever B design question in Amendment 1. B1 and
  B2 are on `wip/048-mapk-pools`.
- PDPK1 (+7 experimental) is recoverable only by a pathogen rule that spares
  sensor ligands. That rule would need a new registration.
