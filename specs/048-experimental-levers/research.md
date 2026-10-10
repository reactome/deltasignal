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
- **Case by case**
  (`~/deltasignal-catalogs/analysis/048/ddx58_trace.py`, routes in the
  `_path` network). All 12 broken cases fail in one of two ways:
  - **The readout itself is held (4 cases).** IKBKB and CASP8, up and down,
    at readout 933478. That readout is
    `dsRNA:RIG-I/MDA5:...:IKK complex`, which contains the viral ligand, so it
    is pinned at 1.0 whatever the perturbation does.
  - **Every route passes a held node (8 cases).** CYLD up and down at
    909690, 177673 and 877351, and MAP3K1 up and down at 177673.
    - Shortest routes before the pin: 19–31 nodes for CYLD, 9 for MAP3K1.
    - Routes that avoid held nodes: none.
    - On the CYLD routes, the first held node is `2x DDX58 ligand:2x
      DDX58:2xATP` (R-HSA-168906), RIG-I bound to viral RNA. On the MAP3K1
      route it is `dsRNA:DDX58/IFIH1:MAVS:TRAF2/TRAF6:MEKK1` (R-HSA-933482).
  - **The one fix.** CREBBP knockout at 877351 has a 7-node route with no
    held node. It moves from 3.5× to 0.5× because the `CREBBP:NS1`
    antagonist is held rather than knocked out with CREBBP.
- **Neo4j confirms the curation** (`ddx58_neo4j.py`, release 97):
  - "DDX58 ligand" (R-NUL-9013905) holds only viral RNAs: Rotavirus RNA, HCV
    5'-ppp poly-U/UC RNA, and the Influenza A dsRNA intermediate.
  - "IFIH1 ligand" (R-NUL-9038422) holds Measles, RSV A and Influenza A
    dsRNA.
  - Neither set has a host member. Reactome curates no form of the RIG-I/MDA5
    platform without the virus, so this is not a generator artefact.
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
B3 is adopted as a faithfulness fix. `LNG_BIND_STOICH` now defaults on in
generator PR #110 (commit 0123fbb), which Adam merges. Canonical numbers do not
change.

## What this leaves

- The RAF knockout→UP class (12 of 21 RAF failures) needs pools that share
  intermediates. That is the lever B design question in Amendment 1. B1 and
  B2 are on `wip/048-mapk-pools`.
- PDPK1 (+7 experimental) is recoverable only by a pathogen rule that spares
  sensor ligands. That rule would need a new registration.

## Amendment 2, arm `pathp` (pathogen proteins only): ADOPT

- **Build:** `20261010-1155_f58ea98_pathp`, with `LNG_PATHOGEN_PROTEIN=1`
  and `LNG_BIND_STOICH=0`.
  - It is complete, with 213,061 nodes and 516,244 edges, the same as
    `_path`. Node and edge counts match in every pathway, and `drugs.csv` is
    identical.
  - Only DDX58/IFIH1's `pathogens.csv` changes: 633 rows → 13, none of them
    new. What remains is SARS-CoV-2 M, N, 9b and nsp13 and their host
    complexes, and RSV NS1 with IRF3:NS1 and CREBBP:NS1.
- **Solver:** commit `e2daa7e`; `src/` is unchanged since `4cb152d`.
- **Results:** `results/e2daa7e/{pathpctl,pathp}`.
- **Held:** 1,754 curator node-solves (14,642 under `path`) and 128
  experimental, the same as `path`.

**Control.** `pathpctl` reproduces `pathctl` on every prediction, validity
and category. 10 curator values shift by more than 1e-6:
- the same DSB Repair readout, 5693527, by up to 0.026 (RNF168, 71.34 →
  71.32);
- the same EGFR value at 0.2965.

These are the shifts that canonical → `pathctl` showed. Each rebuild of
identical networks moves them, so something in that DSB component still
depends on labels, contrary to specs/013. It is small and far from any
threshold, but not traced; it is an open item.

| axis | split | `pathpctl` → `pathp` | net (fixed/broke) | McNemar p |
|---|---|---|---|---|
| experimental | all (tuning) | 69.11% → 69.94% | **+7** (7/0) | 0.016 |
| curator | held-out | 88.99% → 89.00% | **+1** (1/0) | 1 |
| curator | tuning | 76.14% → 76.42% | +14 (14/0) | 0.00012 |
| curator | all | 86.24% → 86.31% | +15 (15/0) | 6.1e-05 |

Two pathways moved, on two perturbations: PDPK1 knockout in PIP3, and
CREBBP knockout in DDX58/IFIH1. Nothing breaks on either axis.

**Predictions:**
- **P2-1 met.** Experimental net is +7 (7/0), the same 7 PDPK1-knockout
  DOWN readouts as `path`.
- **P2-2 met.** None of the 12 DDX58/IFIH1 cases that `path` broke changes.
  The CREBBP fix at 877351 stays (3.5× → 0.5×), and DDX58/IFIH1 net is +1.
- **P2-3 met.** Held-out net is +1, tuning +14, and no other pathway moves.

**Decision rule:**
- experimental > 0: yes;
- held-out ≥ −15: yes;
- no pathway loses more than 10: yes, nothing loses.

**Adopted.** This means:
- `LNG_PATHOGEN_PROTEIN=1` and `DS_PATHOGEN_MODE=inert` become the defaults;
- `bind` (B3) is adopted too, so by the decision rule the combination is
  measured once before it becomes canonical.

**Caveat.** The rule was written after `path`'s trace, so its predictions
confirm the trace rather than test the idea blind. The gain rests on one
curated SARS-CoV-2 complex in PIP3, and the experimental +7 is one
perturbation.
