# specs/084 research: arm cct

Solver `1bb8803`, `--server DS_CARRIER_COMPONENT_TWIN=on`, canonical `20261020-2038_23ef3e1`,
against `results/8190d83/dflt`, results `results/1bb8803/cct`. The logs report 2,860 curator
and 2,496 experimental component edge-solves.

## Results

| split | net (fixed / broken) | p |
|---|---|---|
| curator held-out | **−8** (0 / 8) | 0.0078 |
| curator tuning | +6 (6 / 0) | 0.031 |
| experimental (RAF set aside) | **+4** (4 / 0) | 0.125 |

By pathway:
- **Cell Cycle Checkpoints +6 curator and +4 experimental:** ATM knockout and overexpression,
  nothing broken.
- **Signaling by PDGF −8:** PDGFRB knockout, at 8 readouts.

**The registered gate fails** on held-out −8 (below 0).

## Predictions

1. At least 6 of the 8 ATM cases fixed: **holds**, 6 curator plus 4 experimental.
2. Fewer than 2 pathways lose more than 5: **holds**, only PDGF.
3. No held-out prediction: −8.

## The PDGF loss (traced, `pdgf_trace.jl`)

The readout is PDGF:Phospho-PDGF receptor dimer (R-HSA-186811), with 15 variant copies.
Under a PDGFRB knockout, the copies holding the PDGFRB member read 0.25× under `off`. Their
supply is the specialization bridge from the generic receptor complex (specs/067). specs/068
counts that bridge as identity, not outside supply, so the generator lists the loop as
closed.

The rule added AND edges from the ligands' starting copies (PDGFA, PDGFB) to the copies
with no starting copy of their own entity. Three of them went from 0.25× to 1.0×, and the
readout mean went from 0.75× to 0.85× (DOWN → NORMAL). The receptors, which have no
starting copy, are not among the components read.

The rule fired on forms that already had a supply, the bridge. As registered it applies to
"a closed carrier free form with no starting copy of its own entity". It should apply only
where the pathway gives the form no supply at all.

## Not adopted as registered

**Candidate amendment 1 (post hoc, not run):** skip a closed free form that keeps any
activating in-edge after its loop edges are dropped, such as a specialization bridge. The
ATM complex has none, so it would still be fed.

## Amendment 1: arm cct2

Solver `1b52a97` (amendment 1), same build and baseline, results `results/1b52a97/cct2`.
Logs: 2,780 curator and 2,496 experimental component edge-solves (80 fewer curator than
`cct`: the bridged PDGF copies are skipped).

| split | net (fixed / broken) | p |
|---|---|---|
| curator held-out | **0** (0 / 0) | 1 |
| curator tuning | **+6** (6 / 0), all ATM in Cell Cycle Checkpoints | 0.031 |
| experimental (RAF set aside) | **+4** (4 / 0), all ATM | 0.125 |

**The gate passes:** held-out 0 ≥ 0, experimental +4, no pathway below −10. Post hoc.

Predictions:
1. ATM as in `cct`: **holds**, +6 / +4.
2. PDGF moves by at most 2: **holds**, 0.
3. Held-out ≥ 0: **holds**, 0.

Only Cell Cycle Checkpoints moves, and only ATM, on both axes. The other 14 fed entities
(Base Excision Repair, NOTCH, MHC class I and the rest) change no scored case.
