# specs/081 diagnostic probe (2026-10-20, not gated, before any registration)

Mitotic G1 (R-HSA-453279), canonical build, current defaults, all 850 valid curator
cases. `probe.jl` adds a plain OR edge from a protein copy a reaction produces to the
starting (boundary-leaf) copy of the same entity, which feeds a root complex. Run it
in the dev container with the case list (gene, direction, readout, expected, predicted,
pins, readout uuids) as `ARGS[1]`. The unperturbed network stays at baseline in both
versions, and the baseline network reproduces `results/3e96e7c/dflt` on every case.

| version | edges added | fixed | broken | net |
|---|---|---|---|---|
| every produced copy (`DEADEND=0`) | 48 (7 entities) | 117 | 115 | +2 |
| only produced copies that lead nowhere (`DEADEND=1`) | 2 (RBL2, RBL1) | 95 | 69 | +26 |

- `DEADEND=0` breaks through loops on the cyclin:CDK copies. MAX knockout and
  overexpression read the same, and so do CCNE1 knockout and overexpression.
- `DEADEND=1` fixes the CDKN1B, CDKN2A and PPP2R1A cases, but gives the **wrong sign** for
  CCND1 and CDK4 overexpression. CCNA1 goes from 17.5× (UP, right) to 0.065× (DOWN, wrong):
  more kinase means more p-p130, then more free p130, then more of the repressing complex.
  It also breaks MYC/MAX knockouts and CCNE1 cases. A plain edge reads the cycle's flux as
  an amount. Only a conserved reading gets the sign right: the kinase moves p130 *out of*
  the complex, and total p130 is set by its gene.

Conclusion: the return step needs its own edge type, handled differently by the solver
and ignored by root detection (Adam, 2026-10-20). A plain supply edge is not enough.

## Structural count (`count_recycle.py`, 2026-10-20, before any registration)

Candidates: a produced protein copy that leads nowhere, a starting copy of the same
entity feeding a root complex, and a reaction route from that complex back to the
produced copy (at most 4 reactions, through forms carrying the protein, by reference
entity from Neo4j).

| filter | canonical | faithful (`dsub`) |
|---|---|---|
| any such route, proteins only, ≥ 2 reactions | 20 entity pairs, 10 pathways | — |
| + the protein is present only in a **modified** form somewhere on the route (a modification cycle, not a partner released unchanged) | **2: p130 and p107, both Mitotic G1** | 4: the same two, plus R-HSA-418996 in R-HSA-446728 |

The other 18 are partners released unchanged after a multi-step event: SMURF1/2 and
NEDD4L, BRAF from p-SPRY2, CK1α and AMER1 from the destruction complex, TLS polymerases,
and nuclear envelope proteins. That is the carrier-release class (specs/077), not a
modification cycle.

Mitotic G1 is a tuning pathway (one of the paper's ten), so on the canonical build this
rule cannot move the held-out split. It can move tuning and experimental.
