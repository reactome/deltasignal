# specs/136 research

## R1. Control and arm P (2026-11-05)

Code `d370c52`, frozen build `20261026-1349_8b6527e`, `run_arm.sh --axis curator --axis experimental --axis korkut`.

**Control `k136_ctl`** reproduces its references case for case: korkut `k134_abund_a2` (every class and value), curator
and experimental `ce631f0/a124_c100` (the frozen run; 24,100 and 849 cases, 0 class differences). Convergence: curator
1,724 of 1,725 solves, experimental 244 of 244, korkut 1,493 of 1,493.

**Arm P `k136_P`** (`DS_METABOLITE_POOLS=on`): no class change on any axis. Values move only in PIP3 (curator 393
cases, experimental 180, korkut 145). SRC knockout, AKT pS473 in PIP3: 0.212 → 0.334, still DOWN.

**Parity census, P.** With the conversion reactions R-HSA-2316434 (PI3K) and R-HSA-199456 (PTEN) taken out of the
iterated graph, as the pool solve does, PIP3's strongly connected component falls from 425 nodes to 51 and the cycles
through PIP3 (≤ 16 edges) from 5,997 to 45: the 39 even MKRN1 cycles and the 6 odd THEM4/TRIB3 ones. The 372 gain-1
conversion cycles and the 5,580 cycles in which PTEN raised PIP3 are gone.

Q1 holds. Q2 holds (the switch stays railed; no call moves). Q5 for P holds (curator 0).

## R2. Arms L and PL, and the gate (2026-11-05)

Same code, build and axes. Convergence identical to the control in every arm.

| arm | curator held-out | curator tuning | experimental | korkut class changes | SRC KO AKT pS473 (PIP3) |
|---|---|---|---|---|---|
| control | | | | | 0.212 |
| P | 0 | 0 | 0 | 0 | 0.334 |
| L | 0 | 0 | 0 | 0 | 0.306 |
| PL | 0 | −28 (0 / 28, p < 0.0001) | −5 (3 / 8, p 0.23) | 0 | 0.798 |

L flags 1,256 slots across the 92 pathways that specs/083 did not (868 of them in Interferon α/β, 6 in PIP3;
`count136.jl`), and moves values in 351 curator cases (PIP3 210, Interferon α/β 100), none across a cutoff.

**PL's losses are all KIT and VAV1 knockouts in PIP3.** Both read AKT readouts DOWN in the control (0.171, 0.005) and
NORMAL under PL (0.91–0.98). KIT and VAV1 are each in one of about 30 PI3K activator complexes, like SRC, so they were
DOWN through the same switch. Curators expected DOWN for every one; the experiments agree for KIT (7 cases) and say
NORMAL for 3 of 4 VAV1 cases (the 3 fixed). With the switch damped, one activator out of thirty no longer moves PIP3.
Whether curators' DOWN is the biology (a cell driven by that receptor) or the switch's artefact is not decided here;
specs/108's criterion would need the 0.5× / 2× rail test on these cases.

**Gate: fails for all three arms.** korkut has no class change under any arm, so the balanced-accuracy difference is 0
and its interval is not above 0. The curator held-out guard passes for all three (no case moves). No arm goes on to
CPPA.

Q2 holds. Q3: L alone does not unrail the switch (holds); "fewer than 1,000 slots" is refuted (1,256). Q4 holds: PL
moves more than either and the gate fails. Q5: held-out 0 for every arm (holds); PL's −28 is tuning.

## What this says

Adam's parity question separated two things the PIP3 component mixed. The 5,952 conversion cycles (372 even, 5,580
odd) were the PIP2 ⇄ PIP3 interconversion read as signal, and in 5,580 of them PTEN raised PIP3; the lipid pool removes
them, and the remaining loops through PIP3 have the sign the biology has: the MKRN1 loop is positive feedback. That
loop is curated, and it is a switch. Taking PIP3's repeated counting out (PL) moves the switch most of the way off its
low state (SRC knockout 0.21 → 0.80), but no korkut call changes, because 0.80 is still under the 0.85 cutoff. The
curator and experimental calls it does change are knockouts of single PI3K activators that were DOWN only through the
switch. Nothing adopted. The PIP3 overcalls on korkut are the switch plus the 15% call cutoff, the same cutoff question
specs/135 left to Adam.
