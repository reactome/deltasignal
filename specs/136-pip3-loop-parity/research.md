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
