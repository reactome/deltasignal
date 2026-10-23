# specs/092: `no_path` on the network the solver runs

**Status:** analysis only, 2026-10-23. No rule is proposed or measured here.

Adam, 2026-10-23: "and then do I3 after you are done" (candidate I3 of the 2026-10-22 review,
finding 31: specs/079's `no_path` ceiling is read from the generator's network, not the solver's).

## Method (`solved_reach.jl` → `solved_reach.tsv`, `solved_reach_cut.tsv`)

Canonical `20261022-2011_233665c`, main's defaults (`results/8a94f2e/dflt`). For every wrong
valid case, a directed path from a pinned uuid to a readout uuid is looked for on two networks:
- **generator**: the logic network plus boundary edges, as the bench reads it
  (`build_adjacency`);
- **solved**: `solved_carrier_network` (recycle edges stripped, carrier-loop edges dropped, root
  twins, recycled supply, component twins and sequestration edges added), plus each pool path's
  steps as edges, because the pool update moves the protein along them (the last step lands in
  its target state with no logic-network edge).

The generator-network result reproduces the bench's labels exactly: 0 of 1,100 curator and 0 of
98 experimental `no_path` cases have a route there.

## Results

| | curator | experimental |
|---|---|---|
| wrong, labelled `no_path` | 1,100 | 98 |
| … with a route on the solved network | **2** | **0** |
| wrong, labelled routed | 1,494 | 123 |
| … with **no** route on the solved network | **80** | **17** |

1. **The specs/079 ceiling holds today.** The solver's added edges route 2 of 1,100 curator
   `no_path` errors (the review's correction stands as a matter of principle: specs/084's
   component edges did route 10 such cases, which are now right and so no longer counted).
2. **The solver cuts 80 curator and 17 experimental routes**, all by carrier-loop drops: 72
   by the loop drops alone, 8 with a depletion edge as well. They are labelled
   `propagator_missed` or `false_positive_change` but are really "no route in the network the
   solver runs".
3. Without the pool steps as edges, 116 more curator cases looked cut "by recycle edges". They
   are not: the pool update carries them. The bench's own reachability (review finding 30)
   counts the recycle edge, which gives the same answer for the wrong reason.

## Who is cut (curator, by gene)

| pathway | gene | cases | cut by (`which_edge.jl`) | right before specs/077 (`bf4b688/rst`) → after (`24e3c1f/ocs`) |
|---|---|---|---|---|
| Mitotic G1 | RB1 | 30 | open-loop drops (specs/077), several free forms together | 16/50 → 20/50 |
| MET | CBL | 10 | open-loop drops into HGF:p-MET dimer (R-HSA-6806977) | **20/20 → 10/20** |
| TGF-beta | SMAD7 | 9 | closed-loop drops into ZFYVE9 (SARA, R-HSA-206180) | 5/14 → 5/14 |
| FGFR1–4 | PTPN11, SPRY2 | 16 | loop and depletion drops together | unchanged |
| others | CCNE1, MYC, MAX, SRC | 12 | | |

Experimental: RB1 knockout 4, MYC knockout 3, CCNE1 overexpression 2 (Mitotic G1); ERBB2 3;
HDR 2.

## Reading

- The "no route" half of the errors (specs/079) is still out of reach of solver rules as they
  stand: a solver rule now routes 2 of 1,100.
- The real solver-side leak runs the other way: 97 wrong cases whose route the carrier-loop rules
  cut. **CBL in MET is the one class that was right before the cut** (specs/077 cost it 10, a known
  part of that adoption's 19 losses). The release that matters is the HGF:p-MET dimer freed by
  CBL's monoubiquitination step, an open loop under `fed_single`. RB1, SMAD7 and PTPN11/SPRY2
  were mostly wrong before the cut as well, so restoring the route alone would not fix them.
- A candidate rule for CBL/MET would keep an open release whose freed form is a receptor
  dimer the ligase marks rather than degrades; it is not designed here.
