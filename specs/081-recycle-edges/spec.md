# specs/081: a modified protein returns to the complex it starts in (p130, p107)

**Status:** pre-registered 2026-10-20, before any build or arm; amendment 1 (the `recycle` edge) before the catalog build and the gated arm. **Gate passes** (arm `rcp`): held-out 0, tuning +70 (70 / 0), experimental +8 (10 / 2, p 0.039), only Mitotic G1 moves. Results in `research.md`. **Adopted** (Adam, 2026-10-20).
structural count are in `probe.md`.

Adam, 2026-10-20: "we need to make things that are biologically correct. I say A",
and, on how: "maybe this needs a special type of edge in the logic networks so that
deltasignal can handle them differently. for example we would still know what the
root node is."

## The defect (specs/080)

In Mitotic G1:
- p130 (RBL2) starts in RBL2:E2F4,E2F5:TFDP1,TFDP2 (R-HSA-1226089), which represses E2F
  targets through the DREAM complex.
- Cyclin D:CDK4/6 phosphorylates it, and p-p130 leaves the complex (R-HSA-1226094).
- PP2A dephosphorylates p-p130 back to free p130 (R-HSA-1363276).
- Free p130 rebinding E2F4 is the step Reactome does not curate as a reaction, so the
  complex is a root. The generator builds it from its own starting copy of p130, and the
  free p130 that PP2A makes leads nowhere.

The cycle is never closed. Every CDK4, cyclin D, CDK-inhibitor and PP2A perturbation stops
at free p130, and 187 Mitotic G1 cases read NORMAL (181 curator, 6 experimental, all
expected to change). p107 (RBL1) is the same, with R-HSA-1226088, R-HSA-1226095 and
R-HSA-1363274.

MP-BioPath's network handled it with hand edits. It cut the cycle and added regulation
onto DREAM formation: CDK4:cyclin D, SCF and Ub inhibit it, PP2A activates it. On those
networks our solver scores 146 of these 187 cases.

## The rule

**The special link is a pool row, not a network edge.** The generator writes the cycle as
a conserved pool. The return step is a pool transition whose reaction is the
dephosphorylation and whose target is the starting complex. The logic network is
unchanged, so the starting copy of p130 stays a root and is still what a p130
knockout pins.

**Generator** (`LNG_RECYCLE_POOLS`, default off; LNG branch from `5e663d9`). A recycle
pool is found when all of these hold:
- **P:** a copy of a single protein that a reaction produces, leads nowhere (no
  out-edge) and is not a dissociation sink.
- **B:** a starting copy of the same entity that feeds a complex C by an `assembly` edge.
- **Route:** C reaches P through at most 4 reactions (input to output), through forms that
  carry P's protein.
- **Modified on the way:** some form on the route carries the protein with a different
  modification signature from P (the pool finder's own signature). A partner released
  unchanged is excluded; that is the carrier class of specs/077.

The pool:
- **States:** every copy of C that B feeds, all resting (base), plus M, the first
  modified form on the route.
- **Intermediates:** route forms between M and P.
- **Paths:**
  - forward: C copy → … → M, through that copy's reaction copies;
  - return: M → … → C copy, whose last step is the reaction that makes P.
- A candidate any of whose nodes an existing pool already claims is skipped and counted.
- On the canonical build the candidates are the p130 and p107 cycles in Mitotic G1
  (`count_recycle.py`, 2 of 20 route candidates).

**DeltaSignal:** a pool may have several base states. The baseline share is set by
distance from the nearest base, so each resting copy of C starts at full weight. One
base, as in every existing pool, is unchanged. The pool is solved by cycle balance as
before: the total p130 is conserved, and the kinase and phosphatase set its split. No
other solver change. Tests: `test_cycles.jl`, a two-base pool whose kinase drive
up lowers both bases.

**Benchmark:** no change. Roots are read from the logic network, which this does not alter.

## Arms

| arm | build | solver | against | gated |
|---|---|---|---|---|
| `rcp` | variant `--env LNG_RECYCLE_POOLS=1` | the commit that adds multi-base pools | `results/3e96e7c/dflt` (canonical) | yes |

The variant build must equal canonical everywhere except the Mitotic G1 pool tables.
That is checked and committed before the arm.

**Gate (specs/077's):** curator held-out net ≥ 0, experimental net ≥ −5, and no pathway
below −10 (RAF set aside). Mitotic G1 is a tuning pathway, so the held-out split should not
move. Its result is read on tuning and experimental, and reported as such.

## Predictions (written before the build)

1. **Kinase up or phosphatase down: fixed.** CDKN2A and CDKN1B knockout, PPP2R1A knockout,
   and CCND1 and CDK4 overexpression read UP at the E2F targets, as the curators expect.
   At least 60 of the about 85 such wrong curator cases are fixed.
2. **Kinase down or phosphatase up: not fixed.** CCND1 and CDK4 knockout, CDKN2A and
   CDKN1B overexpression, and PPP2R1A overexpression stay NORMAL (expected DOWN). At rest
   only φ₀ = 0.1 of p130 is phosphorylated, so removing the kinase raises the complex by
   at most about 3%. This is the known asymmetry of φ₀, not a defect of the rule. No
   φ₀ change is made here.
3. Mitotic G1 curator net is at least +40, below MP-BioPath's 146 of 187 because of (2).
4. Held-out net is 0. No pathway other than Mitotic G1 moves on canonical.

## Not in scope

- The 18 partner-release candidates (carrier class).
- The φ₀ asymmetry in (2).
- The faithful build (one more candidate, R-HSA-418996 in R-HSA-446728) is reported, not gated.

## Amendment 1 (2026-10-20, before the catalog build and the gated arm): the `recycle` edge

Implementing the rule showed that the registered design ("the logic network is
unchanged") cannot work. The solver solves a pool only when all its nodes are in one
cyclic component of the network. This cycle exists only in the pool table: free p130
leads nowhere. Adam's suggestion is the fix: "a special type of edge in the logic
networks so that deltasignal can handle them differently… we would still know what the
root node is."

- **Generator:** with the pool, `LNG_RECYCLE_POOLS` writes an edge `free copy → starting
  copy` of type `recycle` to `boundary_edges.csv`. It is a derived edge, like assembly.
- **DeltaSignal** (`src/core/recycle.jl`): a `recycle` edge is structural only.
  - It joins the cycle into one component (`index_reactions`, keyword `structural`).
  - It is never an input: `convert_to_reaction_network` drops it, and the solve removes
    it before any other rule reads the network.
  - The solve reports `recycle_edges`.
- **Benchmark:** a `recycle` in-edge never makes a node a non-root
  (`bench/test_pin_scope.py`), so a p130 knockout still pins the starting copy.
- **Tests:**
  - `test_cycles.jl`: two-base pool with recycle edges; no recycle edges means no pool
    solved; starting copies keep their own value.
  - `bench/test_pin_scope.py`.
  - LNG `tests/test_recycle_pools.py`: the p130 shape; a partner released unchanged; a
    free copy that leads on; a claimed node; flag default off.

**Disclosed: a correctness check on Mitotic G1 only** (`g1_check.jl`, after this
implementation, before the catalog build). A single-pathway test build of R-HSA-453279
with `LNG_RECYCLE_POOLS=1` (LNG branch `recycle-pools-081`) gives 2 recycle pools:
p130 with 4 resting copies and p-p130, and p107 with 2 and p-p107. It was solved for
every Mitotic G1 case with pins and readouts mapped by stable id:

| | right → right | wrong → right | right → wrong | wrong → wrong |
|---|---|---|---|---|
| curator (850) | 387 | **70** | **0** | 393 |
| experimental (88) | 24 | **10** | **2** | 52 |

The unchanged network reproduces `results/3e96e7c/dflt` on every case, and the new one
stays at baseline at rest (4 pools solved, 2 recycle edges, no node off baseline).

Fixed perturbations:
- CCND1 overexpression 17 and CDK4 overexpression 17;
- PPP2R1A knockout 10 and overexpression 8;
- CDKN2A knockout 10;
- RBL2 overexpression 3, CCND1 and CDK4 knockout 2 each, RBL1 knockout 1.

CDKN1B does not move. The gate, the arm and predictions 1–4 are unchanged. The check is
not the arm, and its figures are not the result.

## Build check (before the arm)

Variant build `20261020-1517_55b0f1b_rcp` (LNG `55b0f1b`, `--env LNG_RECYCLE_POOLS=1`),
compared by stable id with canonical `20261017-1246_5e663d9` (`compare_build.py`):

- **3 tables differ, all in R-HSA-453279:** 2 `recycle` edges added, 8 pool state rows added
  (p130: 4 resting copies + p-p130; p107: 2 + p-p107), and 228 transition rows added.
- **Nothing is removed**, and no other pathway differs. 91 pathways log 0 recycle pools.
