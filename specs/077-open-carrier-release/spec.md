# specs/077: an open carrier loop's release is not supply when the form has another producer

**Status:** pre-registered 2026-10-19, before the arms; results in `research.md`. Canonical held-out +35 (p 0.0026) and experimental −5, but intrinsic apoptosis −26: **the gate fails** (no pathway below −10). Faithful build DSB +182, IFN α/β 0. Predictions 1–3 hold. Traced: the apoptosis loss is a displacement (BAD frees tBID from BCL-2), not a carrier release. Not adopted.  Amendment 1 (`fed_single`, keep displacer releases) pre-registered 2026-10-19; post hoc.

Adam, 2026-10-19: "ok" to tracing the DSB Repair loss on the faithful build
and pre-registering a targeted rule.

## The defect (traced)

Faithful build `20261018-0815_ddbd828_dsub`, against canonical
`20261017-1246_5e663d9`. Both use solver `bf4b688` with every current default.
The faithful build is 485 curator cases behind canonical: TP53 −223, DSB Repair −195 and RUNX1 −93.

DSB Repair loses 201 cases and gains 6. All the losses are MMEJ (alt-EJ)
perturbations: PARP1, PARP2, FEN1, XRCC1, LIG3 and POLQ, knockout and
overexpression, each reaching 12 HR and SSA readouts. RAD52 and BRCA2
overexpression also lose, and the NHEJ genes reach 1 or 2 readouts each. The
curators say no change. specs/070 §5 traced the cause:

- The signal leaves MMEJ at RBBP8 homotetramer (R-HSA-5684097). MMEJ's last
  step, "LIG3 ligates remaining SSBs in MMEJ" (R-HSA-5687675), releases it.
- In the faithful build the released copy is the free RBBP8 node. That node
  also feeds HR end processing ("MRN complex binds CDK2 and RBBP8",
  R-HSA-5684081).
- Free RBBP8 reads OR(SIRT6 deacetylation, MMEJ release), so stalling MMEJ halves it.

Neo4j (Release 97) confirms that this is how Reactome curates it. RBBP8 homotetramer
is one entity. Its producers are R-HSA-5687675 (MMEJ release) and R-HSA-5685953
("SIRT6 deacetylates RBBP8"). Its consumers are R-HSA-5687465 ("MRN complex
binds RBBP8", MMEJ's entry) and R-HSA-5684081 (HR). The network is faithful,
and the solver reads a carrier released unchanged as new supply.

The generator lists the loop in `carrier_loops.csv` as OPEN, because SIRT6 also
makes the free form. The default `closed_twin` acts only on closed loops.

## Why not `all` (probed, diagnostic arms `allp` and `allpd`)

`DS_CARRIER_LOOP_MODE=all` drops every listed release edge. Solver `4c0e0b9`
(specs/076 adoption commit, same code as main `8888f18`), against `rst` and `rstd`:

- On canonical (`results/4c0e0b9/allp`), curator net is −166 (134 fixed, 300 broken).
  Interferon α/β is −168 (2/170), intrinsic apoptosis −26 and FGFR1–4 −29. MET
  gains +20, EGFR +18 and IL-2 +8. Experimental net is −4.
- On the faithful build (`results/4c0e0b9/allpd`), curator net is +23. The difference
  is DSB Repair, **+182 (182 fixed / 0 broken)**. The IFN α/β, apoptosis and
  FGFR losses are identical to canonical.

The IFN α/β loss was traced here (it had never been traced in specs/049). Every JAK1 and
SOCS1 readout reads exactly 1.0× under `all`. Reactome curates the activated
receptor (R-HSA-909704, p-IFNAR1:p-TYK2) as a scaffold:

1. "Recruitment of STAT2 to p-IFNAR1" (R-HSA-909719) makes R-HSA-909696.
2. "Phosphorylation of STAT2" (R-HSA-909732) makes R-HSA-909697.
3. "Phosphorylation of STAT1" and "Release of p-STAT2:p-STAT1 dimer"
   (R-HSA-909722) return R-HSA-909704.

The generator lists all three receptor complexes as free forms of one open loop.
909696 and 909697 have no producer except a loop step, so under `all` they
lose every producer and the signal stops at the receptor. RBBP8 is the other
case: a free form that keeps a producer outside the loop (SIRT6).

## The rule

`DS_CARRIER_LOOP_OPEN` = `off` (default, unchanged) | `fed`.

- Under `fed`, with any carrier-loop mode except `off` and `all`, a free form
  of an OPEN loop loses its listed loop in-edges when it has an activating
  in-edge from a source that no loop lists for it. The release is then not a supply; its
  other producers set its level, as conservation says.
- A listed form whose every activating producer is a loop step keeps them all.
  These are the scaffold intermediates (IFN α/β), and the rule leaves them alone.
- Closed loops are handled by the mode as before. Only the generator-listed
  activating edge is dropped; depletion and inhibition edges stay.
- The solve reports `carrier_loop_open_rule` and `carrier_loop_open_edges`. The
  dropped edges are also counted in `carrier_loop_edges_dropped`. The bench prints
  `Open carrier loops:` and refuses a `fed` run that dropped nothing.
- Tests: `test/test_carrier_loops.jl` "open loops" (179 assertions). R is made by S and released
  by a binding step that makes RA from R and A, and RA feeds the readout T. Both R and
  RA are listed. Under `fed`, R reads S (3×) and T follows A (6× at A = 2,
  1.5× at A = 0.5). Under `all`, T reads 1× whatever A is.

## Structural count, before any arm

`count_open_fed.py` (this directory) counts the listed open forms with an outside
producer, and the edges `fed` drops (pairs listed per pathway, summed):

- canonical `20261017-1246_5e663d9`: 529 forms, **705 edges**. 354 listed forms
  have only loop producers and are untouched: 320 in Cell Cycle Checkpoints and
  26 in IFN α/β. The largest counts are Cell Cycle Checkpoints 162, R-HSA-168928 108 and
  R-HSA-1227986 53. IFN α/β has 26 edges (two copies of the release into R-HSA-909704).
- faithful `20261018-0815_ddbd828_dsub`: 720 forms, **1,018 edges**, of which DSB
  Repair has 59.

On both builds every outside producer the rule relies on is a reaction `output`
edge (529 of 529 and 720 of 720 forms), never a bridge or set edge.

## Arms

- `ocf`: canonical, `--server DS_CARRIER_LOOP_OPEN=fed`, against `results/bf4b688/rst`
  (current defaults).
- `ocfd`: faithful build, `--server DS_POOL_SUPPLY_REACH=material --server
  DS_CYCLE_CARRIER_RELEASE=reaction --server DS_CARRIER_LOOP_OPEN=fed`, against
  `results/bf4b688/rstd`.

Read with `holdout_report.py --set-aside RAF_MAP_kinase_cascade` (its
"baseline" is the `--cases` file), paired and on all cases valid in both, and
checked case by case.

**Gate (registered).** `fed` is adopted as the default if `ocf` against
canonical gives curator held-out net ≥ 0, experimental net ≥ −5, and no
pathway below −10 (RAF set aside). `ocfd` informs the faithful-network work
and is not gated.

**Predictions:**
1. `ocfd`: DSB Repair gains at least 150 curator cases net (`allpd`: +182).
2. `ocf` and `ocfd`: IFN α/β moves by at most 10 cases either way (`all`: −168).
3. `ocf`: fewer than 3 pathways lose more than 10 curator cases. Under `allp`, 3 do:
   IFN α/β, intrinsic apoptosis and G2/M. No prediction is made for the
   apoptosis and FGFR losses, because they have not been traced to fed or
   scaffold forms.

## Amendment 1: a release with a displacer stays (pre-registered 2026-10-19, before its arms)

Adam, 2026-10-19: "ok" to amendment 1. **Post hoc:** this rule was written after
seeing `ocf` fail the gate, so any held-out gain it shows is selected on the
held-out set and is reported with that caveat.

**Why.** `fed` cut R-HSA-139897 ("BAD displaces tBID from BCL-2 sequestration";
`research.md`), and every BAD, YWHAE and AKT1 perturbation in intrinsic apoptosis
went to 1.0×. That reaction is a displacement, not the return of a carrier: the
bound share is what changes, and a second input (BAD) drives it. The conservation
argument behind `fed` (bound forms negligible, so free = total) does not hold there.

**Rule.** New value `DS_CARRIER_LOOP_OPEN=fed_single`. As `fed`, except that a
listed release is kept when its source reaction has two or more distinct
non-cofactor `input` parents (activating edges of type `input`; cofactors are
`cofactor_uuids(network)`). `fed` is unchanged and still reproduces `ocf`.
Structural count (`count_displacers.py`): this keeps 23 of the 705 canonical edges
and 49 of the 1,018 faithful ones (24 in R-HSA-170834); the RBBP8 release
(R-HSA-5687675) and every other DSB release are single-input and are still dropped.

**Arms** (solver = the commit that adds `fed_single`):
- `ocs`: canonical, `--server DS_CARRIER_LOOP_OPEN=fed_single`, against `results/bf4b688/rst`.
- `ocsd`: faithful build, `--server DS_POOL_SUPPLY_REACH=material --server
  DS_CYCLE_CARRIER_RELEASE=reaction --server DS_CARRIER_LOOP_OPEN=fed_single`,
  against `results/bf4b688/rstd`.

**Gate:** unchanged. `ocs` against canonical: curator held-out net ≥ 0,
experimental net ≥ −5, no pathway below −10 (RAF set aside). `ocsd` is not gated.

**Predictions:**
1. `ocs`: intrinsic apoptosis net within ±5 of `rst` (the 30 BAD/YWHAE/AKT1 cases return).
2. `ocsd`: DSB Repair still gains at least 150.
3. `ocs`: curator held-out net at least +20 (most of `ocf`'s +35, which is MET, EGFR and IL-2).
4. `ocsd`: RUNX2's −21 (CBFB, ESR1 knockouts) is not predicted to change; it is not traced.
