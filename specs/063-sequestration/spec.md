# specs/063: sequestration — a dead-end complex lowers its free partner

**Status:** registered 2026-10-15, before any code or arm. Groundwork in
research.md R0–R3 (no structural or curated signature separates a
sequestering complex from a terminal product), so the question is put to
measurement: a narrow and a broad rule.

## Why

The largest group of remaining curator misses (348 held-out cases; 274
inhibitor-like) is a gene whose route ends at a complex nothing in the
pathway uses: follistatin:activin, USP18:IFNAR2, IL13RA2:IL13, CDK6:RUNX1.
Reactome curates these as binding reactions; the biology is that the
binding removes the partner from circulation. The model propagates fold
changes and never lowers the free partner (no mass conservation outside
pools, carrier loops, phosphatase and ubiquitin depletion). Consumption by
every binding step was refuted on 2026-07-16; this is restricted to complexes
that go nowhere.

## The rule (solver, `DS_SEQUESTRATION`, default `off`)

A **dead-end binding** is a reaction R with at least two non-cofactor entity
inputs and exactly one non-cofactor output C, where C is an entity whose
every out-edge goes into a dissociation sink (nothing in the pathway uses
it). A **partner** is an input of R with another use (an out-edge into
anything but R or a dissociation sink).

- `deadend_sole` (primary): only dead-end bindings with at least one input
  that has no other use (the would-be sequesterer) and exactly one partner.
  C gains a depletion edge C -| partner.
- `deadend` (broad, the record): every dead-end binding; C -| each partner.
- `off`: nothing changes (byte-identical).

The added edges are ordinary depletion edges (divide form, bounded by
`DS_DEPLETION_H_MIN`), added after the carrier and twin rules. Cofactors are
neither inputs counted nor partners. An observation still wins. The solve
reports `sequestration_rule` and `sequestration_edges`. Census on the
canonical build (research.md R1): `deadend` 15,429 edges on 6,295 complexes in
77 pathways; `deadend_sole` 1,477 complexes in 56 pathways, covering 8 of the
24 target complexes.

## Arms (canonical build `20261012-1401_01206ca`, against the control `results/64d23cb`)

| arm | solver |
|---|---|
| `seqsole` (primary) | `--server DS_SEQUESTRATION=deadend_sole` |
| `seqall` (record) | `--server DS_SEQUESTRATION=deadend` |

Measured against the current defaults, not the specs/061 candidate (not
adopted). Note for a later combination: under any
`DS_DEPLETION_OWN_PRODUCT` other than `full`, these edges are own-product
depletions and could not de-repress, so a sequesterer knockout would not
raise its partner.

`--pair-on-valid` is not needed (pins are unchanged); RAF set aside, both
readings reported.

## Gate (as specs/055–061)

Against the control: curator held-out net > 0 at p < 0.01; experimental net
loss (RAF set aside) not significant; no pathway losing more than 10 on
either axis; every experimental loss outside RAF traced. Convergence counts
reported (the edges close loops C -| P -> R -> C).

## Predictions (before any code)

- `seqsole`: gains in the covered targets: USP18 (56), IL13RA2 (24), SOCS6
  and SH2B3 (16 each), DLK1 (14), BH3:BCL2 (14, the double negative), CER1:BMP
  (4), SOCS1 (2): at most about 146 if every case flips, more likely half,
  since de-repression is capped at 10× and readouts may sit behind other
  inputs. Losses: unknown, from the other ~1,450 complexes; I expect some, in
  pathways whose terminal signalling complex is a dead end (ligand binding
  depletes free receptor). I give it even odds of passing.
- `seqall`: larger gains (follistatin, LEFTY, CDK6, PMEPA1 join) and larger
  losses; I expect it to fail the held-out gate, as the July rule did.

## Amendment 1 (2026-10-15, after R4, before any amended code): the sequesterer depletes its partner directly

R4: the complex is built from the partner, so C -| B is an own-product
depletion and the partner is de-repressed by its own fall (DSB −414 from two
edges). Under multiplicative AND the complex's fold is the sequesterer's fold
times the partner's; the part that is not the partner's own is the
sequesterer's. So:

- `DS_SEQUESTRATION=sole_direct` (new value; `deadend_sole` and `deadend` stay
  as recorded): the same dead-end bindings as `deadend_sole`, but the edge is
  A -| B from each sole-use input A (the sequesterer) to the partner B, not
  from the complex. No loop is closed; B's own changes do not feed back.
- Arm `seqdirect`, same control, gate and readings as above.

**Predictions:** the targets as in `seqsole` (USP18 +56, SOCS6, SH2B3,
IL13RA2, CER1 gain; SOCS1 may still lose), and the DSB and most RUNX1 losses
gone, since the false changes came from the feedback. Losses remain where a
sole-use input is not a sequesterer at all (a terminal ligand–receptor
complex: ligand overexpression would now lower free receptor). I expect
held-out between −200 and +60; it passes only if RUNX1 is mostly feedback.
