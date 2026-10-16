# specs/058: a form fed only by depletion is still a root

**Status:** **Adopted as the bench default 2026-10-15** together with specs/059 and specs/061 (specs/061 Decision); alone it was not adopted: registered 2026-10-14, before any code or arm. **Not adopted** (research.md R1): held-out +9 (33 / 24, p 0.29), experimental +1; the gate fails. Gains in MET, CBL, PTCH, BRCA1/BARD1, MDM2; all losses are SMAD2/3 in TGF-beta (−24), whose plain cytosolic form becomes a root that does not reach the signaling route.

## Why

Found tracing Signaling by MET (specs/057): plain MET (R-HSA-419603) has no
producer in the pathway, but its one in-edge is a **depletion** edge, so the
root rule ("no incoming edge") does not count it, and MET is pinned only at
a side form (p-Y1003,4Y-MET). MET then never perturbs the main route.

Depletion edges are our derived edges (phosphatase and ubiquitin-ligase
consumption, generator `_emit_substrate_depletion_edges`). MP-BioPath's
networks have none, so under its protocol ("perturb root inputs that are, or
contain, the gene") such a form is a root. The generator already treats
depletion as not production when it decides what is a root (LNG #97,
`_existing_upstream`), and the specs/024 fallback already reads depletion and
dissociation in-edges as derived. Only the bench's first test counts them.

## Census (canonical build `20261012-1401_01206ca`)

183 entity nodes have depletion as their only in-edges. Joined to the
perturbed genes (via Neo4j geneName), they belong to 28 (pathway, gene)
pairs with cases, e.g. MET 0/20 and LRIG1 0/20 (MET), PTCH 8/26 (Hedgehog),
CBL 1/12 (EGFR), TRIM33 4/14 (TGF-beta), MDM4 52/102 (TP53, tuning), and
pairs scoring perfectly now (AKT1 in PIP3 28/28 and exp 28/28, DLL1, JAG1,
SMAD2/3 TGF-beta, RUNX2). Some may already be pinned through the specs/024
fallback.

## The rule (bench, `DS_PIN_ROOT_IGNORE=depletion`, default unset)

Under `root` and `root_cycle`, a node is a root when it has no incoming edge
**other than depletion**. Everything else is unchanged: the specs/024
fallback still applies to a gene with no root; dissociation still counts as
production. The Protocol line prints the setting.

## Arm

| arm | build | solver | bench | against |
|---|---|---|---|---|
| control | canonical | defaults (`results/64d23cb`) | default | — |
| `deproot` | canonical | defaults, this branch | `--bench DS_PIN_ROOT_IGNORE=depletion` | control |

Pin sets change by design: `holdout_report.py --pair-on-valid`, and the cases
valid in one arm only are reported. Both axes.

## Gate

The specs/049 amendment 5 faithfulness gate with RAF set aside (as specs/055,
056, 057): curator held-out net > 0 at p < 0.01; experimental net loss (RAF
set aside) not significant; no pathway losing more than 10 on either axis;
every experimental loss outside RAF traced. Both readings reported. If it
passes, the setting is recommended as the pin default (Adam's call); else
nothing changes.

## Predictions

- MET overexpression goes right (10): both MET forms pinned at 80× read UP
  (traced, specs/057). MET knockout stays wrong (the GRB2 depletion still
  rails the recycling loop; specs/057).
- No sign predicted for the rest; the census pairs include perfect scorers
  that could lose.
