# specs/055: drug-derived nodes held at baseline, re-measured on the canonical build

**Status:** registered 2026-10-13, before the arm ran. **Result: gate passes (research.md R2); thin evidence (2 perturbations outside RAF, 0 broken). Adoption is Adam's call.** No code change: the
rule is `DS_DRUG_MODE=inert` (specs/032), unchanged since then.

## Why now

specs/032 held drug-derived nodes (a Reactome `Drug`, a complex containing
one, a set that is all drugs) at baseline, for a cell without the drug. It
gained curator held-out **+30** (30 fixed / 0 broke, p 1.9e-9, 4 pathways)
and was declined on experimental **−5, all of it RAF**. specs/049 measured it
again under carrier-loop arms (`kdrug` −5 → +2 with RAF set aside, held-out
+8; `scd` 0 → +7). specs/054 D1 (Adam, 2026-10-13) sets RAF aside from rule
decisions because its readouts are a bistable switch. The objection to
specs/032 was RAF, so the rule is re-measured, once, on the build every
default is now judged on.

## Arms

| arm | build | solver | against |
|---|---|---|---|
| `cj0` (existing) | canonical `20261012-1401_01206ca` | defaults, solver `a807851` | — |
| `dinert` | same | defaults + `DS_DRUG_MODE=inert` | `cj0` |

`dinert` runs by `scripts/run_arm.sh dinert --server DS_DRUG_MODE=inert` from
this branch. `src/` and `test/` are identical to `a807851` (`git diff
a807851 -- src test` is empty), so `cj0` is its control. Both axes.

## Predictions

- **D-P1.** Curator held-out net > 0. Expected range +5 to +40: 032 gave +30
  before variant nodes and `closed_twin`; `kdrug` gave +8 on the K build.
- **D-P2.** Experimental with RAF set aside: |net| ≤ 3, not significant.
- **D-P3.** RAF moves (direction not predicted; set aside, reported).
- **D-P4.** The solve reports `drug_rule = inert` and `drugs_held` > 0 in
  every scored pathway that has drugs.

## Decision rule

The specs/049 amendment 5 faithfulness gate, read with RAF set aside
(specs/054 D1), against `cj0`:
- curator held-out net > 0 at McNemar p < 0.01;
- experimental net loss (RAF set aside) not significant (p ≥ 0.05);
- no pathway (RAF set aside) losing more than 10 cases on either axis;
- every experimental loss outside RAF traced to a named mechanism.

Both readings are reported (with RAF and without). Held-out is SELECT
(specs/051). If the gate passes, `DS_DRUG_MODE=inert` is recommended as the
default, and the change is Adam's call. If it fails, nothing changes.
