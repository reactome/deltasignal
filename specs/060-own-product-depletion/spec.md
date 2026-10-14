# specs/060: an own-product depleter may suppress, not de-repress (re-measured on the canonical build)

**Status:** registered 2026-10-14, before the arm.

## Why

specs/059 R2 (arm `rtwin58s`) added `DS_DEPLETION_OWN_PRODUCT=suppress_only`
to two other rules and, against `rtwin58`, moved only Transcriptional
Regulation by TP53: tuning +111 (111 / 0; AKT1 KO +38, AKT2 KO +38, MDM4 KO
+32, BRCA1 KO +3) and experimental +6 (6 / 0, p 0.031). The AKT → MDM2 →
TP53 knockouts are the case specs/016 traced and every variant since missed.

The switch exists since specs/016 (`full` default): a depletion edge whose
source complex is built from its target (an "own-product depleter", e.g.
MDM2:TP53 -| TP53) may suppress, but its factor is capped at 1, so it cannot
de-repress. It was measured on a much older catalog and not adopted. What it
loses (specs/016): de-repression when the complex fell because of its OTHER
partner (TP53 KO → less MDM2:TP53 → more free MDM2).

This measures it alone, on the build every default is judged on, as specs/055
did for inert drugs. The TP53 gains were seen first (in `rtwin58s`), so they
are not evidence for the rule here; the held-out axis is.

## Arm

| arm | build | solver | against |
|---|---|---|---|
| control | canonical | defaults (`results/64d23cb`) | — |
| `ownsup` | canonical | `--server DS_DEPLETION_OWN_PRODUCT=suppress_only` | control |

Both axes, `--pair-on-valid` (pins do not change), RAF set aside.

## Gate (as specs/055–059)

Curator held-out net > 0 at p < 0.01; experimental net loss (RAF set aside)
not significant; no pathway losing more than 10 on either axis; every
experimental loss outside RAF traced. Both readings reported.

## Predictions

- Tuning: the TP53 gains reappear (about +111), and experimental about +6.
- Held-out: no sign predicted. If held-out is about 0, the gate fails as
  written (it asks for a held-out gain); the tuning and experimental gains
  are then reported for Adam to weigh, not adopted by this spec.
