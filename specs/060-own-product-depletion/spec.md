# specs/060: an own-product depleter may suppress, not de-repress (re-measured on the canonical build)

**Status:** Superseded by specs/061, whose `suppress_only_source` is the default since 2026-10-15. registered 2026-10-14, before the arm. **Gate fails** (research.md R1): held-out 0 (no case moved), tuning +46, experimental +4 (p 0.125), all in TP53. Not adopted; Adam's call whether tuning and experimental evidence alone justify it. Amendment 1 (`suppress_only_entity`, R2): both arms fail (experimental −17 / −19, the MDM2 knockout); it fixes SMAD2/3 but cannot tell which partner fell.

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

## Amendment 1 (2026-10-14, after R1, before any code): own product by entity

specs/059 R2 traced why `suppress_only` missed the SMAD2/3 knockouts: the
own-product test is per node copy (the depleter is reached from THIS target
copy along substrate/assembly chains). SMAD2:SMURF2 is built from SMAD2's root
copies and depletes its recycled copy, so the test says "not own product",
and the knockout (depleter at 0) de-represses the recycled copy.

**Rule:** a new value `DS_DEPLETION_OWN_PRODUCT=suppress_only_entity`: as
`suppress_only`, and a depletion edge also counts as own product when the
depleter's entity CONTAINS the target's entity (the bundle's containment
table, as specs/022 reads it). Default stays `full`.

**Arms:**

| arm | solver | bench | against |
|---|---|---|---|
| `rtwin58e` | `DS_RECYCLED_TWIN=on`, `DS_DEPLETION_OWN_PRODUCT=suppress_only_entity` | `DS_PIN_ROOT_IGNORE=depletion` | canonical, and `rtwin58s` |
| `ownent` | `DS_DEPLETION_OWN_PRODUCT=suppress_only_entity` | default | canonical, and `ownsup` |

Same gate. **Post hoc**: chosen after seeing the SMAD2/3 loss.

**Predictions:**
- `rtwin58e`: the 12 SMAD2/3 knockout cases recover, so TGF-beta ≥ 0 and
  held-out about +33 against canonical; tuning and experimental as `rtwin58s`
  (about +112, +7), unless the entity test also flags depleters elsewhere.
- `ownent`: no sign predicted beyond `ownsup`'s TP53 gains; the entity test
  widens the rule to every depleter containing its target, catalog-wide.
