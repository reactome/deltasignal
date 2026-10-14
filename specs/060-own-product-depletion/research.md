# specs/060 research

## R1. Arm `ownsup`: held-out unchanged; TP53 only; the gate fails

`scripts/run_arm.sh ownsup --server DS_DEPLETION_OWN_PRODUCT=suppress_only`
at DS `cc9ab80` (solver = `64d23cb` plus the specs/059 rule, off), canonical
build, results `results/cc9ab80/ownsup`. Converged: curator 1,651 of 1,725
(canonical 1,652), experimental 214 of 244 (215).

| axis | net | fixed / broke | p |
|---|---|---|---|
| curator held-out | **0** | 0 / 0 | — |
| curator tuning | +46 | 48 / 2 | <1e-6 |
| experimental (= RAF set aside) | +4 | 4 / 0 | 0.125 |

Every change is in Transcriptional Regulation by TP53: MDM4 KO +34, AKT1 KO
+7, AKT2 KO +7, BRCA1 KO −2 (curator); MDM4 KO +2, AKT1 +1, BRCA1 +1
(experimental). The rule moves no held-out case in 70 pathways.

**Gate: fails** (held-out net 0). Not adopted.

Predictions: tuning "about +111" and experimental "about +6" — **not met**
(+46, +4). Alone, the rule recovers MDM4; the larger AKT1/AKT2 gains
(+38 each in specs/059's `rtwin58s`) need the specs/058 pins and the
specs/059 twin edges as well.

### What this means

- `suppress_only` is inert outside TP53 on this build and right inside it
  (48 / 2 curator, 4 / 0 experimental), on a tuning pathway. The gate asks
  for a held-out gain, which a rule that touches no held-out case cannot
  show. Whether to adopt it on tuning and experimental evidence alone is
  Adam's call; this spec does not.
- For the record, the combination `rtwin58s` (specs/059 R2) is held-out +21
  (p 0.0025), tuning +112, experimental +7 (7 / 0, p 0.016), and fails only
  on TGF-beta −12 (the SMAD2/3 knockouts, not fixed by this rule).
