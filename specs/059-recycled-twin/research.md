# specs/059 research

## R1. Arms `rtwin` and `rtwin58` against the canonical results

Both by `scripts/run_arm.sh` at DS `3d5d0e8` on the canonical build, results
`builds/20261012-1401_01206ca/results/3d5d0e8/{rtwin,rtwin58}`; Protocol lines
checked (`rtwin58` prints `DS_PIN_ROOT_IGNORE=depletion`). Scored with
`--pair-on-valid --set-aside RAF_MAP_kinase_cascade`; no case valid in one arm
only.

| arm | held-out net | fixed / broke | p | tuning | experimental (= RAF set aside) |
|---|---|---|---|---|---|
| `rtwin` | 0 | 0 / 0 | — | −1 (1 / 2) | −1 (0 / 1, MDM2 OE) |
| `rtwin58` | **+21** | 33 / 12 | **0.0025** | +1 (2 / 1) | +1 (1 / 0, MDM2 OE) |
| `rtwin58` vs `deproot` | +12 | 12 / 0 | 0.0005 | −1 | — |

Headline (`rtwin58`, in release): held-out 89.33% → 89.45%.

`rtwin58` by pathway: MET +10, EGFR +11 (CBL), Hedgehog +8 (PTCH), DSB +4
(BRCA1, BARD1), Checkpoints +2 (MDM2 OE, tuning), TP53 −1 (CHEK2 KO, tuning),
**TGF-beta −12** (SMAD2 KO −6, SMAD3 KO −6).

**Gate.** `rtwin`: fails (held-out 0). `rtwin58`: held-out +21 at p 0.0025
and experimental +1 pass, but **TGF-beta loses 12, over the 10-case limit:
the gate fails.** Not adopted.

Predictions: `rtwin` unchanged on SMAD2/3 (held); `rtwin58` keeps the
specs/058 gains (held) and recovers "most of" SMAD2/3's −24: it recovered the
12 overexpression cases, not the 12 knockouts.

### The SMAD2/3 knockouts: a self-depleter de-represses

Under `rtwin58` the SMAD2 KO reads its TGF-beta readouts UP (50.5×, 13×,
10×; p-SMAD2:SMAD4 targets), where canonical read them DOWN (0.5×, 0.73×).
The pins are the two root copies (5f6c57d9, c854b46a). Mechanism: SMAD2 at 0
zeroes SMAD2:SMURF2, which depletes the recycled SMAD2 copy (7632c373);
with its depleter at 0 the copy is de-repressed to the 10× ceiling, and the
regeneration cycle (now anchored by the twin edge, which reads 0) carries it
up. The depleter is built from the target: a cell with no SMAD2 has no
SMAD2:SMURF2, and that absence cannot mean more SMAD2. The existing switch
`DS_DEPLETION_OWN_PRODUCT=suppress_only` (specs/016: an own-product depleter
may suppress, not de-repress) is the rule for exactly this.

`rtwin`'s one experimental loss (MDM2 OE, Cell Cycle Checkpoints) is fixed
again under `rtwin58`; not traced.
