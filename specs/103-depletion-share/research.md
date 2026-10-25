# specs/103 research

## Arm `share` (canonical `20261024-2211_7f4b28a`, results `results/b8616ac/share`; control `results/23a526c/dflt102`)

`scripts/run_arm.sh share --server DS_DEPLETION_COMBINE=share`.

| axis | net | fixed / broke | p | RAF set aside |
|---|---|---|---|---|
| curator held-out | 0 | 2 / 2 | 1 | 0 |
| curator tuning | +13 | 29 / 16 | 0.07 | +13 |
| experimental | 0 | 4 / 4 | 1 | 0 |

By pathway (curator): TP53 +26, PIP3 **−13**, EGFR 0, RUNX2 0. **Gate fails** (PIP3 below −10).

### Predictions

1. **Mechanism holds, case count fails.** Trace (`trace_mdm4b.jl`, `DS_DEPLETION_COMBINE=share`):
   under MDM4 knockout free p-MDM2 reads 1.45× (was 8.3×), the ligase pool 6804885::pool 0.80 (was
   2.73), TP53 [nucleoplasm] 1.10–1.15. Curator MDM4 knockout: 9 fixed. Experimental: 1 of the 5
   agreed MDM4 errors fixed (6803386, 0.24 → 1.42); 140215 and 6798008 read 1.1485, just under
   the 1.15 cutoff (were 0.56).
2. **Fails.** Experimental net 0 (641/849).
3. **Holds.** Only TP53, PIP3, EGFR and RUNX2 move, each holding a shared-depletion node.
4. **Holds** in the case it names. The EPS15 knockout at 179882 does move (10.0 → 1.10, broken), but
   through the de-repression side, not the floor.

### What moved

- **VAV1 knockout in PIP3: curator −13, experimental +2 / −1.** PTEN [cytosol] (R-HSA-199420) has
  four ligase depleters: p-S109-MKRN1, PTEN:TRIM27, and the XIAP,NEDD4 and NEDD4,STUB1,WWP2,XIAP set
  pools. VAV1 knockout acts through AKT → MKRN1 (the specs/082 loop). Under `product` the MKRN1 drop
  alone de-represses PTEN; under `share` it is a quarter of the consumption, and the readouts go
  0.17 → 0.87. The curators say DOWN on all 13. The experiments say NORMAL on 2 of 3. Equal shares
  are an assumption Reactome cannot check (no rates).
- **ATM and CHEK2 knockout in TP53: curator +14 / −3, experimental −3.** They act on TP53 through
  MDM2/MDM4 phosphorylation. They now read near NORMAL, as the curators say; the experiments say DOWN.
- **AKT1, AKT2 knockout in TP53: curator +6, experimental +1** (threshold crossings, 1.10 → 1.16).
- EGFR: SRC +1, EPS15 −1. RUNX2: PPM1D +1 / −1.

**Not adopted** as registered (gate fails). `DS_DEPLETION_COMBINE` stays in the solver, default
`product`, as the record. A per-depleter weight (MKRN1's share of PTEN turnover) would need
measured rates; none is fixed a priori, and fitting one to these cases would be circular.
