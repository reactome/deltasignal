# specs/102 research

## Arm `sum` (build `20261024-2055_67574bf_sum2`, results `results/5a7616a/sum`; control `20261024-1634_f2c6a71/results/bc67b32/dflt101`)

Build: only Mitotic G1 and S Phase tables change. No reaction gains or loses copies, the p130/p107
recycle pools are intact, and 72 dissociation sinks of RB1:E2F and the E3 product go (they are no
longer dead ends).

| axis | net | fixed / broke | p | RAF set aside |
|---|---|---|---|---|
| curator held-out | 0 | 0 / 0 | 1 | 0 |
| curator tuning (Mitotic G1 only) | **+173** | 206 / 33 | < 1e-4 | +173 |
| experimental | **+6** | 11 / 5 | 0.21 | +6 |

Experimental 635 → **641/849** (75.50%; MP-BioPath 643). Curator all 88.48% → 89.20%. Converged:
1,724 / 1,725 curator and 244 / 244 experimental (unchanged). **Gate passes.**

### Predictions

1. **Holds.** RB1: curator 30 fixed / 2 broken; experimental 4 fixed (CCNA1 157443, TK1 197984,
   68891, 73500) / 3 broken (NORMAL readouts 68439, PCNA 68487, 73639, as MP-BioPath's network
   breaks them). FBXO5 stays wrong.
2. **Fails.** AKT1: no experimental case fixed; curator 2 broken. AKT1 now has a route and it
   inverts. Its knockout removes p-T-p27/p21, so more p27 binds CDK4/6:CCND (R-HSA-8941915). In
   Release 97 that complex is the substrate of the activating tyrosine phosphorylation
   (R-HSA-8942607, E3's source), which raises CDK4/6 activity. That is the curated
   assembly-factor role of p27 (specs/081's CDKN1B inversion). The other E2 target, CycE:CDK2
   inhibition (R-HSA-69562), makes a complex that nothing uses.
3. **Holds.** SRC: 3 experimental fixed (overexpression); curator 32 fixed.
4. **Holds.** Only Mitotic G1 moves (S Phase's tables change, its predictions do not).

### Out-of-sample reading (perturbations other than RB1, AKT1, SRC)

Fixed: ABL1 32, JAK2 32 and PTK6 31. All three are the other tyrosine kinases named in E3's
summation and catalysts of the same reaction, not chosen by error. Also CDKN2A
overexpression 14, CCND1 10, CDK4 8, CCNE1 5, MYC/MAX 4, CDKN1B 2, RBL1/RBL2 4, PPP2R1A 2.
Experimental CCND1 1, CCNE1 3.

Broken, each traced:
- **CCNE1, 24 curator / 2 experimental.** Cyclin E:CDK2 phosphorylates RB1 (R-HSA-188390, curated
  in this pathway), so CCNE1 now moves every E2F target (knockout 0.02–0.15, overexpression
  10–90×). The curators call these NORMAL. The experiments agree with the change on 3 of the 5
  CCNE1 readouts that move. A curator–biology disagreement, not a defect.
- **CDKN1B, 2, at readout R-HSA-8942622 (E3's regulator itself).** p27 knockout now reads it 1.33
  (was 0.67): it is no longer dissociated into free subunits, and the readout reads the complex.
- **CDKN2A knockout, 1 (8848441 1.42 → 1.00), and PPP2R1A overexpression, 2 (157443 0.79 →
  0.89).** Threshold shifts at the 0.85/1.15 cutoffs.

### Circularity

The +173 / +6 are mostly in-sample: Mitotic G1 is a tuning pathway, and RB1 and SRC were chosen
by error. ABL1, JAK2 and PTK6 (95 curator cases) were not chosen and move as the summation says.
No held-out pathway contains a source reaction, so there is no held-out reading.

**Adopted** (Adam, 2026-10-24: "yes adopt it"), as registered, E2 included: `LNG_SUMMATION_REGULATIONS` on
by default (LNG PR #127). Headlines are reported with the supplement and with Release 97 alone.
