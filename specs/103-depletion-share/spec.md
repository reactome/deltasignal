# specs/103: several consumers share one pool (pre-registration)

Adam, 2026-10-25: "I want you to figure out the last two" (experimental 641/849 against
MP-BioPath's 643).

## Why (traced, `trace_mdm4b.jl`, canonical `20261024-2211_7f4b28a`, MDM4 knockout in TP53)

Free p-S166,S188-MDM2 [nucleoplasm] (R-HSA-6793669) is consumed by two bindings, "MDM2 forms homo-
or heterodimers" (R-HSA-6804741), so it carries two depletion edges: from the MDM2:MDM4
heterodimer (6804932 family) and from the MDM2 homodimer (6804933). Under the MDM4 knockout the
heterodimer reads 0. The solver multiplies the depleters' factors, `H_dep = Π b/x_k`, so the
heterodimer alone de-represses free p-MDM2 by the full ceiling (10×), as if it consumed the whole
pool. Free p-MDM2 reads 8.3×, the homodimer 8.3–100×, the ligase pool 6804885::pool 2.73×, and
TP53 [nucleoplasm] 0.56: the TP53 targets read DOWN where curators and experiments say UP (5
agreed experimental errors; specs/101 T2).

Mass action says otherwise. A substrate S supplied at rate s and consumed by n reactions with rates
k_i·E_i sits at S = s / Σ k_i E_i. With equal baseline shares (the same assumption as the OR mean),
its fold is S/S₀ = 1 / mean(E_i/E_i₀). The product treats every consumer as the only one. With
the mean, losing one of two consumers at most doubles S.

## The rule

Solver flag `DS_DEPLETION_COMBINE` = `product` (default, today's solver) | `share`. Under `share`,
a reaction's depletion factor is `H_dep = (b+ε) / (mean_k x̃_k + ε)` with the same clamp
[`DS_DEPLETION_H_MIN`, `DS_DEPLETION_H_MAX`]. x̃_k is the depleter's value, except where the own-
product rules (specs/060–061) cap that depleter's factor at 1; there x̃_k = max(x_k, b). One
depleter: identical to `product`. Closure reads (specs/018) unchanged.

Footprint: 44 nodes with two or more depletion in-edges, in 15 pathways (EGFR 8, DSB repair 6,
DNA repair 6, PIP3 5, …; TP53 1).

## Gate (as specs/077)

Curator held-out net ≥ 0, experimental net ≥ −5, no pathway below −10; read with RAF set aside
too.

## Predictions

1. MDM4 knockout in TP53: free p-MDM2 reads ≤ 2× (was 8.3×), the ligase pool falls below 1, and at
   least 3 of the 5 agreed experimental MDM4 knockout errors are fixed.
2. Experimental net ≥ +2 (643/849, MP-BioPath's figure).
3. Only pathways that hold a node with two or more depletion in-edges move.
4. specs/011's EPS15 case (two EGFR:CBL complexes at 100× depleting GRB2-1) still reads GRB2-1 at
   the floor: under `share` the mean is 100×, so the factor is 0.01, clamped to 0.1 as before.

## Circularity

The rule came from tracing one agreed error class (MDM4). The out-of-sample reading is every other
case that moves.
