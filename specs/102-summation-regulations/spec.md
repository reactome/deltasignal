# specs/102: regulations stated in Reactome summations (pre-registration)

Adam, 2026-10-24: "close the gap" (experimental 635/849 against MP-BioPath's 643).

## Why

specs/101 T1 and T3: the agreed errors left in Mitotic G1 (RB1, AKT1 and SRC; 14 experimental, 118
curator cases wrong) have no route in Release 97's event graph. Each gene's reaction makes a
product that no reaction uses. Each reaction's own summation says what the product does.
MP-BioPath's network got these routes through about 1,800 hand edits (specs/034).

## The rule

Generator flag `LNG_SUMMATION_REGULATIONS` (default off). Within a pathway that contains the
source reaction, the source reaction's product becomes a regulator of each target reaction that
the pathway also contains. Every edge is written from the summation's own words, as a structural
rule, not as a list of cases:

| id | source reaction → regulator | sign | targets (structural) | summation |
|---|---|---|---|---|
| E1 | R-HSA-9018017 → RB1:E2F1/2/3:TFDP1/2 (R-HSA-68644) | − | reactions positively regulated by an entity containing E2F1, E2F2 or E2F3 (15 in Mitotic G1) | "RB1 binding inhibits transcriptional activity of E2F1/2/3:TFDP1/2 complexes" |
| E2 | R-HSA-198613 → p-T-CDKN1A/B (R-HSA-198605) | − | reactions with an input containing CDKN1A or CDKN1B and another input containing CDK2, CDK4 or CDK6 (R-HSA-8941915, 69562 in G1; 187934 in S Phase) | "retention in the cytoplasm, segregating these CDK inhibitors from cyclin-CDK complexes" |
| E3 | R-HSA-8942607 → CDK4/6:CCND:p-Y-Cip/Kip (R-HSA-8942622) | + | reactions whose catalyst contains CDK4 or CDK6 with no other active unit (R-HSA-69227, 1226094, 1226095) | "results in increased catalytic activity of CDK4 or CDK6" |

The table is written by a script from these three queries and committed with the generator.

## Circularity (read this before the numbers)

These three reactions were chosen BECAUSE their genes' cases are wrong. Gains on RB1, AKT1 and
SRC cases are in-sample and prove only that the edges carry the effect. The out-of-sample reading
is every other case that moves. Both readings will be reported. If adopted, results are reported
with and without the supplement, as a curation supplement to Release 97, not as Release 97.
The edges are also candidates for Reactome curation (Adam's call).

## Gate (as specs/077)

Curator held-out net ≥ 0, experimental net ≥ −5, no pathway below −10. Read with RAF set aside
too. Also: no case outside RB1, AKT1 or SRC perturbations may break unless traced and recorded.

## Predictions

1. E1: RB1 knockout reads the E2F targets UP in Mitotic G1. At least 20 of the 30 wrong curator
   RB1 cases fixed, and at least 3 of the 5 experimental ones (FBXO5, expected DOWN, stays wrong).
   Some curator NORMAL readouts break (PCNA, RRM2: MP-BioPath's network breaks them).
2. E2: AKT1 knockout lowers the E2F targets in Mitotic G1. At least 3 of its 5 wrong experimental
   cases fixed.
3. E3: SRC knockout lowers and SRC overexpression raises the E2F targets. At least 3 of its 5
   wrong experimental cases fixed.
4. Only Mitotic G1, S Phase and PIP3 activates AKT (the pathways holding a source reaction) move.

## Amendment 1 (before any measurement): how the regulator is wired

The first build (`20261024-2015_26ac7a2_sum`, never benchmarked) added each regulator as a curated
participant. Its set choices then entered the variant enumeration:
- every E1 target got three times the copies;
- the RBL1/RBL2 phosphorylations (R-HSA-1226095/1226094, E3 targets) went over the variant cap
  (120 → 4 copies);
- the capped fallback removed all 174 recycle transitions of the specs/081 p130/p107 pools.

A regulator the release does not curate should not reshape the reaction it regulates. LNG
67574bf (branch `summation-regulations-102`) instead adds it after enumeration, to every copy, as
ONE pool of its variants (`<stId>::pool`, the generator's existing OR pool node). Rules,
targets, gate and predictions are unchanged.
