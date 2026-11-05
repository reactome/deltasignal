# specs/138 research

## R1. The engaged cases, traced (2026-11-05)

Drug sets in which every drug engaged its target, each set's own markers left out: 1,571 cases, 159 wrong (model 0.899
accuracy, always-NOT 0.912). They split as:

| group | cases | what it is |
|---|---|---|
| missed changes with no route | 81 of 83 missed | the readout is not reachable from the drug's target in the scored pathway: Cyclin B1 falling (cell-cycle arrest, 16), MAPK pT202 falling under PI3K/AKT drugs (crosstalk, 12), total AKT falling under Temsirolimus (12), Caveolin rising (9) |
| Temsirolimus → AKT pS473 / pT308 | 18 | predicted DOWN, measured UP, in PIP3, VEGF and CD28 signalling |
| AKT inhibitor in *Mitotic G1* | 21 | RB pS807, cyclin E1, PCNA predicted UP at the rail (34–100×), measured DOWN or NOT: the p27 assembly-factor inversion recorded in specs/102 |
| other false change calls | ~40 | small and spread: SRC pY416 under ZSTK474 (10), TSC2 pT1462 under AKT drugs (11, measured unchanged although TSC2 T1462 is a direct AKT site), GSK3 under Temsirolimus pairs |

**Temsirolimus, checked in Neo4j before any rule.** The bench removes all of MTOR, so mTORC2's S473 phosphorylation goes
and AKT falls. A rapalog blocks mTORC1 only, and in cells AKT phosphorylation rises because mTORC1 → S6K → IRS1 negative
feedback is relieved. Reactome does not curate that feedback: the reactions S6K catalyses are S6, eIF4B, eIF4G and eEF2K
phosphorylation only, and no reaction phosphorylates IRS1 on serine. So an mTORC1-only Temsirolimus would move these 18
from DOWN to NOT, never to UP: no score changes. Not registered. **Curation-review item:** the S6K → IRS1 (S307/S636)
negative feedback, and its loss under rapalogs, is absent from Reactome.

**What this says.** On the cases korkut can judge, the model's remaining errors are mostly changes it has no route to in
one pathway (cell-level responses and crosstalk), one missing curated feedback, and one known curated-role inversion.
No structural rule is indicated.
