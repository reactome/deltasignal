# specs/153 research: rating the unrated arrows on decryptM

## S1 — coverage (free screen): 6 of 16 arrows coverable

`s1_coverage.py` / `s1_coverage.tsv` (arrow entities expanded to Neo4j leaf
descendants; reachability inside the sound Release-98 source networks).
Coverable: the three p-AKT arrows (ERBB2 9 drugs / Intrinsic Apoptosis 1 /
Mitotic G1 2) and three unrated ones — PTEN Regulation -> PTEN -> PIP3
(AZD8055, Dactolisib, MK2206), Hedgehog 'off' -> GLI3R -> RUNX2 (Tideglusib
via GSK3B), DSB Response -> p-S456-ABL1 -> RUNX1 (Dasatinib, Imatinib,
Nintedanib). Attribution unambiguous: in every covered centre exactly one
arrow's source carries the covering drug's targets. The other ten arrows have
no decryptM drug upstream: closed as unratable on the current record.

## S2 — the measurement: decryptM CANNOT rate them (a clean negative)

Arms `r98fix_dm153` (control, sound 92-build) and `arrows98fix_dm153`
(composites), both on `cases_153_centres.tsv` (658 cases, 9 centres, 15
drugs; identical join of 658/658). Across all covered (centre, drug) pairs
there are **103 control-untouched sites, and the composites move exactly 2**
(both Dasatinib: one in ERBB2, one in Mitotic G1; both bridged "down" calls,
both measured "not"). Every other bridged fold stays exactly 1.0 — including
all of MK2206's (an AKT knockout feeding the p-AKT bridges) and all 29
untouched PIP3 sites under the PTEN arrow.

Why this differs from korkut (where the same bridges moved 156 of 2,715
untouched answers): decryptM's readouts are phospho-sites, and the ones a
drug's own targets leave untouched in a centre are overwhelmingly also not
downstream of the roots the bridges feed; korkut's total-protein readouts
were. The 0-of-2 agreement is consistent with the korkut rating (bridged
calls mostly false) but is no rate: n = 2.

## Verdict

- No arrow gains a rate class; the registered n >= 20 bar is nowhere near met
  (largest possible n was 103 untouched, realized bridged calls 2).
- The p-AKT arrows keep their korkut rate (0.24) untouched.
- The served "unrated" statement is corrected for the three coverable arrows:
  a dataset in the record DOES perturb upstream of them now; its readouts are
  not bridge-reachable. Saying "no dataset perturbs upstream" would be false.
- decryptM remains FIT; this was a fitting-side use (rates are honesty
  metadata, the specs/146 footing).
