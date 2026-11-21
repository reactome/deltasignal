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

## Amendment 2 (2026-11-20): S2 re-run on the rebuilt composites — the verdict stands, now on working bridges

specs/149 amendment 1 found that the composites S2 scored (`20261113-arrows98fix-10c`) dangled 21 of 34 bridges. Two
of the three newly coverable arrows were affected: PTEN Regulation → PTEN → PIP3 and Hedgehog 'off' → GLI3R → RUNX2.
So their "not bridge-reachable" was the defect's zero. Re-run as registered (amendment 1: same case file
`cases_153_centres.tsv`, same `s2_rate.py`, F discipline). Arm `arrows98src_dm153` on
`20261120-arrows98fix-10c-src` (0 of 34 bridges dangling). Control `r98fix_dm153_h` on
`20261113-1311_baf608b_r98fix`, both at `206e3bbb`. The control is identical to `r98fix_dm153` (552ba89) on all 658
cases.

- **The same result**: 658/658 joined; 103 covered control-untouched sites; **2 bridged change calls** (Dasatinib,
  ERBB2 and Mitotic G1, both through p-AKT, both measured "not"). No arrow reaches n ≥ 20, and none gains a rate.
- **PTEN → PIP3**: the bridge is now live. It changes 138 of the 179 decryptM site cases PIP3 already routes, under
  every drug in the table, and F keeps PIP3's own answer there. Of PIP3's 166 untouched site cases, **0 move**,
  including the 25 under AZD8055, Dactolisib and MK2206 (the drugs that reach PTEN Regulation). "Its untouched
  readouts are not bridge-reachable" now holds as a measurement.
- **GLI3R → RUNX2**: 0 of 32 untouched RUNX2 site cases move (8 under Tideglusib). The statement holds.
- **p-S456-ABL1 → RUNX1** was intact before and is unchanged: 0 of 5.
- One routed Mitotic G1 call differs from the old arm (Imatinib, down → not, measured not). F keeps the control's
  answer there, so it does not enter S2.

Served statements: PTEN and GLI3R keep "unrated" with this amendment cited (`src/api/leads.json`, specs/149
amendment 2).
