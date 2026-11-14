# specs/155 tasks

- [x] T001 Research: pair sources (R1), decryptM numeric availability (R1a: the
      paper's fitted effect sizes — admitted), stratum reuse (R2), table/lookup
      pattern (R3). — research.md
- [x] T002 P1: `bench/analysis/155/fit_155.py` — pairs, LODO, bootstrap, gates;
      the circular tied-fold band bug caught in-sample and fixed before any
      gate reading. First (mixed-population) measurement recorded. — e5632ab
- [x] T003 P2: loader with provenance refusal, stratum-walk lookup, arm-(b)
      band rule (src/api/readouts.jl); handler fields and top-level status
      (src/api/server.jl); docs/API.md. — 06911bd
- [x] T004 P3 tests: test_fold_ranges.jl (233 asserts; both mutations proven
      caught), registered bare in asserting_suites.txt and the CLAUDE.md
      table. — 06911bd, 06fe83a
- [x] T005 P3 parity: the whole 98 record (24,356 folds) replayed through the
      served stack with the table serving — max diff 5e-7; two replay traps
      recorded. — 10c831f
- [x] T006 Amendment 1 (pre-serving): gates re-measured on the SERVED
      population (untouched pairs excluded) — korkut (a), decryptM (b);
      table, docs and suite updated. — 06fe83a
- [ ] T007 P4: PR, adversarial review + green CI, merge; dev API recreated and
      live-smoked on the merged main.
- [ ] T008 WebsiteAngular consumption (its own spec, the specs/154 → WA
      specs/026 pattern).
