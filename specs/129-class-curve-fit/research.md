# specs/129 research

## §0 specs/128 re-read by balanced accuracy (post hoc, 2026-11-03, before this spec was registered)

specs/128 compared arms by cases right. The Lun test (specs/127) and this stage use three-class balanced accuracy, so
the specs/128 arms (solver `110bae4`, frozen build) were re-read that way. Experimental excludes nothing here (845 valid
cases); korkut is its 11,062 scored cases (always-NOT 0.333).

| arm | curator held-out | curator tuning | experimental | korkut |
|---|---|---|---|---|
| `ec_off` | 0.8681 | 0.8681 | 0.6885 | 0.4133 |
| `ec_inh` | 0.8681 | 0.8669 | 0.6925 | 0.4117 |
| `ec_sub90` | 0.8565 | 0.8454 | 0.7020 | 0.4129 |
| `ec_all90` | 0.8524 | 0.8101 | 0.6977 | 0.4132 |
| `ec_all75` | 0.7996 | 0.6801 | 0.6671 | 0.4194 |

The +353 korkut cases of `ec_all75` are +0.006 balanced accuracy: they came from calling fewer changes on a set that is
91% "no change". specs/128's reading 1 ("measured drug data reward attenuation") holds for raw accuracy only; by
balanced accuracy the measured sets favour at most mild contraction (experimental best at `ec_sub90`).
