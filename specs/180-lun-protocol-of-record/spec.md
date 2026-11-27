# specs/180: the Lun 2019 protocol and report of record

**Status:** decided 2026-11-27 (Adam: "ok", on the recommendation after specs/177–179). A decision record plus the
report script; nothing new is measured or chosen.

## Decisions

1. **Protocol of record: adopted.** Lun 2019 is scored as the frozen model under the abundance rule (specs/175) plus
   the tyrosine-kinase rule (specs/177: a tyrosine-specific kinase at a readout with a tyrosine site reads UP where
   the model calls no change), on the Release 98 build (specs/178), with the specs/173 measured call.
   - The rule is a bench-side call on the Lun axis. It is not a solver rule and no other axis uses it.
   - **Every arm is reported side by side**: model + rule, model alone, reachability, rule alone. The gain is the
     rule's (rule alone 0.426 on half B, model + rule 0.420; model alone 0.357, reachability 0.354), so the headline
     is never worded as the model improving.
2. **The tyrosine-kinase statement is not served.** One screen, one cell line, about 80× over-expression, and kinase
   promiscuity cannot be told from phospho-tyrosine antibody cross-reaction. A second dataset, ideally not
   antibody-based, comes first.
3. **Lun is reported by first-hop engagement beside the raw score**, as korkut is by target engagement (specs/138):
   proteins that raise a measured direct substrate, proteins that do not, and proteins with no measured direct pair
   (specs/179). Direct pairs are left out of the split rows. The raw rows stay the headline and gates stay on the raw
   cases; the split never withholds a call (that gate failed on recall, specs/179 E3).

Both halves are FIT. No Lun figure is a held-out test after specs/177's one scoring of half B.

## The report

`bench/analysis/180/lun_report.py PREDICTIONS A_FILE B_FILE [CLASS_DIR] > report.json`. Release 98 predictions
(`20261126-2322_d6d026b_top98`, `results/906698aa/lun178_r98`), class table `bench/analysis/178`. Three-class
balanced accuracy, with change calls and their exact precision.

| pooled, Release 98 | pairs | proteins | model + rule | model | reachability | rule alone |
|---|---|---|---|---|---|---|
| all (raw) | 6,762 | 428 | 0.428 | 0.371 | 0.369 | 0.429 |
| engaged | 762 | 36 | 0.452 | 0.383 | 0.383 | 0.462 |
| not engaged | 1,090 | 49 | 0.328 | 0.328 | 0.307 | 0.332 |
| no measured direct pair | 4,776 | 343 | 0.415 | 0.351 | 0.353 | 0.417 |

| precision of change calls (calls) | model + rule | model | reachability | rule alone |
|---|---|---|---|---|
| all (raw) | 33% (1,235) | 22% (1,012) | 21% (1,487) | 78% (389) |
| engaged | 46% (252) | 36% (195) | 36% (261) | 75% (124) |
| not engaged | 6% (207) | 6% (207) | 5% (289) | 0% (4) |
| no measured direct pair | 34% (642) | 18% (477) | 17% (803) | 81% (211) |

Per half (raw): half A 0.433 / 0.382 / 0.381 / 0.431; half B 0.420 / 0.357 / 0.354 / 0.426.

Reading: for the 49 proteins that did not move their own direct substrate, almost every change call from any arm is
wrong; the model is not above reachability in any row.
