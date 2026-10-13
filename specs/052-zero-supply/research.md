# specs/052 research: the zero-supply arm

## R1. Result: not adopted (2026-10-13)

**Runs.**
- Arms `zs0` (control, code defaults) and `zs` (`DS_ZERO_SUPPLY_MODE=hold`),
  both by `scripts/run_arm.sh` at **f65618a** on build
  `20261012-1401_01206ca`.
- Results: `builds/20261012-1401_01206ca/results/f65618a/{zs0,zs}`.
- Report: `bench/analysis/zero_supply_report.py BUILD zs0 zs`.
- `zs` held 40,996 node-solves (curator) and 8,076 (experimental). The
  rule ran as registered.

**A void first pair, recorded.** The first pair, at 6273941, changed
predictions (94 curator, 11 experimental). But the solve response did not
return `zero_supply_rule`, so the arm could not show the rule had run. It
was re-run after f65618a added the field, and the bench now refuses an
unknown rule. The predictions of the two `zs` runs were not compared; the
f65618a pair is the record.

| axis | n | control → arm correct | fixed / broken | net | McNemar p | macro-F1 |
|---|---|---|---|---|---|---|
| curator tuning (FIT) | 5,100 | 3,974 → 3,959 | 8 / 23 | **−15** | 0.011 | 0.7645 → 0.7600 |
| curator held-out (SELECT) | 19,000 | 16,848 → 16,828 | 1 / 21 | **−20** | 1.1e-5 | 0.8526 → 0.8512 |
| experimental (FIT) | 849 | 589 → 593 | 6 / 2 | +4 | 0.29 | 0.6159 → 0.6197 |

- By-construction cases (every readout node held by the rule): curator
  +0 / −2, experimental none.
- Knockout→UP subgroup: no wrong UP call was fixed; 10 right ones were
  broken (4 tuning, 5 held-out, 1 experimental).

**Verdict under the registered gate: not adopted.** Curator held-out is
significantly negative. `DS_ZERO_SUPPLY_MODE` stays `off`.

## R2. Traced (as registered: the largest gainer and every pathway losing more than 5)

**HDR, −18 (tuning): an artefact of root pinning, not biology.**
- ERCC1's root forms here are 48 `set_variant` complexes (R-HSA-9980005,
  -7, -8) and one `simple_complex`. The knockout pins them all to 0, and so
  pins every *partner* protein in them to 0 as well.
- The rule cannot tell the knocked-out protein from its partners. It found
  partner leaves with no other source and held 128 downstream nodes.
- 80 of those 128 are produced by reactions with no ERCC1-carrying input
  at all.
- *Corrected (code review):* these counts first read 172 and 120. They were
  computed on the raw network rather than the network as solved, which
  drops closed carrier-loop edges. The mirror now drops them too, and
  matches the solver.
- Restricted to the leaves common to every pin (ERCC1 and ERCC4), it holds
  **0** nodes. The whole HDR loss is the partner-leaf artefact.

**PTK6, −12 (held-out): the hold changes almost-zero to exactly zero.**
- DOK1 knockout pins free DOK1 (a `simple_entity`). The rule holds one node,
  DOK1:PTK6 (R-HSA-8848773), made only from DOK1 and a PTK6 form by
  R-HSA-8848774.
- Solved directly in 052 code:
  - under `off` the complex is already 3.8e-5 (UI scale), a residual of
    the damped iteration;
  - under `hold` it is exactly 0, and the released PTK6 form R-HSA-8848007
    moves from 1.0 to 0.889. The readouts move with it.
- So the loss is not the rule's biology. It is a solver sensitivity: some
  operator, or the pin's effect on how that cyclic component is iterated,
  answers differently at exactly 0 than at 3.8e-7 internal. **Not yet
  isolated.** Candidates: the pin removing the node from its strongly
  connected component (specs/014 knife-edge), or an exact-zero branch in
  the AND or depletion operators.

**RAF, the largest gainer: +5 tuning, +4 experimental, all NF1.**
- NF1 is a RAS GAP. Knocked out, the control reads exactly 1.0. Under the
  rule, the NF1-containing GAP complexes are held at 0, and RAS-GTP
  readouts go UP, as curators and the experiment say.
- This is the open specs/049 item ("an enzyme knockout should raise its
  substrate"). It is the right direction, but the magnitudes rail
  (41–100×).

## R3. decryptM 1 (FIT; it motivated the rule, so it is not evidence)

Served from a pinned worktree at 4bf4d44 (solver identical to f65618a) with
`hold`. Results: `results/decryptm_kinobeads_zs_4bf4d44/`. 230 solves
reported `zero_supply_rule=hold`.

| held-out | kinobeads map, off | hold |
|---|---|---|
| regulated, DS / reachability | 86 / 115 of 278 | 94 / 115 |
| DS − reach, 95% by drug | [−0.139, −0.048] | [−0.108, 0.000] |
| macro-F1, DS / reachability | 0.359 / 0.368 | 0.367 / 0.368 |

**Predictions as registered.**
1. The traced case (gefitinib, EGFR Y1092) reads DOWN: 21.9× → 0. **Met,
   but by construction.** The readout node is an EGFR form, which the rule
   holds.
2. Signaling by EGFR: of the 13 measured-DOWN cases predicted UP, **8 move
   to DOWN and 5 stay UP. Partly met.** 3 measured-unchanged cases also flip
   UP → DOWN.
3. **Met, by construction:** candidates need a pin at exactly 0, so this
   confirms the code rather than the model. Of 12,050 curator and 393
   experimental overexpression cases, 0 changed prediction and 0 changed
   `pred_ui`.

**All 11 changed decryptM cases are EGFR's own autophosphorylation sites
under EGFR inhibitors.** The readout is the drug's target. So the decryptM
gain is by construction, and says nothing about propagation.

It also exposes a data-definition gap in the decryptM bench (specs/021).
Under "drug = knockout", a site on the drug's own target is decided by the
pin, not by the network. Before decryptM is used to fit anything, these
cases should be registered as a separate stratum (`site_gene ∈ targets`).

**"Protein" is loose.** A leaf is any `member_leaves` entry, small
molecules included (there are about 30,000 `R-ALL-*` occurrences in the
canonical node tables). A ligand bound only inside a pinned root complex
is affected the same way as the HDR partners.

## R4. What would be worth trying next (not registered, not run)

1. **Own leaves only.** Candidate leaves are those of zero-pinned *free
   protein* nodes (`simple_entity`), or those common to every zero pin.
   That removes the HDR artefact. It would not touch PTK6.
2. **Isolate the PTK6 sensitivity first.** An exact-zero discontinuity, or
   a pin that restructures a cyclic component, would distort any rule that
   holds nodes at 0. The drug and pathogen rules hold at baseline, so they
   never meet it. Isolating it belongs with specs/014, not here.
3. **NF1 is the real signal.** A GAP or phosphatase knockout raising its
   substrate is the specs/049 open item. A rule aimed at it directly
   (enzymes that convert an active form back) would be narrower than
   zero supply, and could be registered on its own.
