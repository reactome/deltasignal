# specs/107: rules for material-closed positive loops

**Status:** pre-registered 2026-10-25, before any arm.

Adam, 2026-10-25, after specs/106 and the class table: "okay make the three rules". The decision
recorded in specs/106 applies: these rules are not chosen by tracing experimental cases.

## The loops (specs/106 remnants labelled `positive:material`, traced)

| class | sample (101): cycles / pathways | benchmark (92) | examples |
|---|---|---|---|
| metabolite interconversion | 11 / 9 | 0 | LACT ⇄ PYR (LDH), R5P ⇄ RU5P, GSH ⇄ GSSG, PI ⇄ PI3P, OA ⇄ PYR |
| one protein carried round a multi-step path | 3 / 3 | 19 / 6 | BAP31 (MHC I loading), microtubule, SHC1, TRAF3:MAVS, PAK, tBID:BCL-2 |
| mixed | 6 / 6 | 2 / 2 | GLI:SUFU (Hedgehog), GRIP2 membrane ⇄ cytosol, LC3:PE, SRP, MT3:Zn, cyclin B1 nucleus ⇄ cytoplasm |

Each closes a gain-1 loop: a species reads its own return. Both ideas needed already exist in the
solver: the pool balance (specs/039, a conserved quantity split by the drives of its transitions)
and the carrier rule (specs/049, 077, 088, a species bound and released unchanged reads its supply,
not its own release). The generator recognises only a protein's modification forms as a pool and a
one-step bind/release (E → E:X → E) as a carrier loop.

## Rules

**Rule 1, metabolite pools** (`DS_METABOLITE_POOLS=off|on`, solver). A reaction node with exactly
one non-cofactor input and one non-cofactor output, both small-molecule nodes (`node_kind`
`simple_entity`), is a conversion. Small-molecule nodes in one strongly connected set of the
conversion graph (at least two) form a pool solved by the existing balance: states are the
molecules, each conversion a one-step path (enzyme-driven if it has a catalyst edge). Rationale:
the conserved part (a carbon skeleton, a lipid head group) is redistributed between forms at rates
the enzymes set; its total is set by supply. There is no active form, so every state is a resting
state with an equal baseline share (declared, not fitted). A node already in a generator pool, a
cofactor, drug or pathogen is not used.

**Rule 2, multi-step carriers** (`DS_CARRIER_MULTISTEP=off|strict|wide`, solver; `strict` is rule
2). A free form F (a non-reaction node, not cofactor, drug or pathogen, not managed by a pool) is
carried when a walk from F along material edges (input into a reaction, output out of one, and the
identity edges `variant_pool`, `variant_split`, `specialization`) through nodes that carry F's core
reaches a reaction node R with an `output` edge into F within `CYCLE_MAX_STEPS` (10) reactions.
Under `strict` the core is F's protein leaves (F must have one) and carrying means containing
every core leaf. R is F's release. The (F, R) pairs join the generator's carrier-loop table; a loop
is closed when every activating in-edge of F comes from its releases, and the default carrier rules
(`closed_twin`, `fed_single_ub`) then apply unchanged.

**Rule 3, the mixed loops** (`DS_CARRIER_MULTISTEP=wide`): rule 2 with three relaxations, each
taken from a traced loop:
- a non-cofactor small-molecule F has its small-molecule leaves as core (PE through LC3:PE; Zn2+
  through MT3:7Zn2+);
- leaves compare at the reference-entity level, so a copy in another compartment carries the same
  core (GRIP2 at the membrane and in the cytosol; cyclin B1:CDK1 between nucleus and cytoplasm).
  This needs the generator to emit leaf → reference entity (`leaf_refs.csv`, a generator PR);
- set members: a set-variant node carries F's core if its chosen members do (GLI:SUFU, SRP).

## Evaluation (pre-registered)

No ground truth exists for the sample, so each rule is judged by behaviour on the networks it
touches, against the same networks with the rule off:

- **B1 coverage:** loops (specs/106 audit, extended to the rule) dissolved; reported, not a gate.
- **B2 convergence:** share of solves converged.
- **B3 rails:** share of (perturbation, readout) pairs that read ≥ 99× or ≤ 0.01× under a 2× or
  0.5× root perturbation.
- **B4 dose monotonicity:** share of (root, readout) pairs whose values over root doses
  0, 0.5, 1, 2, 10 are monotone (relative tolerance 1e-6).
- **B5 label invariance:** one perturbation per affected pathway re-solved on a relabelled copy;
  values identical within 1e-9.

Perturbations: every root node from which a node of an affected loop is reachable; readouts: the
affected loop's nodes and the terminal nodes they reach. A rule **passes** if B2 and B4 are not worse
and B3 not higher than with the rule off, and B5 holds. As a no-harm guard on the benchmark,
curator held-out must not lose significantly (McNemar p < 0.05 against); tuning and experimental
are reported, not gates. Adoption is Adam's.

## Predictions

- P1: rule 1 dissolves at least 8 of the 11 sample metabolite cycles (two-substrate exchanges such
  as GOT2's transamination are not one-in/one-out conversions).
- P2: rule 2 dissolves at least 4 of the 6 benchmark pathways' carrier cycles; curator held-out
  moves by fewer than 20 cases either way.
- P3: rules 1–3 together leave material-closed remnants in at most 8 of the 15 sample pathways.

## Amendment 1 (2026-10-25, before any run): rule 3b, `DS_CARRIER_MULTISTEP=compartment`

Rule 3 failed the no-harm guard (held-out −50, research.md), and both traced causes came from the
reference match, not the multi-step walk. Adam chose not to adopt rules 1–3 and to pursue a narrower
rule 3. This variant was chosen after reading curator held-out traces, so its curator figure is read
only as the no-harm guard, not as evidence for the rule.

**Rule 3b** is rule 3 (`wide`) with two changes:

1. Two leaves count as one entity only when they have the same reference entity **and the same
   modified residues**: copies of one protein in different compartments join; SMAD2 and
   p-S465,S467-SMAD2 do not, because a modification changes what the protein does. The generator's
   `leaf_refs.csv` gains a `modifications` column (each leaf's modified residues, sorted; empty when
   it has none). A bundle without the column reports `compartment: no modification column` and adds
   no reference-matched rows (the bench refuses it).
2. Free ubiquitin, the nodes `DS_UBIQUITIN_MODE=inert` would hold (specs/098), is never a carrier
   form. Whether it is held stays `DS_UBIQUITIN_MODE`'s question.

`wide` and `strict` are kept unchanged so the earlier arms reproduce.

**Evaluation**, as registered above, with these choices fixed now:
- Arm `mc107`: `DS_CARRIER_MULTISTEP=compartment` alone, rule 1 off (rule 1 is not adopted, so the
  comparison is against the defaults). Control `ctl107`.
- The behaviour checks reuse the perturbations and readouts already selected (`*_selection.tsv`,
  chosen with every rule on). Rule 3b touches a subset of the nodes `wide` touches, so the
  selection covers it and nothing is re-chosen.
- Pass: B2 and B4 not worse, B3 not higher than rules off, B5 holds, on both catalogs; curator
  held-out not significantly against (McNemar p < 0.05). Tuning and experimental reported.

**Predictions:**
- A1: rule 3b adds no row for SMAD2/3 in TGF-beta and none for free ubiquitin in Mitotic G1 or
  Fanconi; the curator guard passes with held-out net between −10 and +10.
- A2: the behaviour checks pass on both catalogs. Rails fall less than under `wide`: on the
  benchmark to at most 20% (`wide` 9.6%, rules off 32.3%), on the sample to at most 10% (`wide` 3.8%).
  The cyclin B1:CDK1 loop (G2/M), which changes CDK1's modification, is no longer dissolved.
- A3: experimental moves by fewer than 5 cases either way.
