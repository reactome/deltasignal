# specs/050: topology motifs as versioned rules, re-evaluated each release

**Status:** plan, 2026-10-12. Nothing is registered or measured here yet.
**Origin:** Adam, 2026-10-12, after specs/049. "After we work through the
tests, we go through all the pathways and look for unique topology and do
similar adjustments so it works with other pathways as well. Each release
we could re-evaluate the rules based on newly curated data."

## The claim this is meant to support

Not "we tuned a model to Reactome", but: **recurring topological motifs in
Reactome have a kinetic meaning, and a propagator that reads each motif by
that meaning predicts perturbations better, on pathways it was not tuned on
and on releases that did not exist when the rule was written.**

That makes the motif catalogue a contribution in its own right, to curators
as well as to modelling. It also gives a defence against the forking-paths
critique that any one rule invites (specs/049 took three amendments on one
benchmark).

## Rules already in the register

Each rule reads a motif by structure, never by pathway or protein name.

| motif | rule | spec | status |
|---|---|---|---|
| interconversion cycle (bind, modify, release) | conserved Markov pool, π = πP | 039 | default |
| inhibitor that contains its own reaction's input | shared fold kept at weight 0.1 | 022, 040 B | default |
| closed carrier loop (E → E:X → E, no outside producer) | free form reads its gene's root copy | 049 am. 5–6 | default |
| metabolic cofactor | held at baseline | 007 | default |
| pathogen protein | held at baseline | 048 | default |
| set-valued catalyst or regulator | one OR pool node | 033 | default |
| complex holding sets | one node per consistent combination | 046 | default |
| drug-derived node | held at baseline | 032 | measured, not adopted |
| self-fed loop input | reads its entry value | 040 A | refuted |
| loop as a conserved pool (blanket) | product of entries | 017 | refuted |

Refuted rules stay in the register. Reporting a rule that stopped helping,
or never did, is part of the result.

## Process

**Per rule (as now, made explicit):**
1. Find the motif by tracing failing cases, not from aggregates.
2. Confirm in Neo4j that the shape is curated, not a generator artefact. If
   it is an artefact, fix the generator instead of writing a rule.
3. Give a biochemical argument for what the motif means kinetically, and
   state the assumption it rests on (for carrier loops: bound share
   negligible).
4. Pre-register: rule, arms, predictions, gate.
5. Measure on both axes. Split the results into causal-route cases and
   no-change cases (`route_class.py`), and by pathway.
6. Check against every earlier rule: the canonical bench under all defaults.
7. Adopt or refute, and record either way.

**Motif-level hold-out (new).** A rule is derived from traced cases in a
few pathways. Its evidence of generality is the score on pathways that
contain the motif but were never inspected while the rule was written.
Every spec records which pathways were traced, so the split can be made
afterwards.

**Per Reactome release (new):**
1. Freeze the rule set of release N (rule, defaults, solver commit) in this
   spec before release N+1 is built.
2. Build N+1 with the generator unchanged. Score it with the frozen rules,
   and record that as the prospective result.
3. Run the motif census on N+1. Report new motifs, changes in prevalence,
   and rules whose motif disappeared because curation changed.
4. Only then open new rules for N+1, by the per-rule process.
5. Track the number of rules over time next to accuracy, so growth in the
   rule count is visible.

## The ground-truth gap

A new release brings new curation but no new expected answers. The curator
answers are MP-BioPath's from 2019, on a fixed set of pathways. So a new
release alone tests stability and coverage, not accuracy. Prospective
accuracy needs answers that grow over time. Candidates, to decide early:
- **Curator annotation:** a small set of new perturbation/readout pairs per
  release, written before the model is run on them.
- **The phospho-site empirical axis** (specs/021), which grows with the
  literature, not with us.
- **Published perturbation datasets** (knockout or overexpression with a
  measured readout), mapped to Reactome entities.

Until one exists, the release check is a check of non-regression and
stability, and is reported as that.

## First step

A **motif census** script: per pathway across the catalog, count each
registered motif (carrier loops closed and fed, root twins, interconversion
pools, self-containing inhibitors, drug and pathogen nodes, set pools, and
so on). For each pathway, record whether it was inspected when the rule was
written. Outputs:
- the prevalence table for the catalogue;
- the motif-level hold-out split for every adopted rule;
- the pathways where a motif is common and nothing has been traced, which
  is where to look next.

Analysis only. It changes no prediction.

## Order of work

1. Finish the open specs/049 thread: an enzyme knockout should raise its
   substrate (SMAD7 → phospho-TGFBR1, PTPN11 → SPRY2). It is a motif rule
   and goes through the per-rule process above.
2. Motif census and the motif-level hold-out of the adopted rules.
3. Freeze the rule set for Release97, and prepare the Release98 run.
4. Decide on the ground-truth source for prospective scoring.

## Not in scope

- Rules keyed on pathway or protein names.
- Re-tuning the rules on each release against the same answers. That repeats
  the forking-paths problem every release, which this plan exists to avoid.

## Release97 freeze (2026-10-12)

The rule set scored prospectively on Release98. It must not change before
Release98 is built and scored under it.

| | |
|---|---|
| solver | deltasignal `main` at **8757909** (DS #103). `src/`, `cli/` and the bench are identical to `1f58a1b`, which produced the results below |
| generator | logic-network-generator `main` at **01206ca** (LNG #111), default environment |
| catalog | `20261012-1401_01206ca`, Reactome Release97 (2026-06-17), 92 pathways |
| results | `builds/20261012-1401_01206ca/results/1f58a1b` |
| protocol | `DS_PIN_SCOPE=root_cycle`, UP pin 80, knockout pin 0, call thresholds of `bench/benchmark_vs_mpbiopath.py` at 8757909 |
| scores | curator held-out in release 89.29% / 0.8612; curator all 86.40% / 0.8316; experimental 589/849 (69.38%) / 0.6159 |

**Rules in force** (code defaults at 8757909; the register above has the
specs):
- cycle balance (039);
- self-inhibitor weight 0.1 with leaf pairs (022, 040 B);
- `closed_twin` carrier loops (049);
- cofactors inert (007);
- pathogen proteins inert (048);
- set pools in product mode (033);
- variant complex nodes (046, generator);
- depletion bounds 0.1–10 (011);
- `hill_sat` AND, `divide` inhibition, OR mean (009, 010).

**Release98 procedure.**
1. Build Release98 with the generator at 01206ca and the default
   environment, and benchmark it with the solver at 8757909. Record that
   as the prospective result, on whatever ground truth exists for
   Release98 (see "The ground-truth gap").
2. Run `motif_census.py` on the Release98 build and diff it against R2.
3. Only then may generator or solver changes made after the freeze be
   scored on Release98, each as its own arm against step 1.

Rules adopted after this freeze are scored against Release97 as usual. They
enter the prospective record only at the next freeze.
