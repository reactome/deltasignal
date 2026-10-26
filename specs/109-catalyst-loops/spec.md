# specs/109: positive loops closed by a catalyst

**Status:** pre-registered 2026-10-26, before any arm.

Adam, 2026-10-26: rule 3b not adopted; move to the next unhandled loop type, and judge loop rules by
a new gate (below), because a three-class score can hardly credit a rule whose main effect is removing
rails (specs/107, "Why the behaviour gain does not show in accuracy"). Rules are not chosen by tracing
experimental or curator cases (specs/106).

## The loops (specs/106 remnants labelled `positive:catalyst`, traced)

specs/106 counts 20 such loops in 17 benchmark pathways and 35 in 13 sample pathways. Every cycle in them
passes through a catalyst or regulator edge: the node a reaction produces catalyses or regulates
(directly or a few steps on) the same reaction. The small ones, read from the networks with names
from Neo4j:

| shape | examples |
|---|---|
| product catalyses its own formation (autocatalysis) | active caspase-8 catalysing "Caspase-8 processing in the DISC"; MMP2(67-660) catalysing "Autocatalytic activation of proMMP2" |
| activation cascade that loops back | MMP1 activates proMMP2, MMP2 completes MMP1 activation (Activation of MMPs) |
| catalyst also listed as an output | TIMM22 inserting proteins; SCF-beta-TRCP ubiquitinating IkB-alpha; PAK1 binding Rac1-GTP |
| transcription factor switching on its own gene, or two switching each other on | CEBPA; HNF1B ⇄ ONECUT1; LHX1 / OSR1 / PAX2 |
| a complex catalysing the modification of one of its own parts | TGFBR:SARA:SMAD2/3 phosphorylating SMAD2/3; the TNFR1 complex and cIAP ubiquitinating RIPK1 |

Under the AND, the catalyst's fold multiplies the substrate's, so round the loop the gain is exactly
1 and the loop is a knife-edge (specs/014): a small push rails it to 100× or to 0. At a steady state,
autocatalysis cannot make more product than its substrate supplies, and a self-activating gene is
switched on or off by its outside inputs. The model has no saturation for either, so the closing edge
is what rails.

## Rule

`DS_SELF_FED_MODE=catalyst_loop` (solver; `off` stays the default, `entry` stays specs/040 rule A,
refuted as a blanket rule). This is rule A, but applied only where this loop type lives and only to
the edges that close it:

1. **Scope.** A cyclic component the solver iterates qualifies when (a) it holds no inhibitor or
   depletion edge between two of its members (a positive loop), (b) it holds no cycle-balance pool node
   (pools have their own rule), and (c) with its internal catalyst and regulator edges removed, no cycle
   is left (the loop is closed only through catalysts or regulators).
2. **Hold.** In a qualifying component, rule A's self-fed test runs unchanged (a node fed by nothing
   except what it feeds, given this solve's signal-carrying entries), and a self-fed node's **catalyst
   and regulator** out-edges inside the component read the component-entry value, which is baseline.
   Its material edges stay live. A pinned node is never held.

So the loop level follows its outside inputs: active caspase-8 follows procaspase-8 and the DISC, PAK1
follows its own supply, and CEBPA follows RUNX1, PU.1 and the others. A perturbation entering the loop
still runs round it once, up to the point where it would feed itself. The rule is label-free (the
test is a dominator set) and exact at baseline. The solve reports `self_fed_rule`, `self_fed_nodes`,
`self_fed_edges_held` and `self_fed_catalyst_components`.

Things the rule gives up: real positive feedback (a bistable switch that, once on, holds itself on) is
read as following its inputs. That is a statement about the steady state the model computes, not a
claim that such switches do not exist.

## Evaluation (pre-registered)

**Behaviour checks,** with specs/107's harness (`behaviour.jl`, `score_behaviour.py`) and definitions
of B2–B5, on the benchmark (`20261024-2211_7f4b28a`) and the 101-pathway sample
(`20261025-1318_7f4b28a_sample101`). Selection, made once by structure before any solve: per
pathway, the nodes of components that qualify under the scope test (with the rule off, no pins); perturbations are root nodes that
reach them (at most 8, by stable id), readouts are those nodes and up to 10 terminal nodes they reach.

**Benchmark arm** `cat109`: `run_arm.sh cat109 --server DS_SELF_FED_MODE=catalyst_loop` on the canonical
build, against `dflt102` (the defaults on the same build; the commit changes no default). RAF set aside
for the gate, both figures reported.

**Gate for loop rules (Adam, 2026-10-26; applies from this rule on, not to earlier rules):**

- G1. B2–B5 pass on both catalogs: convergence and dose monotonicity not worse, rails not higher than
  with the rule off, relabelling identical within 1e-9.
- G2. Of the scored cases the rule takes off a rail (`pred_ui` ≥ 99 or ≤ 0.01 under `dflt102`, inside under
  `cat109`), at least 90% keep their predicted direction (DOWN / NORMAL / UP).
- G3. Every curator held-out loss is put through specs/108: it is **justified** only if the `dflt102`
  answer is rail-dependent by `rail_audit.py` against `mild108`, and a trace ties the change to an edge
  this rule holds, and the new value fits the stated rationale (a loop following its outside inputs).
  No new artefact class is added.
- G4. Unjustified held-out losses are at most 4, and unjustified losses against gains are not
  significant (McNemar p ≥ 0.05). Experimental not significantly worse (p ≥ 0.05), RAF set aside.
- G5. Raw scores are reported on every axis, beside the justified and unjustified counts.

Adoption remains Adam's.

## Predictions

- C1. The scope test qualifies components in at least 10 of the 17 benchmark pathways and 8 of the 13
  sample pathways that specs/106 lists with catalyst-closed remnants (components differ from the
  remnants: the solver sees pool and self-regulator edges that the audit removed).
- C2. G1 holds on both catalogs, and rails among the selected pairs fall by at least a quarter on each.
- C3. Curator held-out moves by fewer than 15 cases either way; experimental by fewer than 5.
- C4. G2 holds.

## Amendment 1 (2026-10-26, after the first arm, before any re-run): a pair that is also an input is not a catalyst slot

Tracing the large loops (NHEJ, TGF-beta) for the next step found a defect in the scope test and in
specs/106's classification. Reactome often curates one complex as both the input and the catalyst of a
reaction: the complex modifies itself (R-HSA-170868: TGFB1:p-TGFBR:ZFYVE9:SMAD2/3 phosphorylating SMAD2/3;
R-HSA-5693604, R-HSA-5687360 and R-HSA-5693575 in NHEJ). The generator writes both an `input` and a `catalyst`
edge for that pair (7,898 of 22,397 catalyst pairs on the canonical build). `catalyst_edge_pairs` keyed by pair,
so removing "catalyst slots" also removed the input edge, and material cycles were counted as catalyst-closed.
The rule could then hold an input edge.

Fix: a (source, target) pair is a catalyst slot only if every edge between them is `catalyst` or `regulator`
(set-pool `set_member` edges as before). Nothing else changes. Re-run, unchanged in every other respect:
behaviour checks on the same selections (the selection used the same pair test, so it is re-made with the fix
before any solve, and both sets are kept), and arm `cat109b` against `dflt102`, judged by G1–G5.

Predictions: A1. The large components in NHEJ and TGF-beta no longer qualify. A2. The gate still passes, with no
curator class change in either direction, more than 2. A3. Rails among the re-made selection fall by at least a
tenth on both catalogs (the selection now holds only catalyst-closed loops).
