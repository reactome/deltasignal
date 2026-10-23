# specs/097: loop audit (structural coverage)

**Status:** analysis, 2026-10-23. Part of the validation plan (specs/096, section 2).

Adam, 2026-10-23: "look through all the rest of the pathways and make some sort of observation
on which loops exist and see if we are handling that type of loop."

## Method (`loop_audit.jl`)

For every network, on the network the solver runs (`solved_carrier_network`): find the cyclic
components; remove the edges a rule handles (pool-internal transitions, specs/039/081/086;
self-contained inhibitor and activator slots, specs/022/040/083/095); the cyclic components left
are **unhandled**, labelled `positive` (gain 1 at baseline, the knife-edge), `self_depleting` (a
complex depleting a node of an entity it contains: a Goodwin ring, specs/066) or `negative`.
Positive ones are sub-labelled by what closes them: `catalyst` (no cycle without catalyst and
regulator edges), `ubiquitin` (the material cycle runs through free ubiquitin) or `material`.
With a case dump, each case is assigned the unhandled loop types on its route.

Two catalogs: the 92 benchmark pathways (canonical `20261022-2011_233665c`) and a sample of 101
other Reactome pathways (`sample_pathways.py` → `sample_pathways.tsv`: every human pathway with
its own diagram and ≥ 5 reactions outside the 92, Disease and Drug ADME, drawn evenly across
top-level areas; build `20261023-1648_c15e27d_sample101`, same generator, no ground truth).

## Results (`bench92_*.tsv`, `sample101_*.tsv`)

| | benchmark (92) | sample (101) |
|---|---|---|
| cyclic components | 213 | 191 |
| unhandled | 97 | 79 |
| nodes in cyclic components | 4.4% | 4.4% |
| nodes in unhandled loops | 2.7% | 3.1% |
| pathways with any unhandled loop | 41 (45%) | 34 (34%) |
| … positive | 30 (33%) | 31 (31%) |
| … closed by catalyst | 19 | 14 |
| … closed by material steps | 7 | **17** |
| … closed by free ubiquitin | 10 | 2 |
| … self-depleting | 21 (23%) | 4 (4%) |
| … other negative | 1 | 2 |

Benchmark cases by route (current defaults, `results/b61b365/dflt095`):

| curator route | cases | accuracy |
|---|---|---|
| no route | 14,652 | 91.4% |
| routed, no unhandled loop | 5,387 | 88.1% |
| through a positive loop | 1,278 | 80.8% |
| through a self-depleting ring | 1,766 | 79.6% |

Experimental: no unhandled loop 83.3% (191), positive 76.4% (55), self-depleting 88.7% (382).

## Reading

- **The benchmark is structurally representative:** the same share of nodes in cyclic
  components (4.4%) and a similar share in unhandled loops (2.7% against 3.1%); positive loops
  in a third of pathways in both. By structure, accuracy should carry over roughly. This is
  coverage, not correctness.
- **Unhandled loops cost about 8 points** on the routes through them (curator 88.1% → 80–81%).
- **Where the sample differs:** fewer self-depleting rings and ubiquitin loops (the benchmark is
  rich in signalling and ubiquitination), and **more material-closed positive loops** (17
  pathways against 7: Macroautophagy, Hedgehog's off state, neurotransmitter release,
  respiratory electron transport, senescence). That class matters more for generalisation than
  the benchmark shows.
- **Against the stop rule (specs/096):** every unhandled type above 5% of pathways needs an
  adopted rule or a recorded refutation: catalyst-closed and material-closed positive loops, and
  self-depleting rings (handled numerically by specs/066's fallback, not structurally). Free
  ubiquitin is under test (specs/098). A constant loop elasticity is refuted (specs/014, 094);
  the knife-edge needs a per-loop or structural answer.
