# specs/140 research: the parity gate

## FR-003: engine pins ≡ bench pins on the frozen catalog

**Setup.** Catalog: the frozen build `20261026-1349_8b6527e`
(`~/deltasignal-catalogs/current`). Bench side:
`bench/analysis/140/dump_pins.py` resolves every benchmark (pathway, gene)
through the bench's own `pin_occurrences` (gene → stable ids via Neo4j, the
same `GENE_NAME_CORRECTIONS`, default env: `DS_PIN_SCOPE=root_cycle`,
`DS_PIN_ROOT_IGNORE=depletion`). Engine side: `bench/analysis/140/pins_140.jl`
resolves the same queries through `parse_complete_network` +
`build_pin_index` + `resolve_pins`, the exact path `/api/pins` serves.
Gate: `bench/analysis/140/compare_140.py` — pins equal as sets, basis and
occurrence counts equal, every difference listed, non-zero exit on any.

**Result (2026-11-06, commit f179409 + this one):**

```
877 (pathway, gene) cases compared; agreeing: 877 (root 806, root_cycle 58, none 13)
no differences — the gate passes
```

82 pathways with a network and a curator file; 877 gene queries. 806 resolve
at true roots, 58 through the specs/024 catalytic-cycle fallback, 13 have no
pinnable form (`none` — returned honestly, never approximated). **The gate
passes with zero differences**, so an API client pinning what `/api/pins`
returns perturbs exactly the nodes every published benchmark figure was
measured at.

## What this does NOT establish

- Readout-side parity. The bench also aggregates readouts
  (`DS_KO_AGG=mean_entity`, specs/105) and reads total-protein antibodies by
  the abundance rule (specs/134). A client reading single node activities is
  not yet reading what the benchmarks scored; that is the next feature of the
  usability programme (spec.md, "The usability programme").
- Dose parity. The benchmarks perturb at 0 (knockout) and 80 (strong
  over-expression, `min(80, 0.8 × ceiling)`); `/api/pins` resolves targets
  only, and a client remains free to pin other doses, where no published
  figure applies.

## Decision

- Decision: `/api/pins` serves the ported protocol; the port is held equal to
  the bench by the committed dump-and-diff harness, to be re-run whenever
  either side's resolution changes.
- Rationale: one protocol, two independent implementations, equality measured
  on all 877 benchmark cases rather than asserted.
- Alternatives considered: resolving pins client-side from the parse response
  (a third implementation to keep in parity, and the parse response lacks
  `member_leaves`/`diagram_entity_id`); having the bench call the API (ties
  every benchmark run to a live server and would have rewired, not verified,
  the protocol).
