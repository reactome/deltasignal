# Contract: cycle structure report

`bench/analysis/cycle_structure.py` is the reproducible measurement required
by FR-010. It is a CLI tool, so the contract is its arguments and output.

## Input

| argument | required | meaning |
|---|---|---|
| `--catalog PATH` | one of these two | a directory of LNG pathway directories, each holding `logic_network.csv` |
| `--mpbiopath PATH` | one of these two | a directory of MP-BioPath four-column TSVs |
| `--classify` | no | add the artifact/feedback classification |
| `--ratio-threshold FLOAT` | no, default 15 | nodes-per-reaction cut; echoed in the output |
| `--json PATH` | no | machine-readable output in addition to the table |

Exactly one of `--catalog` / `--mpbiopath` must be given; supplying both or
neither is an error, not a default.

## Output

Per pathway: node count, edge count, self-loop count, number of components
larger than one node, total cycle-resident nodes, largest component size.
With `--classify`, additionally per component: intra-edge count, negative
intra-edge count, distinct reaction stable ids, nodes-per-reaction, and the
assigned class.

A totals row is always emitted. The threshold in force is always printed,
including when it is the default — a classification whose parameter is
implicit is not reproducible.

## Guarantees

- Reads only; never writes into the catalog.
- Deterministic: component identity depends on the graph, not on iteration
  order, and node sets are reported sorted.
- Exits non-zero on a malformed network rather than reporting partial
  results as complete.
