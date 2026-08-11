# GSoC 2026 DeltaSignal evaluation

This directory contains the frozen, public-facing report for the ten-pathway
evaluation completed on 11 August 2026.

- DeltaSignal and report code: `59e1757697300af6e9c90e062c0b32195f7204a7`
- Logic Network Generator: `7aca90d03307e53135799a20c2aec75552bec1ba`
- Reactome release: 97
- Primary endpoint: exact key-output entity; reaction proxies disabled
- Solver: current propagation defaults with SCC solving
- Evaluation: 847 experimental cases, three development pathways, seven held
  out pathways, and paired diagram-edge on/off catalogs

Start with [`evaluation_report.md`](evaluation_report.md). The two summary JSON
files retain aggregate and pathway-level metrics. `provenance.json` contains
the sanitized software settings and SHA-256 hashes required to identify every
input and graph artifact; absolute workstation and server paths are excluded.

The TCGA LUAD table in the report is a historical external-association
demonstration. It is not perturbation ground truth and is deliberately kept
separate from the 847-case causal benchmark.
