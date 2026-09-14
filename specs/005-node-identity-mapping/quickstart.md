# Quickstart: validating the node identity mapping

```bash
export S=/tmp/claude-1001/.../scratchpad
export PATH=/home/awright/julia-1.10.10/bin:$PATH
export PYTHONHASHSEED=0
export LNG_DIAGRAM_DIR=~/reactome-diagrams/97
```

## 1. Reproduce the problem before fixing it

```bash
# 31.3% of the existing context export is orphaned, all of it catalyst+regulator
python3 - <<'PY'
import csv, glob, os, collections
tot=collections.Counter(); orph=collections.Counter()
for d in glob.glob(os.path.expanduser("$S/cat/*/")):
    live=set()
    for r in csv.DictReader(open(d+"logic_network.csv")):
        live.add(r["source_id"]); live.add(r["target_id"])
    for r in csv.DictReader(open(d+"node_reaction_context.csv")):
        tot[r["role"]]+=1
        if r["context_node"] not in live: orph[r["role"]]+=1
print({k:(tot[k],orph[k]) for k in tot})
PY
```

Expected: `input (3641, 0)`, `output (3491, 0)`, `catalyst (2228, 2228)`,
`regulator (1018, 1018)`. If catalyst orphaning is not 100%, LNG #67 has
already been fixed and the baseline needs restating.

## 2. Regenerate with the mapping

```bash
cd ~/gitroot/logic-network-generator
rm -rf $S/cat5/*/cache      # a populated dir makes regeneration a silent no-op
poetry run python bin/create-pathways.py --pathway-list $S/pathways.tsv --output-dir $S/cat5
ls $S/cat5/*/node_resolution.csv | wc -l    # expect 10
```

## 3. The completeness invariant (SC-001, SC-002, SC-006)

```bash
poetry run python scripts/validate_logic_network.py --catalog $S/cat5 --check-resolution
```

Expected: zero unexplained absences in both directions, and
`node_exclusions.csv` **empty** — every currently-unresolved entity is an
EntitySet (29 of 29 in PIP3, 9 of 9 in Cell Cycle Checkpoints), so once set
membership is mapped nothing should remain.

Then the negative control, which is the part that matters:

```bash
poetry run pytest tests/ -k resolution_negative_control
```

A completeness check that cannot fail is worse than none.

## 4. Glyph identity (SC-004, SC-005)

```bash
python3 bench/analysis/glyph_join_report.py --catalog $S/cat5 --diagrams $LNG_DIAGRAM_DIR
```

Baseline before #67 is fixed: 114 of 156 diagram triples join (73.1%), the
42 misses being 25 catalyst, 9 input, 8 output. **Expected after: ~89%**,
because the 25 catalyst misses are the orphaning. The 411 LNG-only triples
are expected — generation descends below the diagram — and should be
reported, not treated as failures.

Check the duplicate case explicitly: ATP has 9 glyphs in R-HSA-1257604 and
each must resolve distinctly by `(reaction, entity, role)`. 0 of 156
triples resolved to more than one glyph, so any ambiguity here is a
regression.

## 5. Set readouts and the denominator (SC-003)

```bash
python3 bench/benchmark_mpbiopath_cases.py ... --catalog $S/cat5 --output-dir $S/m_sets --port 8331
```

Expected: cases discarded for an unresolvable readout fall **204 → ≤27**.

> **The denominator changed.** Roughly 741 scored cases instead of 564. An
> accuracy figure from this arm is not comparable with 365/564, or with any
> number published before this feature. State both denominators or neither.

## 6. Set-combining rules (SC-007)

Run `mean`, `max` and `mean_reachable` as separate arms, one at a time, on
the same catalog, killing the Julia server between them. Report each with
macro-F1, per-pathway net change and the both-arms-converged count —
including the ones that lose.
