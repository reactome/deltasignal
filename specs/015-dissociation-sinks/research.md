# Where no-path cases die: dissociation sinks, and why bridging them still loses

**Created**: 2026-09-19
**Status**: MECHANISM TRACED to the node; the obvious fix is a door measured shut four times — sizing the one untried variant
**Flag**: none proposed

## 1. How this was found

Tracing `no_path` cases end to end, per Adam's rule (trace cases, not
aggregates). `no_path` is the largest wrong-case bucket on the current dump:
**1,478 of 3,979** wrong cases, every one predicted NORMAL because no route
exists from the perturbed gene to the readout at the uuid level.

Two pathways carry 46%: `Mitotic_G1-G1_S_phases` 430, `Interferon_alpha_beta`
256. In Interferon it is five genes — **JAK1 and IFNAR2 among them, the
pathway's own kinase and receptor** — each failing to reach the same ~25
readouts.

## 2. The trace

For a case, build the network at two resolutions — uuid (what the solver sees)
and stable-id (what the curation says) — find the stable-id path, then walk it
at uuid level and report the first step where the reachable shards of one
entity have no edge into the next.

Mitotic G1, CDKN2A → three different readouts, **identical break**:

```
step 9:  8942622 (1 shard)  -> 68330 (3 shards):   1 edge,  1 reachable
step 10: 68330 (3 shards)   -> 8942803 (30 shards): 12 edges, 0 reachable
!! the reachable shard of 68330 has no edge into 8942803;
   the 12 existing edges leave from 1 OTHER shard
```

`R-HSA-68330` is **CDK4 [cytosol]**, a plain protein. Its three nodes:

| shard | node_kind | in | out |
|---|---|---|---|
| `9d13f1cd` | `dissociation_sink` | 1 (from CDK4/6:CCND:CDKN1A complex) | **0** |
| `e170bccd` | `dissociation_sink` | 1 (from another complex) | **0** |
| `e448aef4` | `simple_entity` | **0** | 23 |

The complex dissociates and drops CDK4 onto a node with no outgoing edges,
while the CDK4 that feeds 23 reactions is a separate node with no incoming
edges. The signal reaches CDK4 and dies. Same shape at CDK6 (76 pairs), TFDP1
(53), IFNAR2 (24), p-STAT2 (21).

Junction tally (one (gene, readout) pair = two cases, both directions):

| pathway | pairs fragmented (stid path exists, uuid path does not) | top 3 junctions cover |
|---|---|---|
| Interferon α/β | 108 of 128 | **all 108** |
| Mitotic G1 | 181 of 215 | 148 |

## 3. Catalog-wide sizing (all 92 pathways, no_path wrong cases = 1,478)

| bucket | cases | share |
|---|---|---|
| **fragmented: stuck on a dissociation sink** | **652** | **44.1%** |
| fragmented: other shard split (frontier: complex 269, entity 94, reaction 2) | 365 | 24.7% |
| disconnected even at stable-id level | 284 | 19.2% |
| readout has no node | 177 | 12.0% |

Sink-stuck by pathway: Mitotic G1 **298**, RUNX1 64, MET 54, SCF-KIT 50,
Interferon α/β 48, AP-2 28. (A first pass mis-mapped Mitotic G1's name and
reported 27.5%; this is the corrected run.)

**The functional node a bridge would feed has fan-out median 112, max 1,976;
≥10 in 646 of 652.** In no case is every functional twin a pure root, so a
"root-only" bridge is not a separable rule. That fan-out number is the whole
story of section 4.

## 4. This is a deliberate design with a measured history — not a bug

`export_nodes` labels `dissociation_sink` post hoc; the construction is in
boundary expansion:

```python
# Fresh per-occurrence readout sink — NOT _leaf_uuid (which would
# share the member's functional node and re-introduce cross-talk).
readout_uuid = str(uuid.uuid4())
```

Adam specified it (2026-05-27): "entities should come out separately so we can
say this part is active, that part isn't." It followed an A/B on the held-out
set:

| boundary edges the solver used | held-out e2e | false-positive change |
|---|---|---|
| assembly + dissociation to **shared** functional nodes | 67.04% | **2,294 (×4)** |
| assembly only / dissociation as **sinks** | **81.98%** | 516 |

Connecting released members to their functional nodes cost **−15pp**. The
sink design was verified at +4.1pp over the pre-boundary baseline.

The solve-time version — bridge the in-1/out-0 receiver to the in-0 feeder —
is `src/core/silo_bridges.jl`, shipped OFF. Measured with a reach cap on 71
pathways, **replicated on an independently regenerated catalog**:

| catalog | correct | net | fixed / broke |
|---|---|---|---|
| 2026-07-16 | 14,960 → 14,887 | −73 | +73 / −146 |
| 2026-09-15 | 18,083 → 18,006 | −77 | +63 / −140 |

126 of 146 losses were NO_CHANGE turning into a change call. Full merge
−204pp; one-bridge-per-silo −4. "Over-coupling, not under-connection, is what
this benchmark punishes." The fan-out of 112 measured above is the mechanism
behind that sentence: the route CDK4 needs is real, but restoring it broadcasts
the complex's activity into ~112 downstream nodes.

**Conclusion: the trace re-derived a known mechanism.** What it adds is the
node-level location and the fan-out number, not a new lever. This door has
been measured shut four times. Do not open it a fifth time with a different
cap.

## 5. A correction to the record

The 2026-09-18 gap anatomy listed "uuid fragmentation, 504 cases,
generation-time fix untried", and it was described in this session as "the
largest untested bucket". That was wrong: the generation-time version is the
May −15pp A/B and the solve-time version is the silo work above. Both tested.

## 6. The one variant that is untried

The composition-gap note (2026-09-17) named it: emit edges along Reactome's
real `hasComponent` hierarchy — cytosolic ISGF3 → ISGF3:KPNA1 →
ISGF3:KPNA1:KPNB1 — rather than flattening a terminal complex to its leaf
subunits. A complex is a component of *few* larger complexes, so this cannot
hub-flood by construction: it is the narrow bridge the silo record asked for.
Whether it covers enough breaks is being sized (does Reactome say `prev` is a
direct component of `next` at each break?). Results: *pending*.

The "other shard split" bucket (365 cases, frontier mostly complex and
entity variants — e.g. `p-STAT1 dimer binds KPNA1` split into variant
reactions whose one-to-one input/output pairing severs set members) is a third
mechanism, unsized and unmeasured.
