# specs/046 research

## Implementation status (2026-10-09, logic-network-generator `feat/variant-nodes`)

| Commit | What |
|---|---|
| 0ff17fc | `src/variant_keys.py`: key, choice enumeration (D2: one variable per set per reaction; D3: homo), parser, leaves |
| (next) | `variant_parts`: one format-tolerant parser; every member-reading site routed through it |
| 60147f3 | `src/variant_emission.py` behind `LNG_VARIANT_NODES`; per-variant `containment.csv`; variant-aware boundary layer |
| ba672e2 | `nodes.csv` `member_leaves` = a variant's terminal leaves; emitter unit tests |

The emitter:
- output-only slots are bound by isoform-level reference identity, ties by
  sorted rank and counted, replacing F9's positional tie-break;
- negative regulators are expanded (D4);
- over the cap there is one copy with plain stable ids (the first step of D6).

**FGFR2 (R-HSA-5654738), against the canonical build:**

| | Control | Variant nodes |
|---|---|---|
| Required inputs beyond the curated ones | 1,156 | **0** |
| Cut hand-offs | 544 | **0** |

- All 1,782 variant input edges are fed, with **0** missed joins (I3b).
- The build is identical under `PYTHONHASHSEED` 0 and 7 (10,382 edges).
- Size: 9,912 curated and 470 boundary edges, against 4,880 and 156.

Tests: 366 pass, flag off.

**Not yet in:** the pooled cap fallback (D5/D6) and the heteromer arm (D3).

## Pre-registration: catalog arm (before it runs)

- **Build:** `--variant vn` from `feat/variant-nodes` (ba672e2),
  `--env LNG_VARIANT_NODES=1`.
- **Control:** the canonical build `20261004-1959_6990015`.
- **Scoring:** solver main (`flow` order), code defaults.

**Structural checks** (must hold before any score is read):
1. I3 cut hand-offs = 0 and I3b missed joins = 0, catalog-wide.
2. Coverage: valid cases and the `Pinned:` count within 1% of the control on
   both axes. A drop means genes no longer resolve to nodes; that is a
   generator or benchmark bug to fix, not a result.
3. The build completes for all 92 pathways.

**Predictions:**
- **RAF:** largely unchanged. Its large steps are capped either way, and this
  build uses the single-copy fallback.
- **The pathways with the most cut hand-offs** (FGFR2/3/4/1, R-HSA-194315,
  Mitotic G1, DDX58, EPH-ephrin) change most. Signal that was cut at a
  dissolved complex now propagates.
- **Single-isoform catalyst knockouts move toward no change.** Losing one
  catalyst member now stops only that member's copies. Previously a
  `product` pool stopped every copy.
- **Direction on held-out and experimental is not predicted.** That is the
  honest position: earlier faithful fixes went both ways.

**Decision rule:**
- adopt if held-out ≥ −15 and experimental ≥ −15, and checks 1–3 hold;
- a loss beyond the floor on either axis is traced before anything else;
- report both axes, held-out and tuning, the pathways and genes moved,
  McNemar p, and convergence.

## First catalog arm (2026-10-09): structural checks 1 and 3 pass, check 2 fails; scores not read

Build `20261009-1058_9408b90_vn`: complete (92/92), 193,724 nodes and
365,241 edges.

**Check 1 passes catalog-wide:**
- required inputs beyond the curated ones: 0 (canonical 13,048);
- cut hand-offs: 0 (canonical 2,504);
- 34,378 of 34,378 variant input edges are fed;
- 0 missed joins.

**Check 2 (coverage) fails.** Arm `vn` against canonical:

| | Canonical | Variant |
|---|---|---|
| Curator cases / valid | 24,100 / 23,511 | 22,264 / 19,847 |
| Experimental cases / valid | 849 / 845 | 592 / 385 |

Per the pre-registration, the scores are not read. Two causes:
1. **TP53 skipped.** It has 43,036 edges, over the bench's `--max-edges 40000`.
   Raise the limit for this arm.
2. **Reactions missing from the network.** In PIP3 the emitter, run on the
   pathway's 89 reactions, gives 1,176 copies over all 89. Inside the build it
   gave 551 copies over 66 reactions. `create_pathway_logic_network` passes the
   emitter a reduced `reaction_connections` (reactions with no preceding or
   following event are apparently dropped by then). The ~23 missing reactions
   take their catalysts with them, for example the RTK set R-HSA-2316432 on
   R-HSA-2316434, which is why ERBB2/ERBB3/EGFR/KIT/PDGFRB resolve to nothing.

   **Fix:** build copies from the full reaction list the canonical pipeline
   decomposes, i.e. the reactome_ids in `decomposed_uid_mapping`. Then rebuild,
   and re-run the arm with a larger `--max-edges`.

## Coverage fixed (arm vn4, 2026-10-09)

Three generator fixes, each traced to a gene that stopped resolving:
1. The emitter is given the pathway's full reaction list
   (`pathway_reaction_ids`), not the reduced connectivity table.
2. Set structure is read straight from Neo4j. The primed connector caches
   answered "no members" for entities outside a pathway's prefetch, so a set
   looked empty and its reaction got zero copies (PIP3: 8 of 89 reactions).
   Zero copies is now an error.
3. The boundary layer reads structure through `variant_keys` as well. The
   primed `get_labels` treated the KIT complex R-HSA-205310 as a leaf, which
   lost KIT and PDGFRB in PIP3.

Build `20261009-1531_19f8ad8_vn` (generator 19f8ad8): 92/92 pathways, 226,730
nodes, 458,355 edges. Required inputs beyond the curated ones 0; 37,942 of
37,942 variant inputs fed; 0 missed joins.

| | Canonical | vn4 |
|---|---|---|
| Curator cases / valid | 24,100 / 23,511 | 24,100 / 23,511 |
| Experimental cases / valid | 849 / 845 | 849 / 845 |
| Perturbations (curator / experimental) | 864 / 122 | 864 / 122 |
| Converged (curator) | 1,664 / 1,725 | 1,642 / 1,725 |
| Converged (experimental) | 218 / 244 | 211 / 244 |

**Correction to check 2.** It listed the `Pinned:` count among the
quantities that must stay within 1%. That was wrong. Expansion is meant to
raise it, since a gene now has a root form in every variant: 1,266 → 1,945
(curator) and 198 → 452 (experimental). The check applies to valid cases and
perturbations, which are identical.

## vn4 scores: NOT adopted, loss far beyond the floor

Solver 0525296, code defaults, against canonical (`results/5979e48`):

| | Canonical | vn4 | Net (fixed / broken) | McNemar p |
|---|---|---|---|---|
| Curator held-out | 88.29% / 0.8485 | 87.47% / 0.8345 | −152 (207 / 359) | 1.7e-10 |
| Curator tuning | 76.50% / 0.7567 | 73.24% / 0.7166 | −164 (148 / 312) | 1.6e-14 |
| Experimental | 70.30% / 0.6036 | 63.67% / 0.5540 | −56 (13 / 69) | 2.3e-10 |

Worst pathways (curator): IFN α/β −129, Transcriptional regulation by TP53
−108, RAF −33, PIP3 −28. Best: Intrinsic apoptosis +28.

Under the decision rule the loss is traced before anything else.

### Where the loss is: disconnection, not the representation

Split by whether the readout stays reachable from the perturbed gene:

| | Curator | Experimental |
|---|---|---|
| Reachability unchanged | **+7** | −11 (−9 is PTEN, below) |
| Newly unreachable (`no_path`) | **−361** | −45 |
| Newly reachable | +38 | 0 |

Where the network stays connected, variant nodes are neutral on the curator
axis. Nearly all of the loss is a cut the build introduces. Two mechanisms
are traced to a line.

### Mechanism 1: the cap fallback is a seam (IFN α/β −202 of the no-path loss)

Example traced: JAK1 knockout → readout R-HSA-1015695, "IRF 1-9 [cytosol]",
a DefinedSet of 9.
- **Canonical:** the readout resolves to 504 set-member nodes, and the
  prediction is correct.
- **vn4:** the readout is one plain node, `no_path`.

"Expression of IFN-induced genes" (R-HSA-1015702) is a black box. Its output
sets and its regulator "ISGF3 bound to ISRE" (R-HSA-1015697, two set slots)
put it over the variant cap. So it took the single-copy plain-id fallback
this build implements (D6's last step only).
- **vn4:** the reaction reads the plain node R-HSA-1015697, which nothing
  produces. Upstream produces only variant keys `R-HSA-1015697::variant::…`.
  Each of those has a producer and **0 out-edges**.
- **Canonical:** the one ISGF3:ISRE node has 90 producers and 504 out-edges.

Every readout of that reaction is cut from every gene upstream of it. That is
8 perturbations (JAK1, PTPN11, SOCS1, IFNAR2; knockout and overexpression)
each losing 25 cases.

Catalog-wide there are 72 such plain nodes, which are consumed, never
produced, and have a produced variant twin: PIP3 17, WNT 8, RAF 3, IFN α/β 1,
and others. D6 already prescribed the fix, and this build skipped it:
1. Pool participants whose choice reaches no output. The ISGF3:ISRE
   regulator is one.
2. Make one copy per output variant.
3. Only then fall back to a single copy, reading one pool per varying
   participant.

### Mechanism 2: a regulator takes the first registry node of its entity (PTEN −9 experimental)

Example traced: PTEN knockout in PIP3. Nine readouts go from correct UP
(100×) to 0.

`append_regulators` gives a catalyst or regulator the uuid of the **first**
`entity_uuid_registry` entry with that stable id (`stid_to_existing_uuid`,
first wins). That is registration order, which is arbitrary.
- **Canonical:** first is the PTEN copy produced by R-HSA-8944497, which the
  knockout zeroes.
- **vn4:** the emitter registers reactions in sorted stable-id order, so first
  is an **unfed** PTEN copy, the input of R-HSA-6807106. Its depleters
  collapse under the knockout, so it reads 10× (the de-repression ceiling).
  Through the catalyst/depletion edges, PI(4,5)P2 reads 0.1, and the readout
  reads 4.5e-7.

The defect predates this feature. Canonical has 97 regulator edges on an
unfed copy that has a fed twin; vn4 has 238 (by stable id). The fix is to
prefer a produced node, and among those the one a preceding reaction
produces.

### Not yet traced

- Transcriptional regulation by TP53, −108. AKT1/AKT2 knockout −49 each and
  MDM4 knockout −43 flip correct UP to DOWN. Reachability is unchanged, so
  this is neither mechanism above.
- WNT −58 and RAF −35 of new `no_path`, presumably mechanism 1 (WNT has 8
  seam nodes, RAF 3), but not traced case by case.

### Fixes (generator `feat/variant-nodes` 01c6add + 5514785), before arm vn5

All three are fixes under `LNG_VARIANT_NODES` only; canonical builds are
byte-unchanged.
1. **Cap fallback per D6.** Over the cap:
   - (1) a participant none of whose slots reaches an output is read as
     `<stId>::pool`;
   - (2) otherwise every input-side participant with a slot outside the
     output variants is pooled, giving one copy per output variant;
   - (3) otherwise a single copy reading pools, outputs by plain id (counted).

   A pool is fed by `variant_pool` OR edges from each variant key's produced
   nodes, or its root node if none is produced. That is D5: the mean, the
   same as expanding.
2. **Variant depleters pooled.** Depletion edges reaching one target from
   several variants of one entity become a single edge from a pool of them.
3. **Regulators read the produced copy.** A catalyst or regulator takes:
   - the copy produced by a curated preceding reaction of the regulated
     reaction;
   - else any produced copy;
   - else the old first-registry node.

   For PTEN on R-HSA-199456, Neo4j's precedingEvent includes R-HSA-8944497,
   PTEN translation, which is the copy canonical used.

Pool ids are never looked up in Neo4j. Their `node_resolution` relation is
`variant_pool`, which no readout or gene resolution reads. The solver needs
no change: `variant_pool` is an OR input, `mean`, under every
`DS_SET_POOL_MODE` (test_set_pool, +36 assertions).

**Probe** (IFN α/β, PIP3, TP53; builds complete):

| | vn4 | probe |
|---|---|---|
| Cap seams | 1 / 17 / 0 | 0 / 0 / 0 |
| Duplicated depleter groups | 270 catalog-wide (per-pathway not measured) | 0 in all three |
| Regulators on an unfed copy with a fed twin | 0 / 1 / 0 | 0 / 0 / 0 |
| JAK1 → IRF 1-9 readout nodes reachable | 0 | 9 of 10 |
| PTEN KO → R-HSA-111910 | 4.5e-7 | 100× (canonical 100×; correct UP) |
| TP53 AKT1 KO → R-HSA-5628829 | 0.59× | 1.02× (canonical 1.69×; expected UP) |

The TP53 sign inversion is gone, but the case reads near baseline. The
producer of TP53 (R-HSA-69488) reads 0.97 where vn4 read 2.05, so the
equilibrium of the TP53 loop has moved, not only the depletion term. Left for
the arm, then traced if it still loses.

## Arm vn5 (build `20261009-1854_5514785_vn5`, solver 0525296): NOT adopted

Structural checks pass catalog-wide:
- cap seams 0; duplicated depleters 0; regulators on an unfed copy 0;
- required inputs beyond curated 0; missed joins 0 in all 92 pathways.

D6 resolved 42 capped reaction-builds at step 1 (19 pathways) and 16 at step 3
(9 pathways), summed over the per-pathway log, so a reaction shared by
pathways counts once per pathway.

Coverage is within the 1% check, but not equal. Curator valid cases are
23,367 against 23,511 (−144); experimental is identical. All 144 are in ROBO:
RPL10, RPL22 and RPL5, knockout and overexpression, now resolve to no node.

| | Net (fixed / broken) | McNemar p |
|---|---|---|
| Curator held-out | **+119** (239 / 120) | 3.3e-10 |
| Curator tuning | −52 (217 / 269) | 0.021 |
| Experimental | **−36** (22 / 58) | 7e-5 |

The held-out gain is concentrated: DSB repair +79, intrinsic apoptosis +28,
IFN α/β +25. Experimental is past the floor, so it is traced.

**All −36 experimental is newly unreachable (no path).** RAF −34, WNT −8,
TP53 +9. Two seams were traced to a line, both from this build's own
fallbacks:
1. **Step-3 outputs.**
   - RAF R-HSA-5672980 is one capped copy that writes plain R-HSA-5672701,
     while R-HSA-5674366 reads its variant keys, which nothing feeds. NRAS
     reaches the reaction and stops there.
   - WNT R-HSA-3965447 writes plain Gβγ (R-HSA-167434), and the 180 copies
     of R-HSA-398040 read Gβγ variants that nothing feeds.
   - **Fix:** a capped plain output feeds each existing node of its variant
     keys (`variant_split`, OR), which is what the expanded copies would
     have produced.
2. **Pools of bare sets** (the ROBO coverage loss).
   - A bare set's variant key is its member's own id, not
     `<set>::variant::…`. The pool looked variants up by prefix, found none,
     and the RPL members got no node at all.
   - **Fix:** a pool enumerates the entity's keys. Each key takes its
     produced node, else an existing node, else a new root, as expansion
     would.

Both fixes are in generator commit `feat/variant-nodes` after 5514785. Next
arm: vn6.

## Arm vn6 (build `20261009-1944_715dff6_vn6`, solver 0525296): meets the decision rule

**Structural checks pass.** Cap seams 0, duplicated depleters 0, regulators on
an unfed copy 0, required inputs beyond curated 0, and 0 missed joins in all
92 pathways.

**Coverage is identical to canonical:**

| | Canonical | vn6 |
|---|---|---|
| Valid cases (curator / experimental) | 23,511 / 845 | 23,511 / 845 |
| Perturbations | 864 / 122 | 864 / 122 |
| Converged (curator) | 1,664 / 1,725 | 1,654 / 1,725 |
| Converged (experimental) | 218 / 244 | 215 / 244 |

| | Accuracy / macro-F1 (canonical → vn6) | Net (fixed / broken) | McNemar p |
|---|---|---|---|
| Curator held-out | 88.29% / 0.8485 → 88.96% / 0.8557 | **+124** (233 / 109) | 1.8e-11 |
| Curator tuning | 76.50% / 0.7567 → 76.14% / 0.7486 | −18 (219 / 237) | 0.43 |
| Experimental | 70.30% / 0.6036 → 69.11% / 0.6127 | **−10** (22 / 32) | 0.22 |

The decision rule (held-out ≥ −15, experimental ≥ −15, checks 1–3 hold) is
met.

**Concentration.**

| Pathway | Net | Readouts | Perturbations |
|---|---|---|---|
| DSB repair | +79 | 16 | 36 |
| Intrinsic apoptosis | +28 | 8 | 12 |
| IFN α/β | +25 | 27 | **1** |
| WNT | +24 | 35 | 12 |
| RAF | −17 | 6 | 11 |
| PIP3 | −14 | 14 | **1** |

- Held-out without DSB repair is +45.
- IFN α/β and PIP3 each move through a single perturbation, so their cases
  are correlated and carry no independent weight.
- On the experimental axis, RAF is −11 and TP53 +9.

**Residual cuts**, small but not zero: curator −11 newly unreachable
against +51 newly reachable; experimental −4 newly unreachable. Not traced
yet.

**Size.** Cellular Senescence (R-HSA-2559583, not scored) creates 1,704 root
nodes for 7 pools, which is what expansion would have created.
