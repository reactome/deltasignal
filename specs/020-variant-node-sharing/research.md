# Re-measure: variant-node sharing on the current catalog

`LNG_SHARE_VARIANT_NODES` (introduced in specs/016, LNG #90, off by default)
keys Phase-1 UUIDs on (entity, **Reactome** reaction, role) instead of
(entity, **virtual** reaction, role), so a participant that does not vary
across a reaction's set-member variants stops being photocopied once per
variant.


> **Reproduction note (2026-09-21).** The `LNG_SHARE_VARIANT_NODES` and
> `LNG_BOUNDARY_LEAF_REUSE` settings in the commands below **no longer exist**:
> variant sharing and downstream-free leaf reuse are unconditional, and setting
> either name is now a hard error (LNG #93, #94). That is deliberate — a stale
> value must not let a run quietly measure the default. It does mean these
> commands will abort rather than reproduce; the recorded numbers are the
> evidence. To re-measure the *old* behaviour you have to check out the
> generator commit named in the arm, not set a flag.

## Why re-measure

The only measurement (specs/016) ran on `cat_fix`, the **v1** loop-fix catalog
that specs/018 superseded, and the adversarial review then showed its one
visible gain (+18 tuning) was 12 cases of cross-catalog relabelling churn plus
6. The decision to flip the default needs a clean measurement on the current
generator with a same-catalog control.

## What the change does to the graph (measured on HDR, R-HSA-5693567)

BCDX2 (R-HSA-5685316) appears in six (reaction, role) contexts. Four of them
are outputs of reactions that expand into 8 variants each, so BCDX2 was
photocopied 8 times per reaction: **33 nodes → 5**. The other participants of
those reactions are photocopied identically (each of R-HSA-5686410's four
feeders also appears ×8) — the whole reaction is duplicated when only one
participant differs between the copies.

Wiring, before and after:

| | before | after |
|---|---|---|
| BCDX2 as output of R-HSA-5686410 | 8 nodes, each with 1 producer, each feeding 4 dissociation sinks | 1 node with **8 producers, OR-clustered**, feeding 4 sinks |
| BCDX2 as input to R-HSA-5685341 | 1 node | 1 node (unchanged) |

### Semantics, case by case

1. **Variants uniform** (the common case): before, 8 nodes each at v; after, one node at mean(v…v) = v. **Identical.**
2. **Perturbation upstream of all variants**: all 8 at the same value → mean is that value. **Identical.**
3. **Perturbation hits one variant's distinguishing set member**: before, one copy at 0 and seven at 1, and readout resolution takes the **max** over copies, so the effect is **masked**; after, one node at mean = 7/8, so a small effect survives. **This is the one genuine semantic change** — max-over-copies becomes mean-over-producers.
4. **AND multiplication of copies** (the dangerous case): does not arise here — the consumer side was already a single node; the copies only fed dissociation sinks. Confirmed by inspection, not assumed.

## Pre-registration (committed before the regeneration)

Regenerate `cat_share2` from the current generator with
`LNG_SHARE_VARIANT_NODES=1`; control is the current production catalog
`cat_fix2` from the same generator without it. Cross-catalog by necessity
(uuids differ), so **both a same-catalog reference and the churn bound apply**:
specs/019 C3 measured cross-catalog churn at up to 96 cases, all in TP53, zero
held-out.

- **P1 (size)**: nodes fall ≥ 25%; oracle root→terminal reachability
  **unchanged in every pathway** (it was identical in 92/92 on the old
  catalog); cyclic component census unchanged.
- **P2 (accuracy)**: held-out within ±20 and not significant (p > 0.05) — the
  claim is neutrality, and a significant move either way falsifies it.
- **P3 (the case-3 risk)**: false change does not rise by more than 20; any
  rise is concentrated in pathways with high variant multiplicity.
- **P4 (experimental)**: conditioned net ≥ −5.
- **P5 (churn discipline)**: TP53 tuning figures are **not** quoted as an
  effect; only held-out and non-TP53 tuning are.
- **Decision**: flip the default ON if P1 holds and P2–P4 hold — i.e. adopt on
  size and faithfulness grounds given accuracy neutrality. If held-out moves
  significantly in either direction, the change is no longer a pure
  correctness fix and gets argued on its merits.

---

## Corrections to the sections above (adversarial review, 2026-09-21)

Three claims in "What the change does to the graph" were wrong. They are left
in place above and corrected here rather than edited away, because the wrong
version is what the pre-registration was written against.

1. **"Each of R-HSA-5686410's four feeders also appears ×8" is false — they
   appear ×1.** The 8 variant reactions share 3 input nodes (24 input edges
   from 3 distinct sources) and 3 catalyst/regulator nodes, *identically in
   both catalogs*. Nothing on the input side changes under sharing at all.
   Only the output side is photocopied, at multiplicities 8/4/2/1.
2. **Phase 1's photocopying is mostly already undone before emission.** Phase 2
   unions `(eid, producer_vr, "output")` with `(eid, consumer_vr, "input")`
   across all variant pairs of connected reactions, and the boundary caches
   share root inputs and terminal outputs per stable id. Catalog-wide only
   **1,431 of 29,655 (4.8%)** `(entity, reaction, role)` contexts end up with
   more than one node. Catalyst and regulator contexts have **zero** multi-node
   groups, which is why handling only input and output in `_register_phase1` is
   not a gap.
3. **Case 4's conclusion holds but its stated reason was wrong.** The collapsed
   copies do *not* "only feed dissociation sinks" — they emit 26,671
   dissociation, 7,068 input, 1,212 assembly, 118 catalyst, 35 regulator and 8
   depletion edges, nearly all AND-typed, and 821 of the 1,431 groups are
   role=input feeding reaction nodes through AND. The correct argument is
   structural: the Phase-1 key already contains the variant reaction id, so no
   single variant reaction ever consumed two copies of the same entity in the
   same role. Merging can therefore only raise a node's *fan-out*, which does
   not multiply in this solver. Verified exhaustively: **0 of 1,431 collapsed
   groups gain an AND input.**

**Case weights were wrong too.** Case 3 is not "the one genuine semantic
change":

| | groups | |
|---|---|---|
| role=input, identical producer sets | 821 | value-preserving |
| role=output, producing variants have identical inputs | 365 | value-preserving |
| role=output, producing variants **differ** | **245** | genuine max→mean |

So 1,186 of 1,431 (82.9%) are exactly value-preserving and 245 (17.1%) are a
real change — 40% of all output-role groups, an order of magnitude more than
"the one". Two further omissions: the merged value **propagates downstream** to
every consumer, not just to the readout; and the max→mean only reaches the
classifier when the group is small. With `DOWN_CUTOFF = 0.85`, mean(0, (N−1)
baselines) crosses the cutoff only for **N ≤ 6** (840 of the 1,431 groups). For
the flagship 8-variant example it changes nothing. The same arithmetic
manufactures **false DOWNs** where the knocked-out member is irrelevant, which
is exactly what P3 exists to bound.

**A fifth case the write-up missed.** Sharing merges within one `(entity,
reaction, role)` context, but the readout aggregation takes the max over *all*
of a stable id's uuids. In HDR, 101 of 259 stable ids still map to more than
one uuid after sharing (down from 110, max 65 → 28). The unmasking case 3
describes therefore requires every surviving context to move together.

## Results

Regenerated `cat_share2` at generator `13847fa` against control `cat_fix2` at
`79feca7`. Those two trees are **byte-identical** (`13847fa` is the merge commit
whose second parent is `79feca7`; `git diff` between them is empty), same
pathway list, same `PYTHONHASHSEED=0`, both `dirty=0`, both regenerated from
scratch — so the generator is not a confound. The arms were scored under
different deltasignal SHAs (`560a85c` for the control's headline, `72015ce` for
the sharing arm); this is checked, not assumed — `sb2_ctrl` at `72015ce`
reproduces the control's 20414/24100 and macro-F1 0.8137 exactly.

### P1 (size) — PASS, but the headline overstates it

| | cat_fix2 | cat_share2 | |
|---|---|---|---|
| nodes | 108,229 | 70,738 | **−34.6%** |
| edges | 285,497 | 224,189 | **−21.5%** |
| cyclic components | 210 | 204 | |
| nodes in cycles | 8,407 | 8,131 | |
| largest component | 675 | 675 | unchanged |

The ≥25% bar passes on nodes. Three caveats keep it honest. **Sharing is a
literal no-op in 47 of 92 pathways** (node count bit-identical), and the median
reduction among the 45 that do move is only **6.7%** — the catalog figure is
carried by a few pathways (Fanconi Anaemia −83.5%, R-HSA-5693532 −74.2%, TP53
−61.8%). About a third of all catalog nodes are dissociation sinks with no
outgoing edges, so node count flatters the structural change; **the −21.5% edge
figure is the more conservative read.** And the ≥25% threshold should not have
been wired into the decision rule for what is argued as a correctness fix.

**The new-cycle risk is checked and does not materialise.** Merging two
variant strands could in principle weld a cycle that neither had. Measured
across all 92 pathways: **0 pathways gained cyclic nodes.** Four pathways lost
a component. This matters because unwelding loops was the largest accuracy win
this project has had (specs/018), and sharing does not undo it.

### P2 (accuracy) — PASS, and the null is unusually strong

| | cat_fix2 | cat_share2 |
|---|---|---|
| curator accuracy | 20414/24100 = 84.71% | 20400/24100 = 84.65% |
| curator macro-F1 | 0.8137 | 0.8125 |

Raw prediction disagreement across the full 24,100 shared cases is **29 cases,
0.12%, in 3 pathways** — TP53 14, RUNX2 10, MET 5. Only TP53 has a nonzero net
(−14); RUNX2 and MET fixed exactly as many as they broke.

Under perturbation-set conditioning, **held-out net is 0 with zero discordant
cases**, and so is tuning. Every discordant case was dropped as non-comparable,
because sharing changes which nodes a gene resolves to (4,094 of 24,100 cases
resolve a different perturbation set).

This deserves a caveat the pre-registration did not carry: **P2's "not
significant" clause was vacuous as written.** At the usual n_discordant ≈ 349,
significance needs |net| ≳ 37, so any |net| ≤ 20 is automatically
non-significant and the conjunction collapses to the ±20 bound; power against a
true net of 20 is about 0.19. P2 was an equivalence claim tested with a
null-hypothesis test, and future pre-registrations should state it as an
equivalence bound. In *this* case the objection is moot for a reason that makes
the result stronger, not weaker: n_discordant is not small, it is **zero**.
There is no movement to be underpowered against.

### P3 (false change) — PASS

1,150 → 1,150, exactly unchanged, loop-heavy 379 → 379 and other 771 → 771. The
false-DOWN mechanism identified above is real but does not fire at catalog
scale.

### P4 (experimental) — PASS

611/849 = 71.97% → 609/849 = 71.73% unconditioned (−2); macro-F1 0.6451 →
0.6432. Conditioned: no comparable case moved, same as the curator axis.

### P5 (churn discipline) — complied with

TP53's −14 is a tuning figure and is not quoted as an effect. Note it is also
not attributable to relabelling churn in the usual way: TP53 is one of the most
heavily collapsed pathways (−61.8% nodes), so a real semantic change there is
expected. The conditioning says those 14 cases are not like-for-like.

### Verdict

**All of P1–P4 hold, so the pre-registered decision is to adopt.** The honest
one-line summary is: *variant sharing removes a third of the nodes and a fifth
of the edges, creates no cycles, and changes 29 of 24,100 predictions with zero
net effect on held-out accuracy on either ground-truth axis.*

It is adopted as a **faithfulness and size** change, not an accuracy change.
Photocopying a participant once per set-member combination was never intended;
one node per (entity, Reactome reaction, role) is what the data says. The
measurement's job was to show the correction costs nothing, and it does not.

## Filed separately: a larger redundancy that sharing does NOT fix

While checking the above, the review found that **exact structural twins
dominate the catalog**. Counting nodes whose incoming edge set *and* outgoing
edge set are byte-identical to a sibling's:

| | nodes with in-edges | redundant twins | |
|---|---|---|---|
| cat_fix2 | 92,900 | 55,608 | 59.9% |
| cat_share2 | 59,276 | 34,371 | 58.0% |

Sharing removes 38% of the redundant twins but leaves the *rate* untouched,
because it does not merge reaction nodes: variant reactions whose
distinguishing members collapse onto the same emitted nodes stay as exact
duplicates. At the reaction level that is 1,074 of 5,030 reaction instances
(21.4%) and 21,357 of 39,913 variant-reaction nodes (53.5%), with eleven
reactions emitting 419 redundant variants out of 420.

**Do not read this as an available win.** The closest intervention already
measured — solver-side activator deduplication — is **load-bearing and
negative** (specs/012: −15 held-out, p = 0.0015), and self-contained inhibitor
removal likewise (−61, p < 0.0001). On these networks, bounding a runaway
operator has paid off and deleting a structurally-duplicated term has not,
twice. A generator-side twin collapse is a different intervention from either,
but the prior is against it. Worth its own spec and its own pre-registration,
not a patch.

## Why R-HSA-5686410 produces 8 variants and not 6

Question resolved. The reaction has two EntitySet participants, both
outputs: R-HSA-5685980 with 3 `hasMember`, and R-HSA-9707299 with 2. The naive
product is 3 × 2 = 6.

The first member of R-HSA-5685980, **R-HSA-5685976 "EXO1:BLM,WRN", is itself a
Complex containing an EntitySet** (BLM,WRN). `_complex_variant_leafsets` and
`_matching_leaves` recurse into exactly that case, so it contributes **2**
leaf-variants rather than 1. R-HSA-5685980 therefore offers 2+1+1 = **4**
alternatives, and 4 × 2 = **8**. Confirmed against the emitted graph: the 8
variants' distinguishing outputs are exactly
`{EXO1+BLM | EXO1+WRN | DNA2:WRN | DNA2:BLM:TOP3A:RMI1:RMI2} × {p-BRCA1:BARD1 |
K6PolyUb-p-BRCA1:BARD1}`.

**This is correct behaviour, not a second bug.** The right denominator for any
future variant-count audit is the product of each EntitySet's *recursive
leaf-variant* count, not its direct `hasMember` count. Against the naive
denominator 1,409 of 4,211 reactions (33.5%) appear to "exceed" their variant
count, median ratio 4.0 — an artefact of the wrong denominator.

The generator docstring saying "33 copies from 5 reactions" should say **6**
reactions (5685318, 5685341, 5686410, 5686469, 5686483, 5693589); four expand
×8, giving 34 context rows on 33 distinct nodes.
