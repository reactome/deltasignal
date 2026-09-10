# Contract: node resolution

Two consumers, one table. The generator writes
`node_resolution.csv` per pathway; DeltaSignal's benchmark and the pathway
browser read it. Shape is in [data-model.md](../data-model.md).

## Producer obligations (logic-network-generator)

1. Emit one row per `(stable_id, uuid, relation, role, reaction_stid)`.
2. Recurse set membership to leaves; record `depth`.
3. Record `glyph_id` **with** `diagram_stid` wherever the node derives from
   something drawn. Never one without the other — a glyph id is unique only
   within its diagram.
4. Record `release` on every row.
5. Emit `node_exclusions.csv` for entities with no node, `reason`
   non-empty. Expected empty.
6. Never emit a row whose `uuid` is absent from `logic_network.csv`. The
   existing export violates this for 100% of catalyst and regulator rows;
   the validator must enforce it.

## Consumer obligations

1. **Do not use a mapping whose `release` differs from the networks'**
   without an explicit override. Version skew has produced a false finding
   on this project before.
2. Treat a partially resolved set as partial. Do not combine over the
   members that happened to resolve (FR-009).
3. Preserve `relation` when reporting. A `set_member` result and a `self`
   result are different claims and must not be flattened into "the node
   for X".

## Resolution queries

| question | operation |
|---|---|
| what nodes represent entity E? | filter `stable_id == E` |
| what does node U stand for? | filter `uuid == U` |
| what nodes does glyph G mean? | filter `glyph_id == G and diagram_stid == D` |
| which glyphs highlight for node U? | filter `uuid == U`, collect `(diagram_stid, glyph_id)` |
| what is the value of set S? | resolve to `relation == set_member` leaves, aggregate per data-model |

An empty result is a valid answer and must be distinguishable from an
error. "Not represented" is information the interface needs in order to
grey a glyph out rather than appear broken.

## Validator obligations (the completeness invariant)

Runs in both directions and **must be demonstrated to fail**:

1. Forward — every entity participating as input or output of any reaction
   in the pathway appears in `node_resolution.csv` or in
   `node_exclusions.csv` with a reason. No silent absence.
2. Reverse — every uuid in `logic_network.csv` appears in
   `node_resolution.csv`.
3. Referential — every `uuid` in the resolution table exists in
   `logic_network.csv`.
4. Glyph — every `glyph_id` carries a `diagram_stid`, and resolves in that
   diagram.
5. **Negative control** — against a fixture with rows deliberately removed,
   checks 1–4 FAIL. A completeness check that cannot fail is the specific
   defect this project has shipped: `_decomposed_ids` passed 11 of 11 while
   masking 18 dropped catalysts.
