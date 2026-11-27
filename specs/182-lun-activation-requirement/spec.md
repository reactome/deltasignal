# specs/182: does the over-expressed protein need activating?

**Status:** registered 2026-11-27, before any protein was classified.
**Requested by:** Adam, 2026-11-27 ("we need to figure out why it is running so poorly and find solutions"); the next
step named in specs/181.

## Why

specs/181: down a route of one to five reaction steps, the readout moves one time in three. specs/179: when the
protein measurably moves its own direct substrate, the model's calls are right 36% of the time, against 6% when it
does not. The screen over-expresses proteins; Reactome records what active ones do. A measurement of activity is not
available to a user asking an over-expression question, so this asks whether **the network's own structure** says
which proteins act on over-expression alone.

## Classes (from the Release 98 build's tables only; no Lun value read)

For each Lun protein, over every network it is in, take its starting forms (nodes containing the protein that no
non-depletion edge feeds: the forms the bench pins).

- **A, the starting form acts:** a starting form is itself a catalyst or regulator of some reaction.
- **B, binds or modifies itself first:** not A, and a starting form is an input of at least one reaction that has no
  catalyst free of the protein (a binding, or a step the protein catalyses on itself).
- **C, another enzyme first:** not A or B: every reaction a starting form enters is catalysed by something that does
  not contain the protein.
- **N:** no starting form, or the starting forms enter no reaction.

A protein is classified once, pooled over networks (A if A anywhere, else B if B anywhere, else C).

## Measured (Release 98 predictions, both halves pooled and per half; both are FIT)

On pairs with a curated route, self pairs out:
1. Per class: pairs, share of readouts that changed, model change calls and their precision; the same restricted to
   routes of at most 10 reaction steps (specs/181's `routes_181.py` lengths).
2. **Registered comparison:** model call precision in class A minus class C, protein bootstrap (10,000 draws, seed
   182). The hypothesis is supported if the interval is above 0 on the pooled data **and** the point estimate is
   positive on each half. B is reported between them, not judged.
3. Confounds, reported: the same without tyrosine kinases (specs/177's table); the class against specs/179's measured
   engagement for the 85 proteins that have it (is A engaged more often than C?).
4. If supported: the size of withholding class-C change calls under specs/161's rule (fixed vs broken, sign test,
   recall may fall at most 0.01). Sized only; both halves are seen, so nothing is adopted from Lun alone.

If not supported, structure does not tell an active over-expressed protein from an inactive one at this resolution,
and that is recorded.
