# Implementation Plan: Faithful OR semantics for set-valued catalysts

**Spec:** [spec.md](./spec.md) · **Status:** prototyped, all knobs default OFF

## Design

Four independent, env-gated knobs across the two repos.

### LNG — representation
| knob | default | effect |
|---|---|---|
| `LNG_SET_MEMBERS_OR` | `0` | in `append_regulators`, when a **positive** regulator decomposes to >1 member and the entity's labels are `EntitySet`/`DefinedSet`/`CandidateSet`, emit `and_or = "or"` instead of `"and"`. Input/output edges untouched. Participates in the cache fingerprint. |

### DeltaSignal — processing
| knob | default | effect |
|---|---|---|
| `DS_INHIBITOR_OR` | `0` | forward `inhibitor_is_and` into `IndexedReaction` (previously parsed then dropped) and aggregate OR-clustered repressors as the dominant one rather than the product. |
| `DS_OR_MODE=capacity` | `mean` | redundant-alternatives aggregator in **fold space**: `fold = w·mean(foldᵢ) + (1−w)·min(foldᵢ)`. |
| `DS_OR_REDUNDANCY` (`w`) | `1.0` | `1.0` credits full redundancy (≡ `mean`); `0.0` makes any lost alternative decisive. |
| `DS_OR_COMBINE=gate` | `max` | combine an OR cluster with the AND result **multiplicatively** (`and · or_fold`) instead of `max(and, or)`. |

## Why the combination rule was the real blocker

`A = max(and_result, or_result)` treats an OR cluster as an *alternative to the
reaction's own inputs*. Since every reaction has AND-clustered input edges sitting
at baseline, `max` falls back to those and any OR-cluster loss is discarded.

Measured on a 3-node probe (one AND input + a 3-member OR catalyst set, KO one
member):

| combine | w=1.0 | w=0.5 | w=0.0 |
|---|---|---|---|
| `max` | fold 1.0 | fold 1.0 | **fold 1.0** |
| `gate` | fold 0.667 | fold 0.333 | — |

Even a *total* OR-cluster kill (`w=0`) leaves the target unchanged under `max`.
So no aggregator can gate a reaction under `max`; the combination rule had to
change. Logically, `max` is a category error — a catalyst does not substitute for
a substrate; the reaction needs its inputs **and** at least one catalyst.

## Aggregator characterisation (KO 1 member, target fold)

| N | w=1.0 | w=0.75 | w=0.5 | w=0.25 | w=0.0 |
|---|---|---|---|---|---|
| 2 | 0.500 | 0.375 | 0.250 | 0.125 | ~0 |
| 3 | 0.667 | 0.500 | 0.333 | 0.167 | ~0 |
| 8 | 0.875 | 0.656 | 0.438 | 0.219 | ~0 |
| 26 | 0.962 | 0.721 | 0.481 | 0.240 | ~0 |

At `w=1.0`, a KO of 1-of-8 (0.875) and 1-of-26 (0.962) sit **above** the 0.85 DOWN
cutoff — the masking mechanism, quantified.

## Measurements (MP-BioPath experimental, 223 scored cases)

Control = 3-pathway Release97 regeneration with `and` (reproduces the shipped
catalog **exactly**: 167/223, macro-F1 0.7276), so comparisons are clean.

| arm | correct | macro-F1 | Δ ctrl | class-0 recall | McNemar vs ctrl |
|---|---|---|---|---|---|
| control (`and`) | 167 | 0.7276 | — | 0.817 | — |
| OR + `max`, w=1.0 | 116 | 0.5183 | −0.2094 | 0.280 | — |
| OR + `max`, w=0.5 | 119 | 0.5335 | −0.1941 | 0.301 | — |
| OR + **gate**, w=1.0 | 165 | 0.7192 | −0.0084 | 0.817 | p=0.500 |
| OR + **gate**, w=0.75 | 166 | 0.7229 | −0.0047 | 0.828 | p=1.000 |
| OR + **gate**, w=0.5 | 169 | 0.7346 | **+0.0069** | 0.849 | p=0.688 |
| OR + **gate**, w=0.25 | 167 | 0.7256 | −0.0020 | 0.849 | p=1.000 |

**Interpretation.** Under `max` the faithful representation is catastrophic
(class-0 recall 0.817 → 0.280; 50 correctly-predicted DOWN cases become "no
change"). Under gating it reaches **parity at every `w`**, and no arm is
statistically distinguishable from control. The result is therefore *"faithful
representation now costs nothing"*, **not** "gating improves accuracy". The
`w=0.5` bump is noise-level and must not be reported as a win.

## Honesty note on `w`

`w=1.0` is the logically clean redundancy model (lose 1 of N routes → lose 1/N of
throughput). `w<1` has **no mechanistic justification**; it is an empirical
correction for the fact that Reactome sets enumerate every paralog while only a
subset is expressed in a given cell line. The principled long-term fix is
expression-aware set membership (weight alternatives by measured expression),
not a global fudge factor.

## Rejected / superseded

- **"The two bugs partially cancel"** hypothesis (duplicate-edge squaring
  compensating for AND-dilution): **disproved** — combining all three fixes is
  *worse* (110/223, 0.4881), not better.
- Marking input/output set members `or`: wrong; VR splitting already encodes it.
