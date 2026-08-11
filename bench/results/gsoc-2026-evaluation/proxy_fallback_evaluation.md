# Reaction-proxy fallback evaluation

## Question

The exact-entity benchmark left 220 of 847 cases unscored. LNG exported an
explicit producing- or consuming-reaction proxy for most missing key-output
entities. This experiment asked whether those proxies can serve as a fallback
readout without silently converting missing predictions to NORMAL.

Exact entities remained preferred. A reaction proxy was used only when the
key-output entity was absent from the generated network.

## Result

| Metric | Exact entities only | Exact plus proxy fallback |
| --- | ---: | ---: |
| Eligible cases | 847 | 847 |
| Scored cases | 627 (74.0%) | 819 (96.7%) |
| Correct cases | 456 | 608 |
| Accuracy among scored cases | 72.7% | 74.2% |
| Correct / all eligible cases | 53.8% | 71.8% |
| Unscored cases | 220 | 28 |

Proxy fallback added 192 predictions, of which 152 were correct (79.2%). On
the three development pathways it added 118 predictions, 102 correct (86.4%).
On the seven held-out pathways it added 74 predictions, 50 correct (67.6%).
Held-out scored accuracy therefore changed from 289/404 (71.5%) to 339/478
(70.9%), while held-out correct coverage increased from 289/503 (57.5%) to
339/503 (67.4%). The main gain is coverage, not a held-out accuracy increase.

On the 819 scored cases, MP-BioPath was correct on 600 (73.3%) and curator
expectations on 645 (78.8%). DeltaSignal was correct on 608 (74.2%). The paired
bootstrap interval for DeltaSignal minus MP-BioPath was -1.8 to +3.8 percentage
points, so the result supports approximate parity rather than superiority.

## Internal proxy-validity control

The ordinary LNG proxy export targets entities missing from the stable-ID
mapping, so it has no direct exact/proxy overlap. For validation, adjacent
reaction proxies were derived from graph structure for exact entities:

- producing proxy: reaction -> entity;
- consuming proxy: entity -> reaction.

Across 606 exact-output cases with an adjacent reaction proxy, entity and proxy
classifications agreed on 592 (97.7%). Exact-entity accuracy was 443/606
(73.1%); proxy accuracy was 436/606 (71.9%). On held-out pathways, agreement
was 379/383 (99.0%), with 276/383 exact and 275/383 proxy predictions correct.
The detailed pathway breakdown is in `proxy_validity.tsv`.

This high agreement supports reaction activity as a practical fallback
readout. It does not prove that every missing entity's proxy is biologically
equivalent: the exact-overlap control and missing-entity cases may differ.
Proxy role and readout semantics should remain explicit in every case record.

## Remaining failures

After proxy fallback, 28 cases remain unscored:

- 25 have gene database identifiers present in Reactome v97 but absent from
  the exported pathway network;
- three key-output identifiers are absent from Reactome v97.

These are mapping/export failures, not model predictions.

## Interpretation

Validated reaction-proxy fallback is the strongest immediate improvement found
in the accuracy sprint. It nearly closes the end-to-end coverage gap while
maintaining comparable scored accuracy. It improves the usable DeltaSignal
system, but it should not be attributed to a new propagation algorithm.

