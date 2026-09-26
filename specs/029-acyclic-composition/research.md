# 029 — Composition edges that bridge severed routes without closing cycles

## Why

In the corrected failure anatomy (specs/026 corrections), 438 held-out misses
have a route in Reactome that our network severs. **390 of them break at a
component → complex step**: 208 where the complex node exists but the
component has no edge into it, and 182 where the complex is a nested
sub-complex the generator never creates.

`LNG_COMPOSITION_EDGES=1` adds exactly these edges, but as a class it loses:
−358 to −1,039 held-out outside IFN α/β. The traced cause is that
composition edges re-close the loops that the specs/018 fix opened.

**The rule (`DS_COMPOSITION_FILTER=acyclic`, bench side).** On the composition
build, drop each composition edge that would close a cycle: added one at a
time in sorted order, each checked against the network plus the edges kept so
far. That drops **1,114 of 5,651**; 4,537 remain. Tested:
- a bridging edge is kept and a closing edge is dropped;
- two edges that close a cycle only together are caught.

## Pre-registration (committed before any arm runs)

**Builds:** `20260926-0124_d4f4f64_ctrl` and `20260926-0126_d4f4f64_comp`
(specs/026). All arms are at the same commit, with current defaults
(`DS_KO_AGG=mean`).

| arm | build | settings |
|---|---|---|
| `ctrl` | ctrl | — |
| `comp_acyclic` | comp | `--bench DS_COMPOSITION_FILTER=acyclic` |
| `comp_acyclic_limit_novel` | comp | `--bench DS_COMPOSITION_FILTER=acyclic --server DS_COMPOSITION_MODE=limit_novel` |

- **Decides:** curator held-out in-release net against `ctrl`, McNemar p,
  pathways moved, perturbations.
- **Also reported:** net outside IFN α/β; the experimental axis; how many of
  the 438 severed routes become correct.
- **Adopt** only if held-out net > +15 with p < 0.05, **positive outside IFN
  α/β**, and experimental no worse than −15.
