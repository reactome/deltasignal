# 026 — The composition-hierarchy gap, re-measured under root pinning

## Why (failure anatomy of the canonical baseline, `184dfd5`)

`bench/analysis/case_oracle.py` asks Reactome itself, per failed case, whether
the perturbed gene reaches the readout. Reactome's graph is built with split
set semantics, so there are no sibling hops, and one-level complex assembly.
For held-out cases:

| our network | Reactome's graph | cases | whose fault |
|---|---|---|---|
| no route | no route in the pathway | 535 | the curator reasoned beyond the pathway |
| no route | only via a complex releasing a component | 350 | dissociation handling |
| **no route** | **route, right sign** | **284** | **the generator** |
| route | route | 625 false change, 262 missed, 121 wrong direction | the propagator or curator judgement |

For the 284 severed routes, the first step our network does not reach is:
- **178: component → complex, where that complex does not exist in our
  network.** 150 of these are in Interferon α/β. The complex is a nested
  sub-complex that no reaction uses directly (IFNAR2:p-JAK1:STAT2 exists only
  inside the full IFN receptor complex; ISGF3:KPNA1 only inside
  ISGF3:KPNA1:KPNB1). The generator makes one node per reaction-used complex,
  so the sub-complex never exists, and unless the outer complex is a root
  nothing joins the protein to it.
- 36: component → complex, where the complex node exists but the edge is
  missing.
- 24: entity → reaction, where the node exists but the edge is missing.
- 40: the route starts from a form of the gene other than our pin.

This is the specs/016 composition-hierarchy gap. `LNG_COMPOSITION_EDGES=1`
adds component → containing-complex edges. Under the old broad pins it was a
net loss outside IFN α/β (−60, then −185 in specs/018). But broad pins set
those very complexes directly, which hid the gap, so that measurement does not
carry over.

## Pre-registration (committed before any build or arm)

**Builds** (`scripts/catalog.sh build --variant`; neither moves `current`):
- `ctrl`: no flag. A fresh regeneration, so arms are paired against a catalog
  with the same uuid draw.
- `comp`: `LNG_COMPOSITION_EDGES=1`.

**Arms** (`scripts/run_arm.sh --catalog`, current defaults otherwise):
- `ctrl` on the ctrl build;
- `comp_assembly` (the default `DS_COMPOSITION_MODE=assembly`),
  `comp_limit` and `comp_limit_novel`, each on the comp build.

**Decides:** curator held-out in-release net against `ctrl`, McNemar p,
pathways moved, perturbations. Also reported:
- the experimental axis;
- how many of the 284 severed routes become routed;
- IFN α/β separately, and everything outside it.

**Adopt** a mode only if held-out net > +15 with p < 0.05, and it is **positive
outside IFN α/β** (the specs/016 lesson), and experimental is no worse than
−15. `ctrl` against the canonical `184dfd5` gives the regeneration noise for
this comparison.

## Result: NOT adopted, and the pattern of specs/016 and 018 holds under root pinning

Builds `20260926-0124_d4f4f64_ctrl` and `20260926-0126_d4f4f64_comp` (+5,651
composition edges). Arms at `e3685f5` (`DS_KO_AGG=max`, the same for all).
Regeneration noise, ctrl against the canonical build: held-out −3 (3 / 6).

| mode | held-out net | p | outside IFN α/β | IFN α/β | tuning | experimental |
|---|---|---|---|---|---|---|
| assembly | −973 (385 / 1,358) | 7e-127 | −1,039 | +66 | −361 | −47 |
| limit | −373 (182 / 555) | 1e-44 | −473 | +100 | −203 | −26 |
| limit_novel | −258 (164 / 422) | 4e-27 | −358 | **+100** | −198 | −10 |

- **What it fixes.** Composition edges close the Interferon α/β gap entirely
  (+100, all of that pathway's severed routes).
- **What it costs.** Everywhere else they lose. DNA Double-Strand Break Repair
  goes −201 to −352, and RUNX regulation −139. These are the pathways whose
  loops composition edges re-weld (specs/018).
- **Next step.** The IFN α/β defect is real: a nested sub-complex that no
  reaction uses directly. The fix has to add that one kind of edge without
  closing cycles elsewhere. That needs its own design, for example only where
  the containing complex is a reaction participant and the component has no
  other route to it.

## Corrections after review of PR #74

**The failure anatomy above undercounted generator defects.** Three defects in
`case_oracle.py`, all fixed and tested:
- a set that is a complex component was wired from its "produced" node, not
  its "used" node, which cut PDGFB → "Active PDGF dimers" → receptor complex;
- catalyst- and regulator-only entities were never expanded into components
  or members;
- small-molecule readouts were unreachable by construction.

**The sign label is weaker than stated.** `signed_parities` finds walks, not
paths, so a route through a cycle containing an inhibition yields both
parities. Those cases are now reported as `both_parities`.

Re-run on the canonical `261c94a/baseline` (held-out, 2,471 failures):

| our network | Reactome's graph | cases |
|---|---|---|
| no route | no route in the pathway | **173** (the first version said 535) |
| no route | only via a complex releasing a component | **568** |
| no route | route: sign matches 214, both parities 208, opposite 16 | **438** |
| route | route | 758 false change, 156 missed, 166 wrong direction |

Of the 438 severed routes, the first step our network does not reach is:
- 208: component → complex, where the complex node exists but our component
  has no edge into it (BCL2 → tBID:BCL-2);
- 182: component → complex, where the complex is absent (a nested
  sub-complex; 150 in IFN α/β);
- 36 start from another form of the gene; 6 are entity → reaction.

**Other corrections to the text above:**
- "Close the IFN α/β gap entirely" was an overclaim. Under `limit` and
  `limit_novel`, 75 of IFN's 150 severed routes become correct, and 208 of its
  448 cases stay wrong.
- The regeneration noise quoted as held-out −3 was incomplete. ctrl against
  canonical is also **tuning +104** (all in TP53) and experimental +6, so TP53's
  answers depend on the uuid draw (the specs/013 label-dependence, still
  present).
- The "DSB loops re-welded" mechanism is inferred from specs/018, not measured
  on this build.


## Second correction (third review of PR #74)

**The route check could still hop between set siblings.** It did so when a
complex released a set component, sending the signal to every member. The
release now goes to the set as used. Re-run on `261c94a/baseline`, held-out
no-path misses:

| Reactome's graph | cases |
|---|---|
| no route in the pathway | **339** |
| only via a complex releasing a component | **402** (568 before this fix; that figure was an upper bound) |
| a route our network severs | **438** (unchanged): sign matches 214, both parities 208, opposite 16 |

**First unreached step of the 438**, from the committed
`bench/analysis/route_breaks.py` (the earlier breakdown came from an
uncommitted script and summed to 432):
- 252: component → complex, where the complex is absent from our network;
- 138: component → complex, where the node exists but the edge is missing;
- 36: the route starts from a form of the gene we did not pin;
- 6: entity → reaction;
- 6: the whole route is reached.

**IFN α/β:** it has **200** severed held-out routes, not 150. Under `limit`
and `limit_novel`, **100** of them become correct.
