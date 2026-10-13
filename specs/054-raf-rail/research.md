# specs/054 research: the RAF readouts rail whichever way they are pushed

**Status:** trace only (2026-10-13). Nothing registered.

Build `20261012-1401_01206ca`, solver at a807851 (= main), canonical case
log `results/a807851/cj0`. Scripts: `bench/analysis/054/`.

## T1. There is no inversion: every perturbation rails p-MAPK

specs/053 R1 read NF1 overexpression as "RAS:GTP down, p-MAPK up", an
inversion. On the canonical build the RAF picture is simpler and worse:
- Every RAS/RAF knockout (ARAF, BRAF, RAF1, HRAS, KRAS, NRAS) reads the
  p-MAPK readouts UP at 30–57× (IL17RD:p-MAP2K:p-MAPK at 100×). Wrong in 13
  experimental and about 30 curator cases.
- The overexpressions read UP at the same magnitudes, and are right.

HRAS swept on its pinned root (`raf_sweep.jl`), p-MAPK dimers (1268261):

| HRAS | 0 | 0.5 | 0.9 | 0.99 | 1.0 | 1.01 | 1.1 | 2 | 80 |
|---|---|---|---|---|---|---|---|---|---|
| p-MAPK | 53.8 | 54.1 | 52.1 | 1.0 | 1.0 | 1.0 | 52.0 | 49.5 | 50.3 |
| converged | no | no | no | no | **yes** | no | no | no | no |

The baseline is an unstable state: a 10% push in either direction lands on
the same railed state. So the overexpression cases are right for the wrong
reason, and the knockout cases wrong for the same one.

RAS:GTP (R-HSA-109783) does not respond to HRAS at all (0.9986 under the
knockout, 0.985 at 80×), so the rail is not carried by RAS:GTP. Not traced.

## T2. The iteration cap does not decide it

`raf_cap.jl`: the readout is 53.8 (knockout) and 50.28 (80×) at every cap
from 250 to 5,000; only cap 100 differs. specs/049 P4 stands. The solves
report non-converged at every cap, so some other component never settles;
it does not move these readouts.

## T3. The drug-complex regulators decide between rail and collapse

`raf_drugreg.jl` removes only the 536 negative regulator edges whose source
is a drug-derived node (288 nodes; drug complexes stay as participants):
- every RAS/RAF perturbation, in both directions, then reads p-MAPK at
  ~1e-6 (collapse), and RAS:RAF complexes at 0.3–0.75;
- 20 curator / 10 experimental cases flip right (the knockouts) and 10 / 10
  flip wrong (the overexpressions).

So neither state is signal. Holding the drug complexes (specs/032) or
removing their regulation only chooses which side of the unstable baseline
the loop falls to. This is why `kcd` fixed knockouts "by collapse"
(specs/049 P8).

## T4. The loop is curated

Neo4j (Release 97): "Phosphorylation of RAF" (R-HSA-5672969) is catalysed
by the set "RAF activating kinases" (R-HSA-5672734), whose members are
MyrG-p-Y419-SRC, **p-2S MAP2K homo/heterodimers**, **p-BRAF, p-ARAF,
p-RAF1**, p-T286-CaMKII:CALM1 and JAK2. Phosphorylated MAP2K and RAF feed
back onto RAF phosphorylation: positive feedback, as curated. All 320 copies
of R-HSA-5672969 in the network take that set as their catalyst. Reactome
curates the pathway's negative feedback (ERK on RAF and SOS) elsewhere, so it
is not in this network to balance the loop. MP-BioPath's networks were
hand-edited to remove such loops (specs/034).

`raf_pool.jl`: `DS_SET_POOL_MODE` = product, geomean, mean or extreme give
identical readouts. The set-pool combination is not the lever.

## What this means

- The RAF experimental and curator scores (28/49 experimental) are a coin
  the unstable baseline flips. Any rule scored on RAF is scored on that
  coin, which explains RAF's history across specs/032, 045, 048, 049, 052
  and 053.
- The lever is the loop's gain at baseline, the specs/014 knife-edge on a
  curated positive feedback, not a drug or pool rule.
- Next, if pursued: measure the loop gain at baseline for the RAF component
  (which edges carry gain > 1), and decide what the model should do with a
  curated positive feedback whose balancing negative feedback lives in
  another pathway. That is a design question for Adam before any rule.
