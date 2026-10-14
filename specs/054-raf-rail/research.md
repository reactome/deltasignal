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

## T5. Loop gain at baseline: a bistable switch, not a linear instability

Scripts: `raf_gain.jl`, `raf_sweep_cap.jl`, `raf_thresh.jl`, `raf_sweep_fine.jl`.

**Raw linearisation.** `raf_gain.jl` takes d log F / d log x of every input
at baseline (carrier-loop rule applied, cofactors and pathogens held; no
input has a kink there, and baseline is an exact fixed point). Two cyclic
components lie upstream of the p-MAPK readouts and both have a leading
eigenvalue above 1:
- **194 nodes, 1.32: the RAS GTP/GDP cycle.** RAS:GTP binds a GAP that the
  same cycle recycles, so "RAS GAPs bind RAS:GTP" (R-HSA-5658435) takes the
  loop twice. RAS:GTP does not move in any solve (T1), so the specs/039
  pool solve evidently holds this component; it is not what rails.
- **685 nodes, 1.02: RAF:scaffold:MAP2K:MAPK.** It has two edges with gain 2,
  both Reactome curating one complex as both input and catalyst: "RAF
  phosphorylates MAP2K dimer" (R-HSA-5672978, by R-HSA-5672720) and "MAP2Ks
  phosphorylate MAPKs" (R-HSA-5672973, by R-HSA-5672723). Multiplying
  folds squares such a complex (the specs/012 double count). The drug-binding
  side reactions (R-HSA-9657599, 9657603, 9657606) close further paths
  through the same complexes. The per-edge eigenvalue shares for this
  component do not sum to the eigenvalue (the leading eigenvalue is
  probably degenerate), so they are not used.

**The solve does not show a linear instability.** At gain 1.02 a damped
sweep would grow a deviation about 1.5% per pass. But HRAS at 0.99, 0.999,
1.001 and 1.01 leaves p-MAPK at 1.0 at caps of 500, 2,000 and 8,000
(`raf_sweep_cap.jl`). At 1% the RAF complexes move only about 4e-5 in log
fold (`raf_thresh.jl`).

**It is a switch with a threshold near 3.5%.** `raf_sweep_fine.jl`:

| HRAS | 0.93 | 0.95 | 0.96 | 0.97–1.03 | 1.04 | 1.05 | 1.07 |
|---|---|---|---|---|---|---|---|
| p-MAPK (default) | 51.9 | 51.7 | 51.6 | 1.0 | 51.6 | 51.7 | 51.8 |
| p-MAPK (`DS_DEDUP_ACTIVATORS=1`) | 2.6 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.9 |

Between a 1% and a 10% push, the response of the RAF:scaffold:MAP2K:MAPK
complexes goes from 4e-5 to 4.61 (= log 100, the cap). That ratio is
1e5–1e6 instead of 10. 119 entities move more than 15% at 10%, none at 1%.
Both directions land on the same railed state (T1).

**The double count lowers the threshold but does not make the switch.**
Removing it (`DS_DEDUP_ACTIVATORS=1`) widens the threshold to about 6%.
Pushes beyond that then drift without converging (2.6× at 0.93; 11× at 0.9
and 3.4× at 1.1), and the values depend on the iteration cap. The
eigenvalues stay above 1 (1.018 and 1.138) with no single edge above 1: the
signal reconverges through several in-loop inputs of one AND reaction.
specs/012 measured that deduplication at held-out −15.

**What this changes.** The RAF baseline is the low state of a bistable
switch. Every benchmark perturbation (knockout 0, overexpression 80×)
crosses its threshold, so RAF predictions read the switch's high state
whatever the input. A fix has to remove the high state for a single
perturbation, not just shift the threshold. The candidates are those already
named in "What this means": bring in the curated negative feedback, or
bound the loop's gain (specs/014). Neither is registered.

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
