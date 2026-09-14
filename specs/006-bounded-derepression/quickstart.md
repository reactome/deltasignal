# Quickstart: bounded de-repression

```bash
export S=/tmp/claude-1001/.../scratchpad
export W=~/codes_and_results
export PATH=/home/awright/julia-1.10.10/bin:$PATH
```

## 1. See the defect before changing anything

```bash
julia --project=. -e '
  bl = 0.01
  for eps in (1e-9, 1e-4, 1e-3, 1e-2)
    println("eps=", eps,
            "  max de-repression=", (bl+eps)/eps,
            "  at x=bl/2=", (bl+eps)/(bl/2+eps),
            "  max suppression=", (bl+eps)/(1+eps))
  end'
```

Expected: at the current `eps=1e-3`, removing one inhibitor multiplies its
target by **11**. At `1e-9` it multiplies by ten million — which is why the
guard cannot simply be lowered without an explicit ceiling.

## 2. The attribution 2×2 (US1)

One arm at a time, killing the Julia server between them.

```bash
# A — current behaviour
DS_INHIBITOR_EPS=1e-3 DS_DEREPRESSION_MAX=11 ... --output-dir $S/r_A
# B — epsilon's non-ceiling effects, ceiling held
DS_INHIBITOR_EPS=1e-9 DS_DEREPRESSION_MAX=11 ... --output-dir $S/r_B
# C — the ceiling, epsilon held
DS_INHIBITOR_EPS=1e-9 DS_DEREPRESSION_MAX=2  ... --output-dir $S/r_C
# D — the original probe
DS_INHIBITOR_EPS=1e-2                        ... --output-dir $S/r_D
```

Each with `--catalog $S/cat7 --resolve-set-readouts`, on both
`--ground-truth experimental` and `curator`.

**How to read it.** A→B *strengthens* suppression (0.0110 → 0.0100). So:

- B worse than A, C recovering it ⇒ the gain was the **ceiling**.
- B better than A ⇒ part of the gain was **suppression**, and the feature's
  framing needs correcting in `research.md` per FR-010.
- Both move ⇒ report the proportion. Do not pick the flattering half.

Arm A must be **byte-identical** to a run with the feature absent. If it is
not, the default does not reproduce current behaviour and nothing else in
this list means anything.

## 3. Prove the guard is only a guard (US2, FR-002)

```bash
for eps in 1e-6 1e-7 1e-8 1e-9; do
  DS_INHIBITOR_EPS=$eps DS_DEREPRESSION_MAX=2 ... --output-dir $S/r_eps$eps
done
```

Expected: **identical predictions across all four**. Any difference means the
guard is still doing modelling work and FR-002 is unmet.

## 4. Compare arms

```bash
python3 bench/analysis/compare_set_rules.py --baseline $S/r_A --arms $S/r_B $S/r_C $S/r_D
```

Reports macro-F1, coverage delta, per-pathway net and the both-arms-converged
count. Read coverage first; read the both-converged count before believing
any small delta — only 9 of 29 changed predictions converged in both arms in
the original probe.

## 5. The targeted subgroup (FR-009)

Report knockouts predicted to increase, separately, for every arm: how many
wrong calls were removed **and** how many correct ones were lost. Baseline is
28 correct of 78. A bound that removes all 78 would look fine on the total
and be wrong.

## 6. Decide the default (US3)

On `$S/cat92`, both ground truths, with the 742-case set as secondary. A
candidate that improves the small set and not the large one is not adopted,
and the discrepancy is recorded.
