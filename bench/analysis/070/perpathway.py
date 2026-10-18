"""specs/070: per-pathway fixed/broken over all cases valid in both runs.

  perpathway.py BASE_CASES ARM_CASES
"""
import collections, csv, sys

def load(p):
    return {(r["pathway"], r["gene"], r["direction"], r["key_output"]): r
            for r in csv.DictReader(open(p), delimiter="\t")}

a, b = load(sys.argv[1]), load(sys.argv[2])
per = collections.defaultdict(lambda: [0, 0, 0])
for k, x in a.items():
    y = b.get(k)
    if not y or x["valid"] != "1" or y["valid"] != "1":
        continue
    ok_a, ok_b = x["predicted"] == x["expected"], y["predicted"] == y["expected"]
    per[k[0]][2] += 1
    if ok_b and not ok_a: per[k[0]][0] += 1
    if ok_a and not ok_b: per[k[0]][1] += 1
tot = [sum(v[i] for v in per.values()) for i in range(3)]
print(f"all valid-in-both: n={tot[2]} fixed={tot[0]} broken={tot[1]} net={tot[0]-tot[1]}")
for p, (f, br, n) in sorted(per.items(), key=lambda kv: kv[1][0] - kv[1][1]):
    if f or br:
        print(f"{f-br:+5d}  {f:4d} / {br:4d}  of {n:5d}  {p}")
