# specs/080: summarise breaks_<sha>.tsv. Usage: python3 summarise_breaks.py breaks_3e96e7c.tsv
import csv, sys
from collections import Counter
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
c = Counter(); wr = Counter()
for r in rows:
    c[r["break_kind"]] += int(r["cases"]); wr[r["break_kind"]] += int(r["wrong"])
print("by kind of step (cases, wrong):")
for k, n in c.most_common(): print(f"  {n:5d} {wr[k]:4d}  {k}")
bp = Counter(); bw = Counter()
for r in rows:
    k = (r["pathway_id"], r["last_reached"], r["first_missed"], r["break_kind"])
    bp[k] += int(r["cases"]); bw[k] += int(r["wrong"])
print(f"distinct break points {len(bp)}; by wrong cases:")
for k, n in sorted(bp.items(), key=lambda x: -bw[x[0]])[:15]: print(f"  {n:4d} ({bw[k]:3d} wrong) {k}")
