# specs/080: summarise reach_<sha>.tsv. Usage: python3 summarise.py reach_3e96e7c.tsv
import csv, sys
from collections import Counter
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
for ax in ("curator", "experimental"):
    a = [r for r in rows if r["axis"] == ax]
    np_ = [r for r in a if r["category"] == "no_path"]
    print(f"== {ax}: {len(a)} valid cases; ours reachable {sum(r['our_reach']=='1' for r in a)}")
    for v in ("v1", "v2", "v3", "v4", "v5"):
        k = f"db_{v}"
        both = Counter((r[k], r["our_reach"]) for r in a)
        print(f"  {v}: db reachable {sum(r[k]=='1' for r in a)}; db yes/ours no {both[('1','0')]}; "
              f"db no/ours yes {both[('0','1')]}; of {len(np_)} wrong no_path, db reachable {sum(r[k]=='1' for r in np_)}")
    if ax == "experimental":
        m = [r for r in a if r["mpb_path"] != ""]
        print(f"  MP-BioPath distance table: path {sum(r['mpb_path']=='1' for r in m)} of {len(m)}; "
              f"mpb yes/ours no {sum(r['mpb_path']=='1' and r['our_reach']=='0' for r in m)}; "
              f"mpb yes/db_v5 no {sum(r['mpb_path']=='1' and r['db_v5']=='0' for r in m)}")
