# specs/077 structural count: open carrier-loop free forms with a producer outside
# their loops, and the release edges DS_CARRIER_LOOP_OPEN=fed drops.
# Usage: python3 count_open_fed.py ~/deltasignal-catalogs/builds/<build>
import csv,sys,glob,os
from collections import defaultdict, Counter
build=sys.argv[1]
names={}
tot=Counter(); per={}
for f in sorted(glob.glob(build+'/R-HSA-*/carrier_loops.csv')):
    d=os.path.dirname(f); pw=os.path.basename(d)
    rows=list(csv.DictReader(open(f)))
    listed=defaultdict(set)   # free -> loop sources
    for r in rows:
        if r['closed']=='False': listed[r['free_uuid']].add(r['loop_source_uuid'])
    if not listed: continue
    pos_in=defaultdict(set)
    for e in csv.DictReader(open(d+'/logic_network.csv')):
        if e['pos_neg']=='pos': pos_in[e['target_id']].add(e['source_id'])
    c=Counter()
    for free,srcs in listed.items():
        outside=pos_in[free]-srcs
        n=len(srcs & pos_in[free])
        if outside: c['kept_supply_forms']+=1; c['edges_dropped']+=n
        else: c['orphan_forms_untouched']+=1; c['edges_all_would_drop']+=n
    per[pw]=c; tot.update(c)
print(dict(tot))
for pw,c in sorted(per.items(), key=lambda x:-x[1]['edges_dropped'])[:12]: print(' ',pw,dict(c))
