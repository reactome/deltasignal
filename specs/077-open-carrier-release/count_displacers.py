# specs/077 research: among the open-loop release edges DS_CARRIER_LOOP_OPEN=fed drops,
# how many come from a reaction with a second non-cofactor input (a displacer).
# Usage: python3 count_displacers.py <build dir> [<build dir> ...]
import csv,sys,glob,os
from collections import defaultdict, Counter
for b in sys.argv[1:]:
    tot=Counter(); per=Counter(); ex={}
    for f in glob.glob(b+'/R-HSA-*/carrier_loops.csv'):
        d=os.path.dirname(f); pw=os.path.basename(d)
        rows=list(csv.DictReader(open(f)))
        cof=set()
        if os.path.exists(d+'/cofactors.csv'):
            for r in csv.DictReader(open(d+'/cofactors.csv')): cof.add(list(r.values())[0])
        nodes={r['uuid']:r for r in csv.DictReader(open(d+'/nodes.csv'))}
        listed=defaultdict(set)
        for r in rows: listed[r['free_uuid']].add(r['loop_source_uuid'])
        pos=defaultdict(set); inputs=defaultdict(list)
        for e in csv.DictReader(open(d+'/logic_network.csv')):
            if e['pos_neg']=='pos': pos[e['target_id']].add(e['source_id'])
            inputs[e['target_id']].append(e)
        for r in rows:
            fr,src=r['free_uuid'],r['loop_source_uuid']
            if r['closed']!='False' or not (pos[fr]-listed[fr]) or src not in pos[fr]: continue
            ins=[e for e in inputs[src] if e['pos_neg']=='pos' and e['edge_type']=='input']
            noncof=[e for e in ins if e['source_id'] not in cof and nodes.get(e['source_id'],{}).get('diagram_entity_id') not in cof]
            k='displacer' if len(noncof)>1 else 'single'
            tot[k]+=1; per[(pw,k)]+=1
            if pw in ('R-HSA-109606','R-HSA-5693532') : ex.setdefault((pw,k),set()).add(nodes.get(src,{}).get('diagram_entity_id'))
    print(os.path.basename(b), dict(tot))
    for (pw,k),n in sorted(per.items(), key=lambda x:-x[1])[:14]: print('  ',pw,k,n)
    for k,v in ex.items(): print('  ex',k,sorted(v)[:8])
