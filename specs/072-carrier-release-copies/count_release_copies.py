import csv,sys,os,collections
b=sys.argv[1]; tot=collections.Counter(); per=collections.Counter()
for pw in sorted(os.listdir(b)):
    d=os.path.join(b,pw)
    if not os.path.exists(d+'/pool_carriers.csv'): continue
    st={r['uuid']:r['stable_id'] for r in csv.DictReader(open(d+'/stid_to_uuid_mapping.csv'))}
    prod=collections.defaultdict(set)
    for r in csv.DictReader(open(d+'/logic_network.csv')):
        if r['pos_neg']=='pos': prod[r['target_id']].add(r['source_id'])
    rel=collections.defaultdict(set)
    for r in csv.DictReader(open(d+'/pool_carriers.csv')): rel[r['carrier_uuid']].add(r['release_reaction_uuid'])
    for c,rs in rel.items():
        ss={st.get(x) for x in rs}
        extra={a for a in prod[c] if a not in rs and st.get(a) in ss}
        tot['carriers']+=1
        if extra: tot['carriers_with_extra']+=1; tot['extra_edges']+=len(extra); per[pw]+=1
print(b,dict(tot)); print(per.most_common(30))
