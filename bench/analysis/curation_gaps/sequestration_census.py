import csv,os,sys,collections,json
sys.path.insert(0,'bench/analysis'); from holdout_report import TUNING_PATHWAYS
C=os.path.expanduser('~/deltasignal-catalogs/builds/20261012-1401_01206ca/')
B=C+'results/64d23cb/curator_cases.tsv'
pid={r[1]:r[0] for r in csv.reader(open('bench/catalog_pathways.tsv'),delimiter='\t') if r and r[0].startswith('R-HSA')}
rows=[r for r in csv.DictReader(open(B),delimiter='\t') if r['valid']=='1' and not r['exclusion'] and r['pathway'] not in TUNING_PATHWAYS]
nopath=collections.defaultdict(list)
for r in rows:
    if r['category']=='no_path': nopath[(r['pathway'],r['gene'])].append(r)
cache={}
def net(p):
    if p in cache: return cache[p]
    d=C+pid[p]; kind={}; 
    for n in csv.DictReader(open(d+'/nodes.csv')): kind[n['uuid']]=n['node_kind']
    cof=set()
    if os.path.exists(d+'/cofactors.csv'):
        st={a:b for a,b in csv.reader(open(d+'/stid_to_uuid_mapping.csv'))}
        cs={r[0] for r in csv.reader(open(d+'/cofactors.csv'))}
        cof={u for u,s in st.items() if s in cs}
    out=collections.defaultdict(list); inc=collections.defaultdict(list)
    for f in ('logic_network.csv','boundary_edges.csv'):
        if os.path.exists(d+'/'+f):
            for e in csv.DictReader(open(d+'/'+f)):
                out[e['source_id']].append((e['target_id'],e['edge_type'])); inc[e['target_id']].append((e['source_id'],e['edge_type']))
    cache[p]=(kind,out,inc,cof); return cache[p]
def reach(out,starts,skip=('dissociation',)):
    seen=set(starts); st=list(starts)
    while st:
        u=st.pop()
        for v,et in out[u]:
            if et not in skip and v not in seen: seen.add(v); st.append(v)
    return seen
tot=collections.Counter(); detail=[]
for (pw,g),rs in nopath.items():
    kind,out,inc,cof=net(pw)
    pins=[u for u in rs[0]['gene_uuids'].split('|') if u]
    R=reach(out,pins)
    sinks=[u for u in R if kind.get(u) not in ('reaction','dissociation') and u not in cof and not any(et!='dissociation' for _,et in out[u])]
    # partners: inputs of the reactions producing each sink, not reached from the pins
    partners=set()
    for s in sinks:
        for rx,et in inc[s]:
            if kind.get(rx)!='reaction': continue
            for p,et2 in inc[rx]:
                if et2=='input' and p not in R and p not in cof: partners.add(p)
    PR=reach(out,partners) if partners else set()
    hit=[r for r in rs if any(o in PR for o in r['output_uuids'].split('|'))]
    tot['cases']+=len(rs); tot['partner_reaches_readout']+=len(hit)
    detail.append((len(hit),len(rs),pw,g))
print(dict(tot))
for x in sorted(detail,reverse=True)[:25]: print(x)
