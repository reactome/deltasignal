"""Simple cycles through a node (bounded length), by the parity of negative edges, with the reactions each passes."""
import sys, csv, collections, itertools
sys.path.insert(0, '/home/awright/gitroot/deltasignal/bench')
import benchmark_vs_mpbiopath as bm
import networkx as nx
pw, node, maxlen = sys.argv[1], sys.argv[2], int(sys.argv[3])
drop_rx = set(sys.argv[4].split(',')) if len(sys.argv) > 4 else set()
pdir = bm.CATALOG_ROOT / pw
nodes = {r['uuid']: r for r in csv.DictReader(open(pdir / 'nodes.csv'))}
G = nx.DiGraph()
with bm.open_network(pdir) as fh:
    for e in csv.DictReader(fh):
        if e['edge_type'] in ('input', 'output') and any(nodes.get(x, {}).get('diagram_entity_id') in drop_rx and nodes.get(x, {}).get('node_kind') == 'reaction' for x in (e['source_id'], e['target_id'])):
            continue
        neg = e['pos_neg'] == 'neg' or e['edge_type'] == 'depletion'
        if G.has_edge(e['source_id'], e['target_id']):
            G[e['source_id']][e['target_id']]['signs'].add((neg, e['edge_type']))
        else:
            G.add_edge(e['source_id'], e['target_id'], signs={(neg, e['edge_type'])})
u = [x for x in nodes if x.startswith(node)][0]
comp = next(c for c in nx.strongly_connected_components(G) if u in c)
print('component size', len(comp))
H = G.subgraph(comp)
par = collections.Counter(); ex = {}
n = 0
for cyc in nx.simple_cycles(H, length_bound=maxlen):
    if u not in cyc: continue
    n += 1
    pairs = list(zip(cyc, cyc[1:] + cyc[:1]))
    # an edge pair with both signs makes the cycle ambiguous
    opts = [sorted({s for s, _ in H[a][b]['signs']}) for a, b in pairs]
    for choice in itertools.product(*opts):
        k = sum(choice) % 2
        types = tuple(sorted(t for (a, b), s in zip(pairs, choice) for ss, t in H[a][b]['signs'] if ss == s and ss))
        par[(('even' if k == 0 else 'odd'), types)] += 1
        ex.setdefault((('even' if k == 0 else 'odd'), types), [nodes.get(x, {}).get('diagram_entity_id', '') for x in cyc if nodes.get(x, {}).get('node_kind') == 'reaction'])
    if n > 200000: print('capped'); break
print('cycles through node', n)
for k, v in par.most_common(15): print(v, k, ex[k][:8])
