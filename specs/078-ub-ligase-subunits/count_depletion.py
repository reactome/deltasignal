# specs/078: depletion edges a variant build removes against canonical, compared by
# entity (uuids are re-minted per build). Also checks every other edge class is unchanged.
# Usage: python3 count_depletion.py <canonical build> <variant build> [<variant build> ...]
import csv, os, sys
from collections import Counter
def edges(b, pw):
    d = f'{b}/{pw}'
    ent = {r['uuid']: r['diagram_entity_id'] for r in csv.DictReader(open(d + '/nodes.csv'))}
    dep, other = Counter(), Counter()
    for f in ('logic_network.csv', 'boundary_edges.csv'):
        if not os.path.exists(f'{d}/{f}'): continue
        for e in csv.DictReader(open(f'{d}/{f}')):
            k = (ent.get(e['source_id'], '?'), ent.get(e['target_id'], '?'))
            (dep if e['edge_type'] == 'depletion' else other)[(k, e['edge_type'])] += 1
    return dep, other
base = sys.argv[1]
pws = sorted(p for p in os.listdir(base) if p.startswith('R-HSA-'))
for v in sys.argv[2:]:
    tot_b = tot_v = 0; per = Counter(); other_diff = 0; added = 0; cks = []
    for pw in pws:
        db, ob = edges(base, pw); dv, ov = edges(v, pw)
        tot_b += sum(db.values()); tot_v += sum(dv.values())
        per[pw] = sum((db - dv).values()); added += sum((dv - db).values())
        other_diff += sum(((ob - ov) + (ov - ob)).values())
        cks += [(pw, k) for k in (db - dv) if k[0][1] == 'R-HSA-187533']
    print(os.path.basename(v), 'depletion edges', tot_b, '->', tot_v, 'removed', sum(per.values()), 'added', added,
          'pathways', sum(1 for n in per.values() if n), 'other-edge differences', other_diff)
    for pw, n in per.most_common(8): print('   ', pw, n)
    print('    free CKS1B (R-HSA-187533) drains removed:', sorted(set((p, k[0][0]) for p, k in cks)))
