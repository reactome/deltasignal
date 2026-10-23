"""specs/094: experimental errors on which the curators and the experiment agree,
grouped by perturbation. Usage: python3 agreed_errors.py <results dir>"""
import csv, sys
from collections import Counter, defaultdict
R = sys.argv[1]
def load(p):
    return {(r['pathway'], r['gene'], r['direction'], r['key_output']): r
            for r in csv.DictReader(open(p), delimiter='\t')}
e = load(R + '/experimental_cases.tsv'); c = load(R + '/curator_cases.tsv')
g = defaultdict(list)
for k, r in e.items():
    if r['valid'] != '1' or k not in c or c[k]['expected'] != r['expected'] or r['predicted'] == r['expected']:
        continue
    g[k[:3]].append((r['expected'], r['predicted'], r['category'], float(r['pred_ui'])))
total = 0
for (pw, gene, d), v in sorted(g.items(), key=lambda x: -len(x[1])):
    total += len(v)
    print(f"{len(v):2d}\t{pw}\t{gene}\t{'KO' if d == '0' else 'OE'}\t"
          f"{dict(Counter((a, b) for a, b, _, _ in v))}\t{dict(Counter(x[2] for x in v))}")
print("total", total)
