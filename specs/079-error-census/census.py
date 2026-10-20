# specs/079: what the remaining errors are, on one results directory.
# Usage: python3 census.py <results dir with curator_cases.tsv and experimental_cases.tsv>
import csv, sys
from collections import Counter
R = sys.argv[1]
def load(ax):
    return {(r['pathway'], r['gene'], r['direction'], r['key_output']): r
            for r in csv.DictReader(open(f'{R}/{ax}_cases.tsv'), delimiter='\t') if r['valid'] == '1'}
cu, ex = load('curator'), load('experimental')
for name, d in (('curator', cu), ('experimental', ex)):
    w = [k for k in d if d[k]['predicted'] != d[k]['expected']]
    c = Counter(d[k]['category'] for k in w)
    print(f'{name}: valid {len(d)}, wrong {len(w)}; ' + ', '.join(f'{k} {v}' for k, v in c.most_common()))
    print(f'  ceiling for a propagator on these networks (every routed case right): '
          f'{100 * (len(d) - c["no_path"]) / len(d):.1f}%, now {100 * (len(d) - len(w)) / len(d):.1f}%')
shared = [k for k in ex if k in cu]
agree = [k for k in shared if ex[k]['expected'] == cu[k]['expected']]
dis = [k for k in shared if ex[k]['expected'] != cu[k]['expected']]
print(f'experimental cases also curated: {len(shared)}; curators agree with the experiment on {len(agree)} '
      f'({100 * len(agree) / len(shared):.1f}%)')
ok = lambda k: ex[k]['predicted'] == ex[k]['expected']
print(f'  where they agree: right {sum(map(ok, agree))} of {len(agree)} '
      f'({100 * sum(map(ok, agree)) / len(agree):.1f}%)')
follow = sum(1 for k in dis if ex[k]['predicted'] == cu[k]['expected'])
print(f'  where they disagree ({len(dis)}): we match the experiment {sum(map(ok, dis))}, the curators {follow}, neither {len(dis) - follow - sum(map(ok, dis))}')
w = [k for k in agree if not ok(k)]
print('  errors where both agree, by category', dict(Counter(ex[k]['category'] for k in w)))
print('  by perturbation', Counter((k[0], k[1], k[2]) for k in w).most_common(10))
