# specs/089: cases whose correctness changes between two runs (valid in both).
# Usage: python3 diff_cases.py BASE_cases.tsv ARM_cases.tsv
import csv, sys
from collections import Counter
def load(p):
    return {(r['pathway'], r['gene'], r['direction'], r['key_output']): r
            for r in csv.DictReader(open(p), delimiter='\t')}
a, b = load(sys.argv[1]), load(sys.argv[2])
net = Counter()
for k in sorted(a.keys() & b.keys()):
    x, y = a[k], b[k]
    if x['valid'] != '1' or y['valid'] != '1': continue
    ok0, ok1 = x['predicted'] == x['expected'], y['predicted'] == y['expected']
    if ok0 != ok1:
        net[k[0]] += 1 if ok1 else -1
        print('FIX ' if ok1 else 'BRK ', *k, x['expected'], x['predicted'], '->', y['predicted'],
              x['pred_ui'][:7], '->', y['pred_ui'][:7], sep='\t')
print(dict(net))
