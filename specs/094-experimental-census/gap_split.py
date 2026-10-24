"""specs/094: the experimental gap to MP-BioPath split by whether the curators agree
with the experiment. Reads bench/analysis/experimental_gap.py's --out file and the
arm's curator_cases.tsv. Usage: python3 gap_split.py GAP.tsv ARM_RESULTS_DIR"""
import csv, sys
from collections import Counter
gap, R = sys.argv[1:]
cur = {(r['pathway'], r['gene'], r['direction'], r['key_output']): r['expected']
       for r in csv.DictReader(open(R + '/curator_cases.tsv'), delimiter='\t')}
for q, name in (('MPB/-', 'only MP-BioPath right'), ('-/DS', 'only DeltaSignal right')):
    c = Counter(); per = Counter()
    for x in csv.DictReader(open(gap), delimiter='\t'):
        if x['quadrant'] != q:
            continue
        ct = cur.get((x['pathway'], x['gene'], x['direction'], x['key_output']))
        tag = 'curators agree' if ct == x['truth'] else ('no curator case' if ct is None else 'curators disagree')
        c[tag] += 1
        if tag == 'curators agree':
            per[(x['pathway'], x['gene'], 'KO' if x['direction'] == '0' else 'OE')] += 1
    print(name, dict(c))
    print('   agreed cases by perturbation:', per.most_common(12))
