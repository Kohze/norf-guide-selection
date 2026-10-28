"""Separately express the new selection calculations using CSV and Python loops."""
from pathlib import Path
import csv, json, statistics, hashlib
from collections import defaultdict
import numpy as np

M=Path(__file__).resolve().parents[1];O=M/'evidence/selection_context_sensitivity'
raw=list(csv.DictReader(next((M/'supplementary').glob('S25_*')).open(encoding='utf-8',newline='')))
by_product=defaultdict(list)
for row in raw:by_product[row['pair_ID']].append(row)
record=list(csv.DictReader((O/'all522_original_and_context_excluded_product_cell_contrasts.csv').open(encoding='utf-8',newline='')))
lookup={(r['rule'],r['pair_ID'],r['cell']):r for r in record}
points=defaultdict(list);max_error=0.
for pid,group in sorted(by_product.items()):
    guide_cells=defaultdict(dict)
    for r in group:guide_cells[r['sgrna']][r['cell']]=(float(r['old_LFC']),float(r['new_LFC']))
    for rule in ['original_all_three_cells','exclude_evaluated_cell']:
        shifts=[]
        for cell in ['A375','A549','HEPG2']:
            train=['A375','A549','HEPG2'] if rule=='original_all_three_cells' else [c for c in ['A375','A549','HEPG2'] if c!=cell]
            ranked=sorted(guide_cells,key=lambda g:(statistics.mean(guide_cells[g][c][0] for c in train),g))
            chosen=ranked[:2]
            shift=statistics.mean(guide_cells[g][cell][1] for g in chosen)-statistics.mean(x[cell][1] for x in guide_cells.values())
            saved=lookup[(rule,pid,cell)]
            assert chosen==[saved['selected_guide_1'],saved['selected_guide_2']]
            max_error=max(max_error,abs(shift-float(saved['newer_selection_shift'])))
            points[(rule,cell)].append(shift);shifts.append(shift)
        points[(rule,'equal_cell_mean')].append(statistics.mean(shifts))
drawfile=np.load(O/'all2000_joint_product_draws.npz');draws=drawfile['draws']
assert list(drawfile['pair_ID'])==sorted(by_product)
saved_boot=list(csv.DictReader((O/'all16000_strategy_bootstrap_estimates.csv').open(encoding='utf-8',newline='')))
boot_lookup={(r['rule'],r['evaluation'],int(r['draw'])):float(r['mean_shift']) for r in saved_boot}
for (rule,cell),values in points.items():
    for i,draw in enumerate(draws):
        actual=statistics.mean(values[int(j)] for j in draw)
        max_error=max(max_error,abs(actual-boot_lookup[(rule,cell,i)]))
assert max_error<1e-12,max_error
for r in csv.DictReader((O/'all8_selection_strategy_summary_estimates.csv').open(encoding='utf-8',newline='')):
    values=points[(r['rule'],r['evaluation'])]
    assert abs(statistics.mean(values)-float(r['mean_shift']))<1e-12
    assert sum(v<0 for v in values)==int(r['negative_products'])
    assert abs(statistics.median(values)-float(r['median_product_shift']))<1e-12
    ordered=sorted(boot_lookup[(r['rule'],r['evaluation'],i)] for i in range(2000))
    def percentile(q):
        position=q*(len(ordered)-1);i=int(position);fraction=position-i
        return ordered[i]+fraction*(ordered[min(i+1,len(ordered)-1)]-ordered[i])
    assert abs(percentile(.025)-float(r['nominal95_lower']))<1e-12
    assert abs(percentile(.975)-float(r['nominal95_upper']))<1e-12
receipt={'status':'pass','all522_selections_and_contrasts_recomputed':True,
 'all16000_bootstrap_estimates_recomputed':True,'all8_summaries_recomputed':True,'all8_intervals_checked_by_manual_interpolation':True,
 'max_abs_difference':max_error,'scope':'Independent expression of calculations within this project; not external scientific validation',
 'source_sha256':hashlib.sha256(next((M/'supplementary').glob('S25_*')).read_bytes()).hexdigest()}
(O/'SEPARATE_IMPLEMENTATION_VERIFICATION.json').write_text(json.dumps(receipt,indent=2)+'\n')
result=json.loads((O/'CONTEXT_EXCLUSION_RESULTS.json').read_text())
result.update(status='verified_by_separate_expression',verification='SEPARATE_IMPLEMENTATION_VERIFICATION.json')
(O/'CONTEXT_EXCLUSION_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(receipt))
