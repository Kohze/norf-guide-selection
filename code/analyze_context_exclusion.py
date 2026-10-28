"""Post hoc, fixed-cohort guide-selection sensitivity. No network or tuning.

Run from any directory: python code/analyze_context_exclusion.py
Reads frozen S25; writes all new outcomes, draws, summaries and source hashes.
"""
from pathlib import Path
from itertools import combinations
import json, hashlib, platform
from datetime import datetime, timezone
import numpy as np
import pandas as pd

M = Path(__file__).resolve().parents[1]
O = M / 'evidence/selection_context_sensitivity'
O.mkdir(parents=True, exist_ok=True)
source = next((M / 'supplementary').glob('S25_*'))
d = pd.read_csv(source, float_precision='round_trip')
cells = ['A375', 'A549', 'HEPG2']
assert len(d) == 1251 and d.pair_ID.nunique() == 87
assert d.sgrna.nunique() == 417 and set(d.cell) == set(cells)
assert not d.duplicated(['pair_ID', 'sgrna', 'cell']).any()
assert np.isfinite(d[['old_LFC', 'new_LFC']]).all().all()
rows = []
expectation_error = 0.
for pid, group in d.groupby('pair_ID', sort=True):
    old = group.pivot(index='sgrna', columns='cell', values='old_LFC')[cells].sort_index()
    new = group.pivot(index='sgrna', columns='cell', values='new_LFC').reindex(old.index)[cells]
    assert len(old) >= 4 and old.notna().all().all() and new.notna().all().all()
    clusters = group.dependence_cluster.unique()
    assert len(clusters) == 1
    for rule in ['original_all_three_cells', 'exclude_evaluated_cell']:
        for cell in cells:
            training = cells if rule == 'original_all_three_cells' else [c for c in cells if c != cell]
            score = old[training].mean(axis=1)
            chosen = sorted(old.index, key=lambda seq: (score[seq], seq))[:2]
            reference = float(new[cell].mean())
            exact = float(np.mean([new.loc[list(pair), cell].mean() for pair in combinations(old.index, 2)]))
            expectation_error = max(expectation_error, abs(exact-reference))
            selected = float(new.loc[chosen, cell].mean())
            if rule == 'original_all_three_cells':
                frozen = set(group.loc[(group.cell == cell) & group.selected_two, 'sgrna'])
                assert set(chosen) == frozen
            rows.append(dict(rule=rule,pair_ID=pid,cell=cell,
                selection_cells=';'.join(training),guide_count=len(old),
                selected_guide_1=chosen[0],selected_guide_2=chosen[1],
                newer_selected_two_LFC=selected,newer_random_two_expectation=reference,
                newer_selection_shift=selected-reference,dependence_cluster=clusters[0]))
t = pd.DataFrame(rows)
t.to_csv(O/'all522_original_and_context_excluded_product_cell_contrasts.csv',index=False)
products = t.groupby(['rule','pair_ID','dependence_cluster'],sort=True).newer_selection_shift.mean().reset_index()
products.to_csv(O/'all174_product_average_contrasts.csv',index=False)
ids = sorted(d.pair_ID.unique())
assert d.groupby('pair_ID').dependence_cluster.nunique().eq(1).all()
assert d[['pair_ID','dependence_cluster']].drop_duplicates().dependence_cluster.nunique() == 87
rng = np.random.default_rng(20261026)
draws = rng.integers(0, len(ids), size=(2000,len(ids)))
np.savez_compressed(O/'all2000_joint_product_draws.npz',draws=draws,pair_ID=np.array(ids))
summary = []
bootstrap = []
for rule in t.rule.unique():
    matrix = t[t.rule == rule].pivot(index='pair_ID',columns='cell',values='newer_selection_shift').reindex(index=ids,columns=cells).to_numpy()
    values = np.column_stack([matrix,matrix.mean(axis=1)])
    boot = values[draws].mean(axis=1)
    for j, cell in enumerate(cells+['equal_cell_mean']):
        v = values[:,j]; ci = np.quantile(boot[:,j],[.025,.975])
        summary.append(dict(rule=rule,evaluation=cell,products=87,dependence_groups=87,
            mean_shift=float(v.mean()),nominal95_lower=float(ci[0]),nominal95_upper=float(ci[1]),
            median_product_shift=float(np.median(v)),negative_products=int((v<0).sum()),
            zero_products=int((v==0).sum()),positive_products=int((v>0).sum()),
            post_hoc_sensitivity=True,bootstrap_seed=20261026,bootstrap_samples=2000))
        bootstrap.extend(dict(rule=rule,evaluation=cell,draw=i,mean_shift=float(z)) for i,z in enumerate(boot[:,j]))
s = pd.DataFrame(summary)
s.to_csv(O/'all8_selection_strategy_summary_estimates.csv',index=False)
pd.DataFrame(bootstrap).to_csv(O/'all16000_strategy_bootstrap_estimates.csv',index=False)
v = products[products.rule == 'exclude_evaluated_cell'].set_index('pair_ID').loc[ids].newer_selection_shift.to_numpy()
pd.DataFrame({'omitted_pair_ID':ids,'mean_shift':[(v.sum()-x)/(len(v)-1) for x in v]}).to_csv(O/'all87_context_exclusion_product_omissions.csv',index=False)
frozen = json.loads((M/'evidence/analysis_extensions/cycle_93_transferred_guide_selection/GUIDE_SELECTION_RESULTS.json').read_text())
original = s[(s.rule == 'original_all_three_cells') & (s.evaluation == 'equal_cell_mean')].iloc[0]
assert abs(original.mean_shift-frozen['mean_newer_selection_shift']) < 1e-12
for cell in cells:
    point = s[(s.rule == 'original_all_three_cells') & (s.evaluation == cell)].iloc[0].mean_shift
    assert abs(point-frozen['cell_specific_newer_shifts'][cell]) < 1e-12
assert expectation_error < 1e-12
result = {'utc':datetime.now(timezone.utc).isoformat(),'status':'computed_pending_separate_verification',
    'source':str(source.relative_to(M)),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
    'protocol_sha256':hashlib.sha256((M/'provenance/NEW_SENSITIVITY_PROTOCOL.json').read_bytes()).hexdigest(),
    'products':87,'guides':417,'cells':cells,'original_point_reproduced':True,
    'exact_pair_expectation_max_abs_error':expectation_error,'seed':20261026,'bootstrap_samples':2000,
    'software':{'Python':platform.python_version(),'NumPy':np.__version__,'pandas':pd.__version__},
    'estimates':s.to_dict(orient='records'),
    'limits':['Post hoc sensitivity; all contexts and underlying outcomes already known',
      'Evaluated cell excluded only from earlier selection score; no unseen biological context',
      'Different evaluated cells may select different pairs',
      'Nominal intervals are not adjusted over the retrospective research program',
      'No optimal strategy, molecular cause or false-positive target established']}
(O/'CONTEXT_EXCLUSION_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':result['status'],'summary':s.to_dict(orient='records')}))
