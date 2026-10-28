from pathlib import Path
import itertools,json,datetime
import numpy as np
import pandas as pd
M=Path(__file__).resolve().parent;N=M.parent if (M.parent/'cycle_87_fresh_guide_cross_screen_nORF_validation').exists() else M.parents[1]/'exploration_2026-10-11'/'novelty_discovery';C87=N/'cycle_87_fresh_guide_cross_screen_nORF_validation'
rows=pd.read_csv(M/'all1251_guide_cell_outcomes_and_older_only_selections.csv');p=pd.read_csv(M/'all87_product_selection_shifts_and_rank_sensitivity.csv');out=json.loads((M/'GUIDE_SELECTION_RESULTS.json').read_text());checks={}
def ck(k,b):
 checks[k]=bool(b)
 if not b:raise AssertionError(k)
old=pd.read_csv(C87/'all3030_older_exact_product_source_guide_cell_outcomes.csv').set_index(['pair_ID','sgrna','cell'])
new=pd.read_csv(C87/'all2160_newer_exact_product_source_guide_cell_outcomes_and_roles.csv').set_index(['pair_ID','sgrna','cell'])
expected=[];cellvals=[];max_ref_error=0
for pid,q in rows.groupby('pair_ID'):
 guides=sorted(q.sgrna.unique());cells=sorted(q.cell.unique())
 train={g:np.mean([old.loc[(pid,g,c),'LFC'] for c in cells]) for g in guides}
 selected=sorted(guides,key=lambda g:(train[g],g))[:2]
 obs=set(q.loc[q.selected_two,'sgrna'])
 ck('older_only_selection_'+pid,set(selected)==obs)
 for cell in cells:
  allvals=np.array([new.loc[(pid,g,cell),'LFC'] for g in guides])
  sval=np.mean([new.loc[(pid,g,cell),'LFC'] for g in selected]);ref=np.mean([np.mean(allvals[list(z)]) for z in itertools.combinations(range(len(guides)),2)])
  max_ref_error=max(max_ref_error,abs(ref-allvals.mean()))
  cellvals.append({'pair_ID':pid,'cell':cell,'shift':sval-ref})
 expected.append({'pair_ID':pid,'shift':np.mean([x['shift'] for x in cellvals if x['pair_ID']==pid])})
expected=pd.DataFrame(expected).set_index('pair_ID');p=p.set_index('pair_ID')
ck('all87_product_shifts_reconstructed_from_pinned_outcomes',np.max(np.abs(expected.loc[p.index,'shift']-p.newer_selection_shift))<1e-12)
ck('exact_random2_reference_mean_for_every_product_cell',max_ref_error<1e-12)
ck('primary_mean_reproduced',abs(expected['shift'].mean()-out['mean_newer_selection_shift'])<1e-12)
plans=np.load(M/'all2000_selection_cluster_draws.npz');boot=[]
for draw in plans['draws']:
 vals=[]
 for j in draw:vals.extend(expected.loc[p.index[p.dependence_cluster.eq(plans['clusters'][j])],'shift'])
 boot.append(np.mean(vals))
stored=pd.read_csv(M/'all2000_selection_cluster_bootstrap.csv').mean_newer_selection_shift.to_numpy()
ck('all2000_target_bootstrap_draws_reproduced',np.max(np.abs(stored-boot))<1e-12)
ck('nominal95_interval_reproduced',np.allclose(np.quantile(boot,[.025,.975]),out['nominal95_cluster_CI'],atol=1e-12))
randomplans=np.load(M/'all1000_random2_selection_plans.npz');randomstored=pd.read_csv(M/'all1000_random2_selection_results.csv');ids=randomplans['pair_ids'];budget=out['rank_sensitivity_descriptive']['priority_budget'];base=set(p.sort_values(['all_newer_LFC','pair_ID']).head(budget).index)
matrices={pid:q.pivot(index='sgrna',columns='cell',values='new_LFC').to_numpy() for pid,q in rows.groupby('pair_ID')};maxdiff=0
for i,plan in enumerate(randomplans['plans']):
 vals=np.array([matrices[pid][pair].mean() for pid,pair in zip(ids,plan)]);shift=np.mean(vals-p.loc[ids,'all_newer_LFC'].to_numpy());maxdiff=max(maxdiff,abs(shift-randomstored.iloc[i].mean_selection_shift))
 ranked=set(ids[np.lexsort((ids,vals))[:budget]])
 ck('random_rank_plan_'+str(i),len(ranked&base)==randomstored.iloc[i].top20pct_reference_overlap)
ck('all1000_random2_mean_shifts_reproduced',maxdiff<1e-12)
report={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'pass','verification':'Separate implementation by same assistant;not external expert review or biological validation.','passed_checks':sum(checks.values()),'all_checks_passed':all(checks.values()),'selection_checks':87,'random_rank_plan_checks':1000,'primary_mean_shift':float(expected['shift'].mean()),'max_bootstrap_difference':float(np.max(np.abs(stored-boot))),'max_random2_shift_difference':maxdiff,'exact_uniform_random2_reference_max_difference':max_ref_error}
(M/'SEPARATE_COMPUTATIONAL_VERIFICATION.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
