from pathlib import Path
import datetime,hashlib,itertools,json,math
import numpy as np
import pandas as pd
M=Path(__file__).resolve().parent
N=M.parent if (M.parent/'cycle_90_within_product_shared_guide_concordance').exists() else M.parents[1]/'exploration_2026-10-11'/'novelty_discovery'
SRC=N/'cycle_90_within_product_shared_guide_concordance'/'all1413_primary_shared_guide_cell_paired_residuals.csv'
d=pd.read_csv(SRC);n=d.groupby('pair_ID').sgrna.nunique();d=d[d.pair_ID.isin(n[n>=4].index)].copy()
assert (d.pair_ID.nunique(),d.sgrna.nunique(),len(d))==(87,417,1251)
old=d.groupby(['pair_ID','sgrna'],as_index=False).old_LFC.mean().rename(columns={'old_LFC':'older_equal_cell_mean_LFC'})
old=old.sort_values(['pair_ID','older_equal_cell_mean_LFC','sgrna'],kind='stable')
old['older_rank']=old.groupby('pair_ID').cumcount()+1;old['selected_two']=old.older_rank.le(2)
d=d.merge(old,on=['pair_ID','sgrna'],validate='many_to_one')
d.to_csv(M/'all1251_guide_cell_outcomes_and_older_only_selections.csv',index=False)
allmeans=d.groupby(['pair_ID','cell']).agg(all_newer_LFC=('new_LFC','mean'),all_older_LFC=('old_LFC','mean'),dependence_cluster=('dependence_cluster','first'),older_source_name=('old_source_name_old','first'),newer_source_ORF_ID=('new_source_ORF_ID_old','first'),guide_count=('sgrna','nunique'))
sel=d[d.selected_two].groupby(['pair_ID','cell']).agg(selected_newer_LFC=('new_LFC','mean'),selected_older_LFC=('old_LFC','mean'),selected_sequences=('sgrna',lambda x:';'.join(sorted(x))))
a=allmeans.join(sel).reset_index();a['newer_selection_shift']=a.selected_newer_LFC-a.all_newer_LFC;a['older_training_shift']=a.selected_older_LFC-a.all_older_LFC
assert a.groupby('pair_ID').size().eq(3).all()
a.to_csv(M/'all261_product_cell_selection_contrasts.csv',index=False)
p=a.groupby('pair_ID').agg(newer_selection_shift=('newer_selection_shift','mean'),older_training_shift=('older_training_shift','mean'),all_newer_LFC=('all_newer_LFC','mean'),selected_newer_LFC=('selected_newer_LFC','mean'),guide_count=('guide_count','first'),dependence_cluster=('dependence_cluster','first'),older_source_name=('older_source_name','first'),newer_source_ORF_ID=('newer_source_ORF_ID','first'),selected_sequences=('selected_sequences','first')).reset_index()
k=math.ceil(.2*len(p))
base=set(p.sort_values(['all_newer_LFC','pair_ID']).head(k).pair_ID)
chosen=set(p.sort_values(['selected_newer_LFC','pair_ID']).head(k).pair_ID)
p['allguide_top20pct']=p.pair_ID.isin(base);p['selected2_top20pct']=p.pair_ID.isin(chosen)
clusters=sorted(p.dependence_cluster.unique());idx=[np.flatnonzero(p.dependence_cluster.to_numpy()==c) for c in clusters]
rng=np.random.default_rng(20261018);draws=rng.integers(0,len(idx),size=(2000,len(idx)))
boot=np.array([p.newer_selection_shift.to_numpy()[np.concatenate([idx[j] for j in draw])].mean() for draw in draws])
np.savez_compressed(M/'all2000_selection_cluster_draws.npz',draws=draws,clusters=np.array(clusters))
pd.DataFrame({'draw':np.arange(2000),'mean_newer_selection_shift':boot}).to_csv(M/'all2000_selection_cluster_bootstrap.csv',index=False)
# Fixed random-pair plans retain whole profiles across the fixed cell contexts.
matrices={pid:q.pivot(index='sgrna',columns='cell',values='new_LFC').to_numpy() for pid,q in d.groupby('pair_ID')}
ids=list(p.pair_ID);rng=np.random.default_rng(20261019);random_results=[];topcount=dict.fromkeys(ids,0);plans=[]
for j in range(1000):
 vals=[];pairs=[]
 for pid in ids:
  mat=matrices[pid];chosen_ix=rng.choice(len(mat),2,replace=False);pairs.append(chosen_ix.tolist());vals.append(mat[chosen_ix].mean())
 vals=np.array(vals);random_top={ids[i] for i in np.lexsort((np.array(ids),vals))[:k]}
 for pid in random_top:topcount[pid]+=1
 random_results.append({'draw':j,'mean_selection_shift':float(np.mean(vals-p.all_newer_LFC.to_numpy())),'top20pct_reference_overlap':len(base&random_top)})
 plans.append(pairs)
np.savez_compressed(M/'all1000_random2_selection_plans.npz',plans=np.array(plans),pair_ids=np.array(ids))
p['random2_top20pct_selection_fraction']=p.pair_ID.map(topcount)/1000
p.to_csv(M/'all87_product_selection_shifts_and_rank_sensitivity.csv',index=False)
pd.DataFrame(random_results).to_csv(M/'all1000_random2_selection_results.csv',index=False)
omissions=pd.DataFrame([{'omitted_pair_ID':pid,'mean_shift':p.loc[p.pair_ID.ne(pid),'newer_selection_shift'].mean()} for pid in ids])
omissions.to_csv(M/'all87_single_product_omissions.csv',index=False)
# Exact enumeration establishes the count-matched reference for every product-cell.
maxerr=0
for (pid,cell),q in d.groupby(['pair_ID','cell']):
 vals=q.new_LFC.to_numpy();expected=np.mean([(vals[i]+vals[j])/2 for i,j in itertools.combinations(range(len(vals)),2)])
 maxerr=max(maxerr,abs(expected-vals.mean()))
assert maxerr<1e-12
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'analyzed_pending_separate_verification','products':len(p),'guides':int(d.sgrna.nunique()),'guide_cell_pairs':len(d),'dependence_clusters':len(clusters),'mean_newer_selection_shift':float(p.newer_selection_shift.mean()),'nominal95_cluster_CI':np.quantile(boot,[.025,.975]).tolist(),'older_training_shift':float(p.older_training_shift.mean()),'cell_specific_newer_shifts':a.groupby('cell').newer_selection_shift.mean().to_dict(),'single_product_omission_range':[float(omissions.mean_shift.min()),float(omissions.mean_shift.max())],'rank_sensitivity_descriptive':{'priority_budget':k,'reference_and_selected2_overlap':len(base&chosen),'reference_and_selected2_replacements':k-len(base&chosen),'random2_reference_overlap_median':float(np.median([x['top20pct_reference_overlap'] for x in random_results])),'random2_reference_overlap_minmax':[int(min(x['top20pct_reference_overlap'] for x in random_results)),int(max(x['top20pct_reference_overlap'] for x in random_results))]},'exact_random2_expectation_max_difference':maxerr,'source':{'path':str(SRC),'bytes':SRC.stat().st_size,'sha256':hashlib.sha256(SRC.read_bytes()).hexdigest()},'limits':json.loads((M/'FROZEN_GUIDE_SELECTION_SCENARIO.json').read_text())['limits']}
(M/'GUIDE_SELECTION_RESULTS.json').write_text(json.dumps(result,indent=2),encoding='utf8');print(json.dumps(result,indent=2))
