"""Independent endpoint replay from released analysis-ready supplementary tables."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
M=Path(__file__).resolve().parents[1]
E=M/'evidence/analysis_extensions'
read=lambda n:pd.read_csv(next((M/'supplementary').glob(n+'_*.csv')),float_precision='round_trip')
expected=json.loads((M/'REPORTED_RESULTS.json').read_text())
checks={}
def check(k,v):checks[k]=bool(v);assert v,k
def close(k,a,b):check(k,np.allclose(a,b,rtol=0,atol=2e-12))
def flag(s):return s.astype(str).str.lower().eq('true')
def conditional(d,x,y,seed):
    d=d.copy();cells=d.cell.nunique();n=d.groupby('pair_ID').sgrna.nunique()
    for source,col in [('x',x),('y',y)]:d[source]=d[col]-d.groupby(['pair_ID','cell'])[col].transform('mean')
    w=1/(cells*d.pair_ID.map(n));xx=w*d.x*d.x;yy=w*d.y*d.y;xy=w*d.x*d.y
    a=pd.DataFrame({'cluster':d.dependence_cluster,'xx':xx,'yy':yy,'xy':xy}).groupby('cluster',sort=True)[['xx','yy','xy']].sum().to_numpy()
    s=a.sum(0);point=float(s[2]/np.sqrt(s[0]*s[1]))
    rng=np.random.default_rng(seed);s=a[rng.integers(0,len(a),size=(2000,len(a)))].sum(1)
    boot=s[:,2]/np.sqrt(s[:,0]*s[:,1]);ci=np.quantile(boot,[.025,.975]).tolist()
    return {'products':int(d.pair_ID.nunique()),'guides':int(d.sgrna.nunique()),'guide_cell_pairs':len(d),'correlation':point,'nominal95_cluster_CI':ci}
identity=read('S03');exact=flag(identity.exact_complete_reference_coding_chain)&flag(identity.exact_complete_mature_coding_DNA)
check('274_relationships_127_exact',len(identity)==274 and int(exact.sum())==127)
main=conditional(read('S12'),'old_LFC','new_LFC',20261015)
close('primary_conditional_r',main['correlation'],expected['within_product_cell_residual_correlation'])
close('primary_conditional_CI',main['nominal95_cluster_CI'],expected['nominal95_cluster_CI'])
follow=conditional(read('S23'),'secondary_LFC','newer_LFC',20261016)
close('additional_assay_r',follow['correlation'],expected['additional_followup_correlation'])
close('additional_assay_CI',follow['nominal95_cluster_CI'],expected['additional_followup_nominal95_cluster_CI'])
fp=E/'cycle_92_secondary_screen_guide_validation/FOLLOWUP_GUIDE_RESULTS.json'
f=json.loads(fp.read_text());f['primary'].update(follow);fp.write_text(json.dumps(f,indent=2)+'\n')

# A fresh reconstruction of the earlier-only selection rule, not stored ranks.
g=read('S25');old=g.groupby(['pair_ID','sgrna']).old_LFC.mean().reset_index().sort_values(['pair_ID','old_LFC','sgrna'],kind='stable')
selected=set(map(tuple,old.groupby('pair_ID',sort=True).head(2)[['pair_ID','sgrna']].to_numpy()))
g['chosen']=[(p,q) in selected for p,q in zip(g.pair_ID,g.sgrna)]
allmean=g.groupby(['pair_ID','cell']).new_LFC.mean();selectedmean=g[g.chosen].groupby(['pair_ID','cell']).new_LFC.mean()
delta=(selectedmean-allmean).groupby('pair_ID').mean()
close('selection_mean',delta.mean(),expected['transferred_newer_LFC_shift'])
check('selection_negative_products',int((delta<0).sum())==79)
clusters=g[['pair_ID','dependence_cluster']].drop_duplicates().set_index('pair_ID').loc[delta.index].dependence_cluster
values=pd.DataFrame({'delta':delta,'cluster':clusters}).sort_values('cluster').delta.to_numpy()
check('selection_one_product_per_cluster',clusters.is_unique)
rng=np.random.default_rng(20261018);boot=values[rng.integers(0,len(values),size=(2000,len(values)))].mean(1)
selection_ci=np.quantile(boot,[.025,.975]).tolist();close('selection_CI',selection_ci,expected['selection_nominal95_cluster_CI'])

# Balanced uncertainty: summarize all released plans/draws; full mode regenerates them.
plans=pd.read_csv(E/'cycle_87_fresh_guide_cross_screen_nORF_validation/all1000_balanced_guide_count_draw_statistics.csv',float_precision='round_trip')
balanced_boot=read('S10');gap=float(plans.reused_minus_fresh_rank_agreement.mean())
gap_ci=np.quantile(balanced_boot.reused_minus_fresh_rank_agreement,[.0125,.9875]).tolist()
close('balanced_point',gap,expected['count_balanced_product_agreement_gap'])
close('balanced_interval',gap_ci,expected['balanced_gap_CI97_5'])

third=read('S29');point=float(spearmanr(third.third_A375_LFC_mean2,third.newer_A375_LFC_mean).statistic)
groups=[a.index.to_numpy() for _,a in third.groupby('dependence_cluster',sort=True)]
rng=np.random.default_rng(20261021);draws=rng.integers(0,len(groups),size=(2000,len(groups)));boots=[]
for draw in draws:
    ix=np.concatenate([groups[k] for k in draw]);q=third.iloc[ix]
    boots.append(float(spearmanr(q.third_A375_LFC_mean2,q.newer_A375_LFC_mean).statistic))
ci=np.quantile(boots,[.025,.975]).tolist()
close('third_study_rank',point,expected['third_study_A375_Spearman_rho'])
close('third_study_CI',ci,expected['third_study_nominal95_cluster_CI'])
check('no_strong_third_study_hits',int(third.third_A375_LFC_mean2.le(-1).sum())==0)
check('196_products_192_groups',len(third)==196 and len(groups)==192)

result={'status':'pass','scope':'Frozen-input statistics; no new genome scan or original-workbook reconstruction',
 'main':main,'additional_assay':follow,'selection':{'products':len(delta),'mean':float(delta.mean()),'CI95':selection_ci},
 'balanced':{'point':gap,'CI97_5':gap_ci,'default_uses_released_frozen_plan_statistics':True},
 'third':{'products':len(third),'rho':point,'CI95':ci},'checks':checks}
(M/'FROZEN_INPUT_REPLAY_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='checks'}))
