from pathlib import Path
import json,datetime,hashlib
import numpy as np,pandas as pd
from scipy.stats import rankdata
O=Path(__file__).resolve().parent;C87=O.parent/'cycle_87_fresh_guide_cross_screen_nORF_validation'
r=lambda p:pd.read_csv(p,keep_default_na=False,float_precision='round_trip',low_memory=False)
meta=r(O/'all_shared_guide_metadata_eligible_in_both_sources_before_residual_outcomes.csv')
den=r(O/'all127_shared_guide_concordance_metadata_cohorts.csv')
a=r(C87/'all3030_older_exact_product_source_guide_cell_outcomes.csv')
b=r(C87/'all2160_newer_exact_product_source_guide_cell_outcomes_and_roles.csv')
a=a[['pair_ID','sgrna','cell','LFC']].rename(columns={'LFC':'old_LFC'})
b=b[['pair_ID','sgrna','cell','LFC']].rename(columns={'LFC':'new_LFC'})
j=meta.merge(a,on=['pair_ID','sgrna'],validate='one_to_many').merge(b,on=['pair_ID','sgrna','cell'],validate='one_to_one')
assert len(j)==len(meta)*3 and np.isfinite(j[['old_LFC','new_LFC']].to_numpy(float)).all()
j=j.merge(den[['pair_ID','dependence_cluster']],on='pair_ID',validate='many_to_one')
j['both_splice_safe']=j.splice_safe_gt50_old.astype(str).str.lower().eq('true')&j.splice_safe_gt50_new.astype(str).str.lower().eq('true')
cells=['A375','A549','HEPG2']
def cohort(d,minn=3):
 counts=d.groupby('pair_ID').sgrna.nunique();d=d[d.pair_ID.isin(counts[counts>=minn].index)].copy()
 assert not d.duplicated(['pair_ID','sgrna','cell']).any()
 for v in ['old','new']:
  d[v+'_residual']=d[v+'_LFC']-d.groupby(['pair_ID','cell'])[v+'_LFC'].transform('mean')
 d['gene_weight_per_guide_cell']=1/(3*d.groupby('pair_ID').sgrna.transform('nunique'))
 return d
def weighted_corr(x,y,w):
 xm=np.average(x,weights=w);ym=np.average(y,weights=w)
 return float(np.sum(w*(x-xm)*(y-ym))/np.sqrt(np.sum(w*(x-xm)**2)*np.sum(w*(y-ym)**2)))
def point(d,x='old_residual',y='new_residual'):
 return weighted_corr(d[x].to_numpy(),d[y].to_numpy(),d.gene_weight_per_guide_cell.to_numpy())
pri=cohort(j);assert pri.pair_ID.nunique()==105 and pri.dependence_cluster.nunique()==105 and pri.sgrna.nunique()==471
assert (pri.groupby('pair_ID').gene_weight_per_guide_cell.sum().round(12)==1).all()
pri.to_csv(O/'all1413_primary_shared_guide_cell_paired_residuals.csv',index=False)
clusters=sorted(pri.dependence_cluster.unique());stats=[]
for cl,g in pri.groupby('dependence_cluster',sort=True):
 w=g.gene_weight_per_guide_cell.to_numpy();x=g.old_residual.to_numpy();y=g.new_residual.to_numpy()
 stats.append([float(w.sum()),np.sum(w*x*x),np.sum(w*y*y),np.sum(w*x*y)])
stats=np.array(stats);rng=np.random.default_rng(20261015);draw=rng.integers(0,len(clusters),size=(2000,len(clusters)));s=stats[draw].sum(axis=1);rho=s[:,3]/np.sqrt(s[:,1]*s[:,2])
np.savez_compressed(O/'all2000_primary_cluster_draws.npz',draw_indices=draw,clusters=np.array(clusters),seed=20261015)
pd.DataFrame({'draw':np.arange(2000),'within_product_residual_correlation':rho}).to_csv(O/'all2000_primary_cluster_estimates.csv',index=False)
estimate=point(pri);ci=np.quantile(rho,[.025,.975]).tolist()
# Conditional identity permutation: one permutation per guide block, shared over all3cells.
blocks=[]
for pair,g in pri.groupby('pair_ID',sort=True):
 x=g.pivot(index='sgrna',columns='cell',values='old_residual')[cells].to_numpy()
 y=g.pivot(index='sgrna',columns='cell',values='new_residual')[cells].to_numpy()
 blocks.append((pair,x,y,1/(3*len(x))))
prng=np.random.default_rng(20261016);null=[]
norm=np.sqrt(stats[:,1].sum()*stats[:,2].sum());perms=[]
for drawid in range(5000):
 value=0.
 for pair,x,y,w in blocks:value+=float(np.sum(x*y[prng.permutation(len(y))])*w)
 null.append(value/norm)
pd.DataFrame({'draw':np.arange(5000),'within_product_guide_identity_null_correlation':null}).to_csv(O/'all5000_within_product_identity_permutation_estimates.csv',index=False)
pvalue=float((1+np.sum(np.abs(null)>=abs(estimate)))/(len(null)+1))
sens=[]
for c in cells:
 z=pri[pri.cell.eq(c)];sens.append({'analysis':'cell_specific','cell':c,'products':z.pair_ID.nunique(),'guides':z.sgrna.nunique(),'correlation':point(z)})
for minn,safe,name in [(4,False,'min4_shared_guides'),(3,True,'both_screen_splice_safe_min3')]:
 z=cohort(j[j.both_splice_safe] if safe else j,minn);z.to_csv(O/(name+'_residuals.csv'),index=False);sens.append({'analysis':name,'cell':'all3','products':z.pair_ID.nunique(),'guides':z.sgrna.nunique(),'correlation':point(z)})
ranked=pri.copy()
for side in ['old','new']:
 ranked[side+'_rank_residual']=ranked.groupby(['pair_ID','cell'])[side+'_LFC'].transform(lambda a:rankdata(a)-np.mean(rankdata(a)))
sens.append({'analysis':'within_product_cell_guide_ranks','cell':'all3','products':105,'guides':471,'correlation':point(ranked,'old_rank_residual','new_rank_residual')})
# Descriptive shared GC adjustment, separately per cell, after product-by-cell centering.
gc=pri.copy();gc['gc_residual']=gc.GC_fraction_20nt_old-gc.groupby(['pair_ID','cell']).GC_fraction_20nt_old.transform('mean')
for side in ['old','new']:
 gc[side+'_GC_adjusted_residual']=np.nan
 for cell,g in gc.groupby('cell'):
  xx=g.gc_residual.to_numpy();yy=g[side+'_residual'].to_numpy();ww=g.gene_weight_per_guide_cell.to_numpy();slope=np.sum(ww*xx*yy)/np.sum(ww*xx*xx)
  gc.loc[g.index,side+'_GC_adjusted_residual']=yy-slope*xx
sens.append({'analysis':'descriptive_GC_by_cell_adjusted','cell':'all3','products':105,'guides':471,'correlation':point(gc,'old_GC_adjusted_residual','new_GC_adjusted_residual')})
pd.DataFrame(sens).to_csv(O/'all_declared_shared_guide_concordance_sensitivity_points.csv',index=False)
loo=[]
for i,cl in enumerate(clusters):
 z=stats.sum(axis=0)-stats[i];loo.append({'omitted_cluster':cl,'correlation':float(z[3]/np.sqrt(z[1]*z[2]))})
pd.DataFrame(loo).to_csv(O/'all105_single_product_omission_estimates.csv',index=False)
# Exact covariance decomposition, equal products and equal cells. Not a causal variance attribution.
total=pri.copy()
for side in ['old','new']:
 for c,g in total.groupby('cell'):
  total.loc[g.index,side+'_cell_centered']=g[side+'_LFC']-np.average(g[side+'_LFC'],weights=g.gene_weight_per_guide_cell)
 total[side+'_between_component']=total[side+'_cell_centered']-total[side+'_residual']
w=total.gene_weight_per_guide_cell.to_numpy();weight=w.sum()
within=float(np.sum(w*total.old_residual*total.new_residual)/weight)
between=float(np.sum(w*total.old_between_component*total.new_between_component)/weight)
tcov=float(np.sum(w*total.old_cell_centered*total.new_cell_centered)/weight)
assert abs(tcov-within-between)<1e-12
decomp={'within_product_covariance':within,'between_product_covariance':between,'total_cell_centered_covariance':tcov,'within_fraction_of_covariance':within/tcov,'not_variance_explained_or_causal_fraction':True}
total.to_csv(O/'all1413_shared_guide_covariance_decomposition_rows.csv',index=False)
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'frozen_analysis_complete_pending_independent_verification','primary':{'products':105,'clusters':105,'unique_shared_guides':471,'guide_cell_pairs':1413,'equal_product_cell_residual_Pearson':estimate,'cluster_bootstrap_CI95':ci,'bootstrap_draws':2000,'within_product_guide_identity_permutation_two_sided_P':pvalue,'permutation_draws':5000},'secondary':sens,'all105_product_omission_range':[min(x['correlation'] for x in loo),max(x['correlation'] for x in loo)],'covariance_decomposition':decomp,'does_not_change_uncertain_C88_balanced_product_level_gap':True,'peptide_specific_causality_or_medical_target_established':False}
(O/'WITHIN_PRODUCT_SHARED_GUIDE_CONCORDANCE_RESULTS.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result))
