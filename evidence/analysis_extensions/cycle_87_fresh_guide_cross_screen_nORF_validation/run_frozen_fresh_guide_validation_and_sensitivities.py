from pathlib import Path
import os,sys,json,gzip,hashlib,datetime,traceback
os.environ.setdefault('OMP_NUM_THREADS','2');os.environ.setdefault('OPENBLAS_NUM_THREADS','2');os.environ.setdefault('MKL_NUM_THREADS','2')
import numpy as np,pandas as pd
from scipy.stats import spearmanr
O=Path(__file__).resolve().parent;N=O.parent;C66=N/'cycle_66_independent_nORF_fitness_benchmark';C83=N/'cycle_83_human_nORF_guide_splice_liability';C35=N/'cycle_35_guide_replicates_and_medical_provenance'
STATE=O/'FRESH_GUIDE_ANALYSIS_JOB_STATE.json'
def state(stage,status='running',**kw):
 d={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pid':os.getpid(),'status':status,'stage':stage,**kw};STATE.write_text(json.dumps(d,indent=2));print(json.dumps(d),flush=True)
def clean(x):
 if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
 if isinstance(x,list):return [clean(v) for v in x]
 if isinstance(x,(np.integer,)):return int(x)
 if isinstance(x,(np.bool_,)):return bool(x)
 if isinstance(x,(float,np.floating)):return float(x) if np.isfinite(x) else None
 return x
def read(p,**kw):return pd.read_csv(p,keep_default_na=False,float_precision='round_trip',**kw)
def rho(a,b):
 if len(a)<3 or not np.isfinite(a).all() or not np.isfinite(b).all() or len(np.unique(a))<2 or len(np.unique(b))<2:return np.nan
 return float(spearmanr(a,b).statistic)
def main():
 plan=json.loads((O/'FROZEN_FRESH_GUIDE_REPLICATION_AND_WITHIN_ASSAY_CONTRAST_PROTOCOL.json').read_text());pins=json.loads((O/'PRE_OUTCOME_JOIN_INPUT_AND_SOURCE_PINS.json').read_text())['pins']
 for p in pins:
  b=Path(p['path']).read_bytes();assert len(b)==p['bytes'] and hashlib.sha256(b).hexdigest()==p['sha256']
 cells=['A375','A549','HEPG2'];den=read(O/'all127_exact_product_cohort_eligibility_and_clusters_before_outcome_join.csv').set_index('pair_ID');og=read(O/'all1010_older_exact_product_guide_eligibility_before_outcome_join.csv');ng=read(O/'all720_newer_exact_product_guide_roles_and_eligibility_before_outcome_join.csv')
 truth=lambda s:s.astype(str).str.lower().eq('true')
 for k in ['primary_eligible','unique_only_eligible','at_least2_each_new_group_eligible','expression_qualified_cohort_eligible','source_AA_clean_eligible']:den[k]=truth(den[k])
 for d in [og,ng]:
  for k in ['primary_eligible','unique_only_eligible']:d[k]=truth(d[k])
 state('join_complete_source_guide_outcomes')
 old=read(C66/'all_original_primary_screen_guide_LFC_rows.csv');old=old[old['Construct Barcode'].isin(og.sgrna)&old['Cell Line'].isin(cells)].rename(columns={'Construct Barcode':'sgrna','Cell Line':'cell','LFC_D21_relative_to_pDNA':'LFC'});assert not old[['sgrna','cell']].duplicated().any()
 newer=read(C83/'all215704_original2026_guide_cell_primary_log_abundance_changes.csv');newer=newer[newer.sgrna.isin(ng.sgrna)&newer.cell.isin(cells)];assert not newer[['sgrna','cell']].duplicated().any()
 oldj=og.merge(old[['sgrna','cell','LFC','source_Excel_row','Gene Symbol']],on='sgrna',how='left',validate='one_to_many');newj=ng.merge(newer.drop(columns=['plasmid_coverage_pass']),on='sgrna',how='left',validate='one_to_many');assert len(oldj)==3030 and len(newj)==2160
 for d in [oldj,newj]:d['LFC']=pd.to_numeric(d.LFC,errors='coerce');d['source_LFC_finite']=np.isfinite(d.LFC)
 oldj.to_csv(O/'all3030_older_exact_product_source_guide_cell_outcomes.csv',index=False);newj.to_csv(O/'all2160_newer_exact_product_source_guide_cell_outcomes_and_roles.csv',index=False)
 counts=read(C35/'all_original_table13_B_guide_measurements.csv');assert counts.sgrna.is_unique;counts=counts.set_index('sgrna').loc[ng.sgrna].copy();pdna=counts[['pDNA_1','pDNA_2','pDNA_3']].to_numpy(float);alternative=[];views={}
 for pc in [1.,.5,5.]:
  view={}
  for cell in cells:
   cols=[cell+'-311Cas9_Rep '+x for x in ['A','B','C']];late=counts[cols].to_numpy(float);effect=np.log2(late+pc).mean(axis=1)-np.log2(pdna+pc).mean(axis=1)
   for q,value in zip(counts.index,effect):view[(q,cell)]=float(value)
   if pc!=1:
    alternative.extend({'sgrna':q,'cell':cell,'pseudocount':pc,'LFC':float(value)} for q,value in zip(counts.index,effect))
  views[pc]=view
 check=np.array([views[1.][(q,c)] for q,c in zip(newj.sgrna,newj.cell)]);assert np.allclose(check,newj.LFC.to_numpy(),atol=1e-11,rtol=1e-11);pd.DataFrame(alternative).to_csv(O/'all4320_newer_guide_cell_alternative_pseudocount_effects.csv',index=False)
 baseindex=pd.MultiIndex.from_product([den.index,cells],names=['pair_ID','cell'])
 def aggregate(flag='primary_eligible',pc=1.):
  a=oldj[oldj[flag]&oldj.source_LFC_finite];b=newj[newj[flag]&newj.source_LFC_finite&newj.guide_group.isin(['reused','fresh'])].copy()
  if pc!=1.:b['LFC']=[views[pc][(q,c)] for q,c in zip(b.sgrna,b.cell)]
  result=pd.DataFrame(index=baseindex);result=result.join(a.groupby(['pair_ID','cell']).LFC.agg(old_mean='mean',old_median='median',old_guides='size'))
  for group in ['reused','fresh']:
   z=b[b.guide_group.eq(group)].groupby(['pair_ID','cell']).agg(LFC=('LFC','mean'),guides=('LFC','size'),GC=('GC_fraction_20nt','mean'));z.columns=[group+'_'+x for x in z.columns];result=result.join(z)
  result['finite_complete_cell']=np.isfinite(result[['old_median','old_mean','reused_LFC','fresh_LFC']]).all(axis=1)&result[['old_guides','reused_guides','fresh_guides']].ge(1).all(axis=1);result['fresh_minus_reused']=result.fresh_LFC-result.reused_LFC;return result
 primary=aggregate();unique=aggregate('unique_only_eligible');primary.to_csv(O/'all381_primary_product_cell_guide_group_statistics.csv');unique.to_csv(O/'all381_unique_only_product_cell_guide_group_statistics.csv')
 finite=primary.finite_complete_cell.groupby('pair_ID').all().reindex(den.index);finiteu=unique.finite_complete_cell.groupby('pair_ID').all().reindex(den.index);den['all3_primary_cell_outcomes_finite']=finite;den['all3_unique_only_cell_outcomes_finite']=finiteu;den['primary_outcome_eligible']=den.primary_eligible&finite;den['unique_only_outcome_eligible']=den.unique_only_eligible&finiteu
 den.to_csv(O/'all127_complete_product_outcome_eligibility_and_denominators.csv')
 def arrays(frame,ids,oldmetric='old_median'):
  take=[(i,c) for i in ids for c in cells];z=frame.loc[take];return z[oldmetric].to_numpy().reshape(len(ids),3),z.reused_LFC.to_numpy().reshape(len(ids),3),z.fresh_LFC.to_numpy().reshape(len(ids),3)
 def stats(A,R,F,retain=[0,1,2]):
  delta=float(np.mean(F[:,retain]-R[:,retain]));rR=[rho(A[:,k],R[:,k]) for k in retain];rF=[rho(A[:,k],F[:,k]) for k in retain];gap=float(np.mean(np.array(rR)-np.array(rF)));return delta,gap,rR,rF
 pids=den.index[den.primary_outcome_eligible].tolist();uids=den.index[den.unique_only_outcome_eligible].tolist();A,R,F=arrays(primary,pids);point=stats(A,R,F);summary=[];cellrows=[]
 variants=[('primary',primary,pids,'old_median'),('older_mean',primary,pids,'old_mean'),('unique_only_no_CDS_exclusion',unique,uids,'old_median'),('prior_source_AA_clean',primary,[i for i in pids if bool(den.loc[i,'source_AA_clean_eligible'])],'old_median'),('expression_qualified_all3',primary,[i for i in pids if bool(den.loc[i,'expression_qualified_cohort_eligible'])],'old_median'),('at_least2_each_new_group',primary,[i for i in pids if bool(den.loc[i,'at_least2_each_new_group_eligible'])],'old_median'),('newer_pc0_5',aggregate(pc=.5),pids,'old_median'),('newer_pc5',aggregate(pc=5.),pids,'old_median')]
 state('declared_point_estimates',primary_products=len(pids),primary_clusters=int(den.loc[pids,'dependence_cluster'].nunique()))
 for name,frame,ids,oldmetric in variants:
  aa,rr,ff=arrays(frame,ids,oldmetric);d,g,rR,rF=stats(aa,rr,ff);summary.append({'variant':name,'products':len(ids),'clusters':int(den.loc[ids,'dependence_cluster'].nunique()),'fresh_minus_reused_LFC':d,'reused_minus_fresh_rank_agreement':g,'old_aggregation':oldmetric})
  for cell,x,y in zip(cells,rR,rF):cellrows.append({'variant':name,'cell':cell,'products':len(ids),'rho_old_vs_reused':x,'rho_old_vs_fresh':y,'rho_difference':x-y})
 pd.DataFrame(summary).to_csv(O/'all_declared_primary_and_sensitivity_point_estimates.csv',index=False);pd.DataFrame(cellrows).to_csv(O/'all_declared_cell_specific_cross_screen_rank_correlations.csv',index=False)
 clusterlabels=den.loc[pids,'dependence_cluster'].to_numpy();clusters=sorted(set(clusterlabels));members=[np.flatnonzero(clusterlabels==c) for c in clusters];rng=np.random.default_rng(20261015);draws=rng.integers(0,len(clusters),size=(2000,len(clusters)));boot=[]
 with gzip.open(O/'all2000_primary_cluster_resampling_plans.json.gz','wt') as out:json.dump({'seed':20261015,'clusters':clusters,'primary_pair_IDs':pids,'draws':draws.tolist()},out)
 for j,draw in enumerate(draws):
  if j%500==0:state('primary_cluster_bootstrap',completed_draws=j,total_draws=2000)
  ix=np.concatenate([members[int(v)] for v in draw]);d,g,rR,rF=stats(A[ix],R[ix],F[ix]);row={'draw':j,'sampled_products':len(ix),'fresh_minus_reused_LFC':d,'reused_minus_fresh_rank_agreement':g}
  for cell,x,y in zip(cells,rR,rF):row[cell+'_rho_old_vs_reused']=x;row[cell+'_rho_old_vs_fresh']=y
  boot.append(row)
 bootdf=pd.DataFrame(boot);bootdf.to_csv(O/'all2000_primary_cluster_bootstrap_statistics.csv',index=False);gate=len(pids)>=10 and len(clusters)>=10 and np.isfinite(point[1]);intervals={}
 for k in ['fresh_minus_reused_LFC','reused_minus_fresh_rank_agreement']:
  vals=bootdf[k].to_numpy();ok=np.isfinite(vals);missing=int((~ok).sum());ci=np.quantile(vals[ok],[.0125,.9875]).tolist() if gate and missing/len(vals)<=.01 else None;intervals[k]={'CI97_5':ci,'undefined_draws':missing,'inference_gate_met':bool(gate and missing/len(vals)<=.01)}
 loo=[]
 for j,cell in enumerate(cells):
  d,g,rR,rF=stats(A,R,F,[k for k in range(3) if k!=j]);loo.append({'omitted_cell':cell,'products':len(pids),'fresh_minus_reused_LFC':d,'reused_minus_fresh_rank_agreement':g})
 pd.DataFrame(loo).to_csv(O/'all3_leave_one_common_cell_out_point_estimates.csv',index=False)
 # Prespecified descriptive adjustment at equal GC and equal group size.
 z=primary.loc[[(i,c) for i in pids for c in cells]].groupby('pair_ID').mean(numeric_only=True).reindex(pids);X=np.column_stack([np.ones(len(pids)),z.fresh_GC-z.reused_GC,np.log2(z.fresh_guides/z.reused_guides)]);rank=int(np.linalg.matrix_rank(X));coef=np.linalg.lstsq(X,z.fresh_minus_reused.to_numpy(),rcond=None)[0] if rank==3 else None;adjust={'design_rank':rank,'columns':['intercept_equal_GC_equal_counts','fresh_minus_reused_GC','log2_fresh_over_reused_guidecount'],'coefficients':coef.tolist() if coef is not None else None,'descriptive_not_causal_or_additional_significance_test':True};(O/'DESCRIPTIVE_GC_AND_GUIDECOUNT_ADJUSTED_WITHIN_ASSAY_DIFFERENCE.json').write_text(json.dumps(adjust,indent=2))
 state('balanced_guide_count_resampling',draws=1000)
 vectors=[];guideinfo=[]
 for i in pids:
  q=newj[newj.pair_ID.eq(i)&newj.primary_eligible&newj.source_LFC_finite]
  qr=sorted(q.loc[q.guide_group.eq('reused'),'sgrna'].unique());qf=sorted(q.loc[q.guide_group.eq('fresh'),'sgrna'].unique());br=np.array([[views[1.][(s,c)] for c in cells] for s in qr]);bf=np.array([[views[1.][(s,c)] for c in cells] for s in qf]);k=min(len(qr),len(qf));assert k>=1;vectors.append((br,bf,k));guideinfo.append({'pair_ID':i,'reused_guides':qr,'fresh_guides':qf,'balanced_count_each':k})
 balanced=[];choices=[]
 for j in range(1000):
  rs=[];fs=[];one=[]
  for br,bf,k in vectors:
   ir=rng.choice(len(br),size=k,replace=False);ifresh=rng.choice(len(bf),size=k,replace=False);rs.append(br[ir].mean(axis=0));fs.append(bf[ifresh].mean(axis=0));one.append([ir.tolist(),ifresh.tolist()])
  d,g,rR,rF=stats(A,np.array(rs),np.array(fs));row={'draw':j,'fresh_minus_reused_LFC':d,'reused_minus_fresh_rank_agreement':g}
  for cell,x,y in zip(cells,rR,rF):row[cell+'_rho_old_vs_reused']=x;row[cell+'_rho_old_vs_fresh']=y
  balanced.append(row);choices.append(one)
 bdf=pd.DataFrame(balanced);bdf.to_csv(O/'all1000_balanced_guide_count_draw_statistics.csv',index=False)
 with gzip.open(O/'all1000_balanced_guide_selection_plans.json.gz','wt') as out:json.dump({'seed':20261015,'RNG_order':'2000clusterplans_then1000balancedguidedraws','guide_inventory':guideinfo,'draw_choices':choices},out)
 bsummary={}
 for k in ['fresh_minus_reused_LFC','reused_minus_fresh_rank_agreement']:
  vals=bdf[k].to_numpy();ok=np.isfinite(vals);v=vals[ok];bsummary[k]={'mean':float(v.mean()),'guide_subset_percentile2_5_97_5':np.quantile(v,[.025,.975]).tolist(),'Monte_Carlo_SE_of_mean':float(v.std(ddof=1)/np.sqrt(len(v))),'undefined_draws':int((~ok).sum()),'not_biological_confidence_interval':True}
 result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'primary_and_all_declared_sensitivities_complete_pending_independent_verification','primary_products':len(pids),'primary_dependence_clusters':len(clusters),'source_shared_cells':cells,'primary_fresh_minus_reused_LFC':point[0],'primary_reused_minus_fresh_rank_agreement':point[1],'primary_intervals':intervals,'primary_per_cell_correlations':[x for x in cellrows if x['variant']=='primary'],'all_declared_sensitivity_points':summary,'all3_cell_omissions':loo,'descriptive_adjusted_delta':adjust,'guide_count_balanced_summary':bsummary,'source_guide_counts_newer_in_primary_cohort':{'reused':int(z.reused_guides.sum()),'fresh':int(z.fresh_guides.sum())},'all_source_original_outcomes_retained':True,'not_Chronos_CNV_or_clinical_endpoint':True,'protein_specific_mechanism_or_new_medical_target_accepted':False}
 (O/'FROZEN_FRESH_GUIDE_VALIDATION_RESULTS.json').write_text(json.dumps(clean(result),indent=2,allow_nan=False));state('all_frozen_analyses_complete_pending_independent_verification','complete',primary_products=len(pids));print(json.dumps(clean(result),allow_nan=False),flush=True)
if __name__=='__main__':
 try:main()
 except BaseException as e:
  (O/'FRESH_GUIDE_ANALYSIS_FAILURE.txt').write_text(traceback.format_exc(),encoding='utf8');state('exception','failed',error=repr(e));raise
