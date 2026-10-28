from pathlib import Path
import json,hashlib,datetime
import numpy as np,pandas as pd
from scipy.stats import ttest_ind
O=Path(__file__).resolve().parent;N=O.parent
R=lambda p,**kw:pd.read_csv(p,keep_default_na=False,float_precision='round_trip',low_memory=False,**kw)
den=R(O/'all145_within_gene_pair_eligibility_before_outcome_join.csv')
meta=R(O/'all1984_multiORF_expression_context_before_outcome_join.csv')
lfcpath=N/'cycle_83_human_nORF_guide_splice_liability'/'all215704_original2026_guide_cell_primary_log_abundance_changes.csv'
lfc=R(lfcpath)
assert not lfc.duplicated(['sgrna','cell']).any()
L=lfc.pivot(index='sgrna',columns='cell',values='LFC')
cells=['A375','A549','A673','CADOES1','HEPG2','JURKAT','RDES','THP1']
assert set(cells)==set(L.columns)
M=meta.set_index(['orf_targeted','cell_line']);assert M.index.is_unique
G=R(N/'cycle_70_ORF_resolved_context_and_host_dependency'/'all17568_source_resolved_context_host_fitness_and_eligibility.csv',usecols=['orf_targeted','cell_line','native_fitness_score']).set_index(['orf_targeted','cell_line'])
guide_fields=['A_guides','B_guides','A_splice_safe_guides','B_splice_safe_guides']
for col in guide_fields:den[col]=den[col].map(json.loads)
out=[]
for p in den.to_dict('records'):
 for cell in cells:
  z={k:p[k] for k in ['pair_ID','source_gene','ORF_A','ORF_B','biotype_A','biotype_B','dependence_cluster','primary_metadata_eligible','min3_metadata_eligible','splice_safe_metadata_eligible','A_guide_count','B_guide_count','A_splice_safe_count','B_splice_safe_count']}
  for role in ['A','B']:
   ori=p['ORF_'+role];v=M.loc[(ori,cell)]
   z[role+'_RNA_TPM']=float(v.rnaseq_tpm);z[role+'_Ribo_TPM']=float(v.riboseq_tpm)
   z[role+'_depletion']=-float(L.loc[p[role+'_guides'],cell].mean()) if p[role+'_guides'] else np.nan
   z[role+'_splice_safe_depletion']=-float(L.loc[p[role+'_splice_safe_guides'],cell].mean()) if p[role+'_splice_safe_guides'] else np.nan
   z[role+'_Chronos_depletion']=-float(G.loc[(ori,cell),'native_fitness_score'])
  assert z['A_RNA_TPM']==z['B_RNA_TPM']
  z['cell']=cell;z['x_log1p_Ribo_difference']=np.log1p(z['A_Ribo_TPM'])-np.log1p(z['B_Ribo_TPM'])
  z['y_depletion_difference']=z['A_depletion']-z['B_depletion']
  z['y_splice_safe_difference']=z['A_splice_safe_depletion']-z['B_splice_safe_depletion']
  z['y_Chronos_difference']=z['A_Chronos_depletion']-z['B_Chronos_depletion']
  out.append(z)
allout=pd.DataFrame(out);allout.to_csv(O/'all1160_within_gene_pair_cell_outcomes.csv',index=False)
def arrays(d,ycol='y_depletion_difference',cellset=cells):
 p=d.drop_duplicates('pair_ID').sort_values('pair_ID').copy()
 x=d.pivot(index='pair_ID',columns='cell',values='x_log1p_Ribo_difference').loc[p.pair_ID,cellset].to_numpy(float)
 y=d.pivot(index='pair_ID',columns='cell',values=ycol).loc[p.pair_ID,cellset].to_numpy(float)
 w=(1/p.groupby('source_gene').pair_ID.transform('count')).to_numpy(float)
 assert np.isfinite(x).all() and np.isfinite(y).all()
 return p,x,y,w
def fit(x,y,w):
 xr=x-x.mean(axis=1,keepdims=True);yr=y-y.mean(axis=1,keepdims=True)
 xr-=np.average(xr,axis=0,weights=w);yr-=np.average(yr,axis=0,weights=w)
 sxx=np.sum(w[:,None]*xr*xr);syy=np.sum(w[:,None]*yr*yr);sxy=np.sum(w[:,None]*xr*yr)
 sd=np.sqrt(sxx/(x.shape[1]*w.sum()))
 total=np.sum(w[:,None]*(y-np.average(y.mean(axis=1),weights=w))**2)
 gain=sxy*sxy/sxx if sxx>0 else 0.
 return dict(beta_raw=sxy/sxx,beta_per_residual_predictor_SD=sxy/sxx*sd,residual_predictor_SD=sd,partial_R2=gain/syy,incremental_total_R2=gain/total,baseline_FE_R2=1-syy/total,full_R2=1-(syy-gain)/total,weighted_residual_outcome_SD=float(np.sqrt(syy/(x.shape[1]*w.sum()))),pairs=len(w),cells=x.shape[1])
primary=allout[allout.primary_metadata_eligible].copy();p,x,y,w=arrays(primary)
assert p.source_gene.nunique()>=10 and p.dependence_cluster.nunique()>=10
point=fit(x,y,w);point.update(genes=p.source_gene.nunique(),clusters=p.dependence_cluster.nunique())
# Dependence-cluster bootstrap. Sufficient statistics preserve all pairs and contexts per draw.
cl=list(sorted(p.dependence_cluster.unique()));stats=[]
xd=x-x.mean(axis=1,keepdims=True);yd=y-y.mean(axis=1,keepdims=True)
for c in cl:
 ix=np.flatnonzero(p.dependence_cluster.to_numpy()==c);cw=w[ix]
 stats.append(np.r_[cw.sum(),np.sum(cw[:,None]*xd[ix]**2),np.sum(cw[:,None]*yd[ix]**2),np.sum(cw[:,None]*xd[ix]*yd[ix]),np.sum(cw[:,None]*xd[ix],axis=0),np.sum(cw[:,None]*yd[ix],axis=0)])
stats=np.array(stats);rng=np.random.default_rng(20261015);draw=rng.integers(0,len(cl),size=(2000,len(cl)))
np.savez_compressed(O/'all2000_cluster_bootstrap_draws.npz',draw_indices=draw,clusters=np.array(cl),seed=20261015)
weights=np.zeros((len(draw),len(cl)),int)
for i in range(len(draw)):weights[i]=np.bincount(draw[i],minlength=len(cl))
S=weights@stats;W=S[:,0];sx=S[:,4:12];sy=S[:,12:20]
xx=S[:,1]-(sx*sx).sum(axis=1)/W;yy=S[:,2]-(sy*sy).sum(axis=1)/W;xy=S[:,3]-(sx*sy).sum(axis=1)/W
boot=pd.DataFrame({'draw':np.arange(2000),'beta_per_residual_predictor_SD':xy/xx*np.sqrt(xx/(8*W)),'partial_R2':xy**2/(xx*yy),'sampled_cluster_count':len(cl),'sampled_pair_weight':W})
boot.to_csv(O/'all2000_cluster_bootstrap_estimates.csv',index=False)
point['beta_95_cluster_percentile_CI']=np.quantile(boot.beta_per_residual_predictor_SD,[.025,.975]).tolist();point['partial_R2_95_cluster_percentile_CI']=np.quantile(boot.partial_R2,[.025,.975]).tolist()
# Store model inputs/residuals for direct replay.
xr=x-x.mean(axis=1,keepdims=True);yr=y-y.mean(axis=1,keepdims=True);xr-=np.average(xr,axis=0,weights=w);yr-=np.average(yr,axis=0,weights=w)
res=[]
for i,pair in enumerate(p.pair_ID):
 for j,cell in enumerate(cells):res.append({'pair_ID':pair,'cell':cell,'source_gene':p.iloc[i].source_gene,'dependence_cluster':p.iloc[i].dependence_cluster,'gene_weight':w[i],'x':x[i,j],'y':y[i,j],'x_residual':xr[i,j],'y_residual':yr[i,j]})
pd.DataFrame(res).to_csv(O/'primary_weighted_FE_model_inputs_and_residuals.csv',index=False)
sens=[]
for c in cells:
 z=fit(x[:,[k for k,v in enumerate(cells) if v!=c]],y[:,[k for k,v in enumerate(cells) if v!=c]],w);z['analysis']='leave_one_cell_out';z['omitted_cell']=c;sens.append(z)
q,cx,cy,cw=arrays(primary,'y_Chronos_difference');z=fit(cx,cy,cw);z['analysis']='same_experiment_source_Chronos';z['omitted_cell']='';sens.append(z)
for flag,ycol,name in [('min3_metadata_eligible','y_depletion_difference','min3_guides_descriptive'),('splice_safe_metadata_eligible','y_splice_safe_difference','splice_safe_min4_descriptive')]:
 q,cx,cy,cw=arrays(allout[allout[flag]],ycol);z=fit(cx,cy,cw);z['analysis']=name;z['omitted_cell']='';z['genes']=q.source_gene.nunique();z['clusters']=q.dependence_cluster.nunique();sens.append(z)
pd.DataFrame(sens).to_csv(O/'all_declared_secondary_point_estimates.csv',index=False)
# Exploratory Ewing cell-context contrast; these cell lines are not patient replicates.
ew=['A673','CADOES1','RDES'];other=[c for c in cells if c not in ew];leads=[]
for pair,g in primary.groupby('pair_ID',sort=True):
 g=g.set_index('cell');a=g.loc[ew,'y_depletion_difference'].to_numpy(float);b=g.loc[other,'y_depletion_difference'].to_numpy(float)
 e=float(a.mean()-b.mean());t=ttest_ind(a,b,equal_var=False)
 row={'pair_ID':pair,'source_gene':g.source_gene.iloc[0],'ORF_A':g.ORF_A.iloc[0],'ORF_B':g.ORF_B.iloc[0],'Ewing_minus_other_pair_contrast':e,'Welch_t':float(t.statistic),'Welch_context_p':float(t.pvalue)}
 role='A' if e>=0 else 'B';partner='B' if role=='A' else 'A'
 row['oriented_product']=g['ORF_'+role].iloc[0];row['partner_product']=g['ORF_'+partner].iloc[0]
 row['oriented_mean_Ewing_depletion']=float(g.loc[ew,role+'_depletion'].mean())
 row['oriented_Ewing_minus_other_own_depletion']=float(g.loc[ew,role+'_depletion'].mean()-g.loc[other,role+'_depletion'].mean())
 row['all3_Ewing_positive_vs_partner']=bool((g.loc[ew,role+'_depletion']>g.loc[ew,partner+'_depletion']).all())
 row['all3_Ewing_RNA_ge10_Ribo_ge5']=bool((g.loc[ew,role+'_RNA_TPM']>=10).all() and (g.loc[ew,role+'_Ribo_TPM']>=5).all())
 row['oriented_min_Ewing_Ribo_TPM']=float(g.loc[ew,role+'_Ribo_TPM'].min())
 row['oriented_min_Ewing_RNA_TPM']=float(g.loc[ew,role+'_RNA_TPM'].min())
 leads.append(row)
leads=pd.DataFrame(leads);idx=np.argsort(leads.Welch_context_p.to_numpy());pv=leads.Welch_context_p.to_numpy()[idx];qv=np.minimum.accumulate((pv*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1];qvals=np.empty(len(pv));qvals[idx]=np.minimum(qv,1);leads['BH_q_all127_pairs']=qvals
leads['statistical_and_effect_lead_filter']=(leads.BH_q_all127_pairs<.05)&(leads.Ewing_minus_other_pair_contrast.abs()>=.5)
leads['all_qualification_before_guide_omission']=leads.statistical_and_effect_lead_filter&(leads.oriented_mean_Ewing_depletion>=.5)&(leads.oriented_Ewing_minus_other_own_depletion>=.5)&leads.all3_Ewing_positive_vs_partner&leads.all3_Ewing_RNA_ge10_Ribo_ge5
om=[]
for lead in leads[leads.all_qualification_before_guide_omission].to_dict('records'):
 pair=den[den.pair_ID.eq(lead['pair_ID'])].iloc[0];sign=1 if lead['oriented_product']==pair.ORF_A else -1
 for role in ['A','B']:
  for sg in pair[role+'_guides']:
   av=-(L.loc[[v for v in pair.A_guides if not(role=='A' and v==sg)],cells].mean(axis=0))
   bv=-(L.loc[[v for v in pair.B_guides if not(role=='B' and v==sg)],cells].mean(axis=0))
   v=av-bv;value=float(sign*(v[ew].mean()-v[other].mean()))
   om.append({'pair_ID':lead['pair_ID'],'oriented_product':lead['oriented_product'],'omitted_role':role,'omitted_sgrna':sg,'oriented_Ewing_pair_contrast':value,'direction_retained':value>0})
pd.DataFrame(om,columns=['pair_ID','oriented_product','omitted_role','omitted_sgrna','oriented_Ewing_pair_contrast','direction_retained']).to_csv(O/'all_medical_lead_single_guide_omissions.csv',index=False)
stable={k:all(r['direction_retained'] for r in om if r['pair_ID']==k) for k in leads[leads.all_qualification_before_guide_omission].pair_ID}
leads['all_single_guide_omissions_retain_direction']=leads.pair_ID.map(stable).fillna(False).astype(bool)
leads['fully_qualified_exploratory_lead']=leads.all_qualification_before_guide_omission&leads.all_single_guide_omissions_retain_direction
leads.sort_values('BH_q_all127_pairs').to_csv(O/'all127_exploratory_Ewing_pair_contrasts_and_qualification.csv',index=False)
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'frozen_analysis_complete','primary':point,'primary_CI_excludes_zero':bool(np.quantile(boot.beta_per_residual_predictor_SD,.025)>0 or np.quantile(boot.beta_per_residual_predictor_SD,.975)<0),'leave_one_cell_out_beta_range':[min(z['beta_per_residual_predictor_SD'] for z in sens if z['analysis']=='leave_one_cell_out'),max(z['beta_per_residual_predictor_SD'] for z in sens if z['analysis']=='leave_one_cell_out')],'secondary_points':sens[-3:],'exploratory_Ewing':{'all_primary_tests':len(leads),'q_lt05':int((leads.BH_q_all127_pairs<.05).sum()),'statistical_and_effect_filter':int(leads.statistical_and_effect_lead_filter.sum()),'fully_qualified_exploratory_leads':leads[leads.fully_qualified_exploratory_lead].to_dict('records'),'minimum_q':float(leads.BH_q_all127_pairs.min())},'claim_limits':['Retrospective reanalysis of published experiments, not independent biological validation.','Within sourcegene contrasts do not prove same RNA molecule or peptide-mediated causality.','Welch/BH screen compares observational cell contexts; not patient-sampled significance.','Novelty of any product or insight remains to be checked against primary literature.']}
(O/'WITHIN_LOCUS_FUNCTIONAL_PARTITION_RESULTS.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result))
