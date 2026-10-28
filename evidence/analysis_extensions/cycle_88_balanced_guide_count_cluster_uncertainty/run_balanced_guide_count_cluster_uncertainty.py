from pathlib import Path
import os,json,gzip,csv,datetime,hashlib,time,traceback
os.environ.setdefault('OMP_NUM_THREADS','2')
import numpy as np
from scipy.stats import spearmanr
O=Path(__file__).resolve().parent;C87=O.parent/'cycle_87_fresh_guide_cross_screen_nORF_validation';C35=O.parent/'cycle_35_guide_replicates_and_medical_provenance'
def rows(p):
 with open(p,encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def state(stage,**kw):
 d={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pid':os.getpid(),'stage':stage,**kw};(O/'BALANCED_COUNT_JOB_STATE.json').write_text(json.dumps(d,indent=2));print(json.dumps(d),flush=True)
def main():
 protocol=json.loads((O/'FROZEN_POSTHOC_BALANCED_GUIDE_CLUSTER_UNCERTAINTY_PROTOCOL.json').read_text())
 bp=json.loads(gzip.open(C87/'all2000_primary_cluster_resampling_plans.json.gz','rt').read());gp=json.loads(gzip.open(C87/'all1000_balanced_guide_selection_plans.json.gz','rt').read());ids=bp['primary_pair_IDs'];cells=['A375','A549','HEPG2']
 agg={(r['pair_ID'],r['cell']):r for r in rows(C87/'all381_primary_product_cell_guide_group_statistics.csv')};A=np.array([[float(agg[(i,c)]['old_median']) for c in cells] for i in ids])
 counts={(r['sgrna'],r['cell']):float(r['LFC']) for r in rows(C87/'all2160_newer_exact_product_source_guide_cell_outcomes_and_roles.csv')};vec=[]
 for i,inv in zip(ids,gp['guide_inventory']):
  assert i==inv['pair_ID'];v=[]
  for group in ['reused','fresh']:
   vv=[]
   for q in inv[group+'_guides']:
    vv.append([counts[(q,c)] for c in cells])
   v.append(np.array(vv))
  vec.append(v)
 R=[];F=[]
 for choice in gp['draw_choices']:
  R.append([vec[j][0][ir].mean(axis=0) for j,(ir,ifr) in enumerate(choice)])
  F.append([vec[j][1][ifr].mean(axis=0) for j,(ir,ifr) in enumerate(choice)])
 R=np.array(R);F=np.array(F);np.savez_compressed(O/'all_verified_balanced_guide_plan_vectors.npz',older=A,reused=R,fresh=F)
 def prep(v):
  order=np.argsort(v,axis=1,kind='stable');sortedv=np.take_along_axis(v,order,axis=1);L=v.shape[1];pos=np.broadcast_to(np.arange(L),v.shape)
  starts=np.concatenate([np.ones((len(v),1),bool),sortedv[:,1:]!=sortedv[:,:-1]],axis=1);ends=np.concatenate([sortedv[:,:-1]!=sortedv[:,1:],np.ones((len(v),1),bool)],axis=1)
  start=np.maximum.accumulate(np.where(starts,pos,0),axis=1);end=np.minimum.accumulate(np.where(ends,pos,L-1)[:,::-1],axis=1)[:,::-1]
  return order,np.argsort(order,axis=1),start,end
 def wrank(p,w):
  order,inv,start,end=p;cum=np.cumsum(w[order],axis=1)
  before=np.where(start==0,0,np.take_along_axis(cum,np.maximum(start-1,0),axis=1))
  after=np.take_along_axis(cum,end,axis=1);rs=(before+after+1)/2
  return np.take_along_axis(rs,inv,axis=1)
 pa=[prep(A[:,c][None,:]) for c in range(3)];pr=[prep(R[:,:,c]) for c in range(3)];pf=[prep(F[:,:,c]) for c in range(3)]
 def corrs(w):
  center=(w.sum()+1)/2;cr=[];cf=[]
  for c in range(3):
   old=wrank(pa[c],w)-center;vr=wrank(pr[c],w)-center;vf=wrank(pf[c],w)-center
   vo=np.sum(w*old*old,axis=1);rr=np.sum(w*old*vr,axis=1)/np.sqrt(vo*np.sum(w*vr*vr,axis=1));ff=np.sum(w*old*vf,axis=1)/np.sqrt(vo*np.sum(w*vf*vf,axis=1));cr.append(rr);cf.append(ff)
  return np.array(cr),np.array(cf)
 den={r['pair_ID']:r['dependence_cluster'] for r in rows(C87/'all127_exact_product_cohort_eligibility_and_clusters_before_outcome_join.csv')}
 ixby=[np.array([j for j,i in enumerate(ids) if den[i]==c],int) for c in bp['clusters']]
 # Frozen explicit duplicate validation, first20 cluster draws x first20 guide plans x3cells x2groups.
 checks=0;maxerr=0
 for draw in bp['draws'][:20]:
  ix=np.concatenate([ixby[k] for k in draw]);w=np.bincount(ix,minlength=len(ids));rr,ff=corrs(w)
  for j in range(20):
   for c in range(3):
    for out,v in [(rr,R),(ff,F)]:
     expected=float(spearmanr(A[ix,c],v[j,ix,c]).statistic);err=abs(expected-out[c,j]);assert err<2e-12;maxerr=max(maxerr,err);checks+=1
 state('explicit_duplicate_validation_pass',checks=checks,max_error=maxerr)
 w=np.ones(len(ids),int);rr,ff=corrs(w);perplan=(rr-ff).mean(axis=0);delta=(F-R).mean(axis=(1,2));oldstat=rows(C87/'all1000_balanced_guide_count_draw_statistics.csv')
 for j,r in enumerate(oldstat):
  assert abs(perplan[j]-float(r['reused_minus_fresh_rank_agreement']))<2e-12;assert abs(delta[j]-float(r['fresh_minus_reused_LFC']))<2e-12
 point={'reused_minus_fresh_rank_agreement':float(perplan.mean()),'fresh_minus_reused_LFC':float(delta.mean())};mean_product_delta=(F-R).mean(axis=(0,2));boot=[];start=time.monotonic()
 for j,draw in enumerate(bp['draws']):
  ix=np.concatenate([ixby[k] for k in draw]);w=np.bincount(ix,minlength=len(ids));rr,ff=corrs(w)
  gap=float((rr-ff).mean());d=float(np.dot(mean_product_delta,w)/w.sum());boot.append({'draw':j,'sampled_products':int(w.sum()),'reused_minus_fresh_rank_agreement':gap,'fresh_minus_reused_LFC':d,'undefined_plan_cell_correlations':int((~np.isfinite(rr)).sum()+(~np.isfinite(ff)).sum())})
  if j%100==0:state('balanced_plan_averaged_cluster_bootstrap',completed=j,total=2000,elapsed_seconds=time.monotonic()-start)
 with open(O/'all2000_balanced_plan_averaged_cluster_bootstrap_statistics.csv','w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=list(boot[0]));writer.writeheader();writer.writerows(boot)
 intervals={}
 for key in point:
  v=np.array([r[key] for r in boot]);ok=np.isfinite(v);intervals[key]={'point':point[key],'CI97_5':np.quantile(v[ok],[.0125,.9875]).tolist() if (1-ok.mean())<=.01 else None,'undefined_cluster_draws':int((~ok).sum())}
 inputs=[C87/'all2000_primary_cluster_resampling_plans.json.gz',C87/'all1000_balanced_guide_selection_plans.json.gz',C87/'all381_primary_product_cell_guide_group_statistics.csv',C87/'all2160_newer_exact_product_source_guide_cell_outcomes_and_roles.csv',O/'FROZEN_POSTHOC_BALANCED_GUIDE_CLUSTER_UNCERTAINTY_PROTOCOL.json']
 pins=[{'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in inputs]
 report={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'complete','posthoc_robustness_not_primary_discovery':True,'products':len(ids),'clusters':len(bp['clusters']),'cells':cells,'guide_plans':len(R),'cluster_draws':len(boot),'estimates':intervals,'explicit_duplicate_verification':{'checks':checks,'maximum_absolute_difference':maxerr,'status':'pass'},'all1000_point_plan_statistics_match_C87':True,'all_bootstrap_undefined_plan_cell_correlations':sum(x['undefined_plan_cell_correlations'] for x in boot),'point_MC_diagnostics':{'first500_gap_mean':float(perplan[:500].mean()),'last500_gap_mean':float(perplan[500:].mean()),'gap_MCSE':float(perplan.std(ddof=1)/np.sqrt(len(perplan))),'delta_MCSE':float(delta.std(ddof=1)/np.sqrt(len(delta)))},'input_pins':pins,'elapsed_seconds':time.monotonic()-start}
 (O/'BALANCED_GUIDE_COUNT_CLUSTER_UNCERTAINTY_RESULTS.json').write_text(json.dumps(report,indent=2));state('complete',status='complete',elapsed_seconds=time.monotonic()-start);print(json.dumps(report),flush=True)
if __name__=='__main__':
 try:main()
 except:
  (O/'BALANCED_COUNT_FAILURE.txt').write_text(traceback.format_exc());state('failed',status='failed');raise
