from pathlib import Path
import json,datetime
import pandas as pd,numpy as np
O=Path(__file__).resolve().parent;N=O.parent
r=lambda p:pd.read_csv(p,keep_default_na=False,float_precision='round_trip',low_memory=False)
z=r(O/'all1413_primary_shared_guide_cell_paired_residuals.csv');report=json.loads((O/'WITHIN_PRODUCT_SHARED_GUIDE_CONCORDANCE_RESULTS.json').read_text())
old=r(N/'cycle_66_independent_nORF_fitness_benchmark'/'all_original_primary_screen_guide_LFC_rows.csv')
old=old[old['Cell Line'].isin(['A375','A549','HEPG2'])].set_index(['Construct Barcode','Cell Line'])
counts=r(N/'cycle_35_guide_replicates_and_medical_provenance'/'all_original_table13_B_guide_measurements.csv').set_index('sgrna')
oldv=np.array([float(old.loc[(q,c),'LFC_D21_relative_to_pDNA']) for q,c in zip(z.sgrna,z.cell)])
newv=[]
for q,c in zip(z.sgrna,z.cell):
 row=counts.loc[q];late=row[[c+'-311Cas9_Rep '+k for k in ['A','B','C']]].to_numpy(float);p=row[['pDNA_1','pDNA_2','pDNA_3']].to_numpy(float)
 newv.append(np.mean(np.log2(late+1))-np.mean(np.log2(p+1)))
newv=np.array(newv);assert np.allclose(oldv,z.old_LFC,atol=1e-12,rtol=1e-12);assert np.allclose(newv,z.new_LFC,atol=1e-12,rtol=1e-12)
# Reconstruct product-by-cell residuals from the independently rejoined source values.
src=z[['pair_ID','dependence_cluster','sgrna','cell']].copy();src['a']=oldv;src['b']=newv
for k in ['a','b']:src[k]-=src.groupby(['pair_ID','cell'])[k].transform('mean')
src['w']=1/(3*src.groupby('pair_ID').sgrna.transform('nunique'))
def corr(d):
 cov=np.cov(d.a.to_numpy(),d.b.to_numpy(),aweights=d.w.to_numpy(),ddof=0);return float(cov[0,1]/np.sqrt(cov[0,0]*cov[1,1]))
point=corr(src);assert abs(point-report['primary']['equal_product_cell_residual_Pearson'])<1e-12
plan=np.load(O/'all2000_primary_cluster_draws.npz');boot=r(O/'all2000_primary_cluster_estimates.csv');members=[src[src.dependence_cluster.eq(c)] for c in plan['clusters']]
err=[]
for i,draw in enumerate(plan['draw_indices']):
 d=pd.concat([members[j] for j in draw],ignore_index=True);v=corr(d);err.append(abs(v-boot.loc[i,'within_product_residual_correlation']))
assert max(err)<1e-12
assert np.allclose(np.quantile(boot.within_product_residual_correlation,[.025,.975]),report['primary']['cluster_bootstrap_CI95'],atol=1e-12,rtol=1e-12)
perm=r(O/'all5000_within_product_identity_permutation_estimates.csv');p=(1+(perm.within_product_guide_identity_null_correlation.abs()>=abs(point)).sum())/5001
assert p==report['primary']['within_product_guide_identity_permutation_two_sided_P']
# Directly replay first100 permutation draws using source residuals and explicit paired rows.
rng=np.random.default_rng(20261016);groups=[]
for pair,g in src.groupby('pair_ID',sort=True):
 av=g.pivot(index='sgrna',columns='cell',values='a')[['A375','A549','HEPG2']].to_numpy()
 bv=g.pivot(index='sgrna',columns='cell',values='b')[['A375','A549','HEPG2']].to_numpy()
 groups.append((av,bv))
pe=[]
for i in range(100):
 av=[];bv=[];ww=[]
 for a,b in groups:
  av.extend(a.ravel());bv.extend(b[rng.permutation(len(b))].ravel());ww.extend([1/(3*len(a))]*a.size)
 cv=np.cov(av,bv,aweights=ww,ddof=0);v=float(cv[0,1]/np.sqrt(cv[0,0]*cv[1,1]));pe.append(abs(v-perm.loc[i,'within_product_guide_identity_null_correlation']))
assert max(pe)<1e-12
# Canonical same-gene residual control is infeasible: positive-control guides supply fewer than3 per mapped source gene.
feas=r(O/'canonical_control_shared_library_metadata_feasibility_all_matches.csv');can=feas[feas.source_target_old.str.startswith('UNASSIGNED_SOURCE_CONTROL::')]
cnt=can.groupby('source_target_old').sgrna.nunique()
res={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'pass','direct_original_guide_cell_pairs_verified':len(src),'source_old_max_abs_error':float(np.max(abs(oldv-z.old_LFC))),'source_new_count_reconstruction_max_abs_error':float(np.max(abs(newv-z.new_LFC))),'independent_weighted_covariance_primary':point,'all2000_explicit_row_bootstrap_max_abs_error':max(err),'first100_explicit_pairing_permutation_max_abs_error':max(pe),'guide_coverage_and_coordinate_gate_inherited_C87':True,'canonical_control_feasibility':{'shared_assigned_positive_control_sequences':len(can),'source_control_gene_labels':len(cnt),'genes_with_atleast3shared_guides':int((cnt>=3).sum()),'status':'no_inferential_canonical_class_comparison_possible'},'novelty_and_source_experiment_scope_review_required':True}
(O/'INDEPENDENT_SOURCE_AND_SHARED_GUIDE_CONCORDANCE_VERIFICATION.json').write_text(json.dumps(res,indent=2));print(json.dumps(res))
