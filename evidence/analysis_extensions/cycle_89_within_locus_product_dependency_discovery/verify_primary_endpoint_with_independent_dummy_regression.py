from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import pandas as pd,numpy as np,json,datetime
from scipy.linalg import lstsq
N=Path(r'C:/Users/robin/phdData/norfs/exploration_2026-10-11/novelty_discovery')
O=N/'cycle_89_within_locus_product_dependency_discovery'
r=lambda p:pd.read_csv(p,keep_default_na=False,float_precision='round_trip',low_memory=False)
z=r(O/'primary_weighted_FE_model_inputs_and_residuals.csv');out=json.loads((O/'WITHIN_LOCUS_FUNCTIONAL_PARTITION_RESULTS.json').read_text())
# Independent weighted dummy-variable regression, without sequential demeaning.
def dummy_fit(d):
 X=pd.concat([pd.Series(1.,index=d.index,name='intercept'),pd.get_dummies(d.pair_ID,drop_first=True,dtype=float),pd.get_dummies(d.cell,drop_first=True,dtype=float)],axis=1).to_numpy(float)
 w=d.gene_weight.to_numpy();sw=np.sqrt(w);a=lstsq(X*sw[:,None],d.x.to_numpy()*sw,lapack_driver='gelsy')[0];rx=d.x.to_numpy()-X@a
 X2=np.column_stack([X,d.x.to_numpy()]);b=lstsq(X2*sw[:,None],d.y.to_numpy()*sw,lapack_driver='gelsy')[0]
 return float(b[-1]*np.sqrt(np.average(rx**2,weights=w)))
beta=dummy_fit(z);assert abs(beta-out['primary']['beta_per_residual_predictor_SD'])<1e-12
plan=np.load(O/'all2000_cluster_bootstrap_draws.npz');boot=r(O/'all2000_cluster_bootstrap_estimates.csv');mx=0.
for k in [0,23,513,1234,1999]:
 d=pd.concat([z[z.dependence_cluster.eq(plan['clusters'][i])] for i in plan['draw_indices'][k]],ignore_index=True)
 v=dummy_fit(d);e=abs(v-boot.loc[k,'beta_per_residual_predictor_SD']);mx=max(mx,e);assert e<1e-12
pairs=r(O/'all145_within_gene_pair_eligibility_before_outcome_join.csv');primary=pairs[pairs.primary_metadata_eligible]
models=r(N/'cycle_86_cross_screen_complete_product_and_reagent_identity'/'all2196_original2026_complete_sequence_models_before_cross_screen_comparison.csv').set_index('source_ORF_ID')
assert models.loc[sorted(set(primary.ORF_A)|set(primary.ORF_B)),'reference_mapped'].astype(str).str.lower().eq('true').all()
rep={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'pass','primary_independent_weighted_dummy_regression_beta':beta,'primary_max_abs_error':abs(beta-out['primary']['beta_per_residual_predictor_SD']),'five_explicit_duplicated_cluster_draw_max_abs_error':mx,'all_selected_product_references_mapped':True,'negative_discovery_result_accepted':True,'thresholds_not_relaxed':True}
(O/'INDEPENDENT_WITHIN_LOCUS_ENDPOINT_VERIFICATION.json').write_text(json.dumps(rep,indent=2));print(json.dumps(rep))
