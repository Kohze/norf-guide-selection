from pathlib import Path
import json,datetime,hashlib
import pandas as pd,numpy as np
O=Path(__file__).resolve().parent;N=O.parent;C87=N/'cycle_87_fresh_guide_cross_screen_nORF_validation'
r=lambda p:pd.read_csv(p,keep_default_na=False,float_precision='round_trip',low_memory=False)
old=r(C87/'all1010_older_exact_product_guide_eligibility_before_outcome_join.csv')
new=r(C87/'all720_newer_exact_product_guide_roles_and_eligibility_before_outcome_join.csv')
den=r(C87/'all127_exact_product_cohort_eligibility_and_clusters_before_outcome_join.csv')
truth=lambda s:s.astype(str).str.lower().eq('true')
a=old[truth(old.primary_eligible)].copy();b=new[truth(new.primary_eligible)&new.guide_group.eq('reused')].copy()
join=a.merge(b,on=['pair_ID','sgrna'],suffixes=('_old','_new'),validate='one_to_one')
assert (join.cut_interbase0_old==join.cut_interbase0_new).all()
assert (join.chrom_old==join.chrom_new).all()
assert join.sgrna.is_unique
join.to_csv(O/'all_shared_guide_metadata_eligible_in_both_sources_before_residual_outcomes.csv',index=False)
rows=[]
for p in den.to_dict('records'):
 g=join[join.pair_ID.eq(p['pair_ID'])]
 rows.append({'pair_ID':p['pair_ID'],'old_source_name':p['older_source_name'],'new_source_ORF_ID':p['newer_source_ORF_ID'],'dependence_cluster':p['dependence_cluster'],'shared_both_eligible_guides':len(g),'both_splice_safe_shared_guides':int((truth(g.splice_safe_gt50_old)&truth(g.splice_safe_gt50_new)).sum()),'primary_ge3shared_metadata_eligible':len(g)>=3,'sensitivity_ge4shared_metadata_eligible':len(g)>=4,'sensitivity_splice_safe_ge3shared_metadata_eligible':int((truth(g.splice_safe_gt50_old)&truth(g.splice_safe_gt50_new)).sum())>=3})
d=pd.DataFrame(rows);d.to_csv(O/'all127_shared_guide_concordance_metadata_cohorts.csv',index=False)
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_exact_products':127,'all_both_eligible_shared_guides':len(join),'cohorts':{}}
for flag in ['primary_ge3shared_metadata_eligible','sensitivity_ge4shared_metadata_eligible','sensitivity_splice_safe_ge3shared_metadata_eligible']:
 z=d[d[flag]];result['cohorts'][flag]={'products':len(z),'clusters':z.dependence_cluster.nunique(),'shared_guides':int(z.shared_both_eligible_guides.sum())}
(O/'SHARED_GUIDE_METADATA_FEASIBILITY.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result))
