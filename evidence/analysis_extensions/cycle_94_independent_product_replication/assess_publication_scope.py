from pathlib import Path
import datetime as dt
import json
import pandas as pd

C=Path(__file__).resolve().parent
p=pd.read_csv(C/'all_primary_product_pairs_and_rank_overlap.csv')
scope={
 'utc':dt.datetime.now(dt.timezone.utc).isoformat(),
 'status':'verified_empirical_extension_not_medical_discovery',
 'posthoc_descriptive_check':True,
 'reason':'A positive rank correlation and a relative top20percent overlap must not be mistaken for replication of strong depletion hits.',
 'third_study_reported_hit_threshold':{'depletion_LFC_at_or_below':-1,'source':'Schlesinger etal2025 Figure1E and CRISPR-screen-analysis description; DOI10.1016/j.isci.2025.111884; inspected search-index primary-source text before live-page access challenge.'},
 'qualified_cohort':{'products':len(p),'third_mean_LFC_min':float(p.third_A375_LFC_mean2.min()),'third_mean_LFC_max':float(p.third_A375_LFC_mean2.max()),'third_mean_LFC_at_or_below_minus1':int(p.third_A375_LFC_mean2.le(-1).sum())},
 'novelty_scope':'New exact-product cross-study rank estimate and uncertainty for this archived comparison. No new mechanism, first-ever concordance metric, global priority claim or validated target.',
 'bounded_primary_prior_check':[
  {'url':'https://www.nature.com/articles/s41586-026-10459-x','finding':'Methods meta-analysis integrates references7 and28 plus its current screen; stored Table13 individual-score columns include older Prensner, current Valls-labelled and Hofman assays, and no Schlesinger field. Live primary page search found no Schlesinger mention. This is not an exhaustive literature/systematic priority audit.'},
  {'url':'https://doi.org/10.1016/j.isci.2025.111884','finding':'The source study already compared its ORF overlap/hits with Chen2020 and Prensner2021. That comparison is prior work, not our discovery. Our source2026 comparison is distinct.'}
 ],
 'manuscript_value':'Moderate addition: external study-level aggregate support and an explicit boundary on relative rankings. It does not close the guide-independence or medical-mechanism gaps.',
 'reagent_independence_unknown':True,
 'peptide_causality_established':False,
 'therapeutic_window_established':False,
 'top_medical_paper_established':False,
 'do_not_repeat':['Do not request the blocked PMC article or publisher S3/S4 files again, or retrieve mirrors as a bypass.','Do not reinterpret top40 membership as a hit threshold.','Do not loosen the frozen primary cohort to obtain stronger hits.']
}
(C/'PUBLICATION_SCOPE_AND_EFFECT_MAGNITUDE_REVIEW.json').write_text(json.dumps(scope,indent=2)+'\n',encoding='utf-8')
print(json.dumps(scope['qualified_cohort'],indent=2))
