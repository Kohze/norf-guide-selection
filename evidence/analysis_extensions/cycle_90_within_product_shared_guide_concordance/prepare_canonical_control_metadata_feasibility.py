from pathlib import Path
import pandas as pd,json
N=Path(r'C:/Users/robin/phdData/norfs/exploration_2026-10-11/novelty_discovery');O=N/'cycle_90_within_product_shared_guide_concordance'
r=lambda p:pd.read_csv(p,keep_default_na=False,low_memory=False,float_precision='round_trip')
a=r(N/'cycle_83_human_nORF_guide_splice_liability'/'all_original2021_library_guide_splice_and_CDS_geometry_before_outcome_join.csv')
b=r(N/'cycle_83_human_nORF_guide_splice_liability'/'all_original2026_library_guide_splice_and_CDS_geometry_before_outcome_join.csv')
oldmodels=r(N/'cycle_86_cross_screen_complete_product_and_reagent_identity'/'all553_original2021_complete_sequence_models_before_cross_screen_comparison.csv')
print('oldmodelcols',oldmodels.columns.tolist())
m=a.merge(b,on='sgrna',suffixes=('_old','_new'))
print('all_shared',len(m))
print('new parental',m.source_target_new.str.endswith('.parental').sum())
print('parental source names',m[m.source_target_new.str.endswith('.parental')][['source_target_old','source_target_new']].drop_duplicates().head(25).to_string(index=False))
m.to_csv(O/'canonical_control_shared_library_metadata_feasibility_all_matches.csv',index=False)
