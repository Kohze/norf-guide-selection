from pathlib import Path
import datetime as dt
import hashlib
import json
import math
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
N = ROOT / 'exploration_2026-10-11/novelty_discovery'

def rho(a,b):
    return float(spearmanr(a,b).statistic) if len(a)>=3 else None

def main():
    started = dt.datetime.now(dt.timezone.utc).isoformat()
    protocol_path = OUT / 'FROZEN_PRODUCT_CONCORDANCE_PROTOCOL.json'
    protocol = json.loads(protocol_path.read_text())
    cp = OUT / 'all_metadata_eligibility_and_dependence_before_outcomes.csv'
    assert hashlib.sha256(cp.read_bytes()).hexdigest() == protocol['cohort_sha256']
    cohort = pd.read_csv(cp, dtype={'source_ORF_ID':str})
    cohort = cohort[cohort.metadata_eligible].copy()
    assert len(cohort)==196 and cohort.dependence_cluster.nunique()==192
    cols = ['new_id'] + protocol['third_outcomes']
    third = pd.read_csv(OUT / 'sources/gene_summaries_final_v4.csv',usecols=cols,dtype=str,keep_default_na=False,encoding='latin1')
    invalid = []
    for c in cols[1:]:
        num = pd.to_numeric(third[c],errors='coerce')
        bad = ~np.isfinite(num)
        for i in third.index[bad]:
            invalid.append({'third_new_id':third.loc[i,'new_id'],'column':c,'original_value':third.loc[i,c]})
        third[c] = num
    pd.DataFrame(invalid,columns=['third_new_id','column','original_value']).to_csv(OUT / 'all_third_study_nonfinite_outcomes.csv',index=False)
    third['third_A375_LFC_mean2'] = third[cols[1:]].mean(axis=1,skipna=False)
    guides = pd.read_csv(OUT / 'all_qualified_source2026_guides_before_outcomes.csv')
    measurements = pd.read_csv(N / 'cycle_83_human_nORF_guide_splice_liability/all215704_original2026_guide_cell_primary_log_abundance_changes.csv')
    a375 = measurements[measurements.cell.eq('A375')].copy()
    assert not a375.sgrna.duplicated().any()
    joined = guides.merge(a375,on='sgrna',how='left',validate='many_to_one')
    joined.to_csv(OUT / 'all1147_qualified_source2026_A375_guide_outcomes.csv',index=False)
    assert np.isfinite(joined.LFC).all(), 'Do not silently skip missing source2026 guide results'
    stats = joined.groupby('source_target').agg(newer_A375_LFC_mean=('LFC','mean'),newer_A375_LFC_median=('LFC','median'),newer_guides_observed=('sgrna','nunique'))
    safe = joined[joined.splice_safe_gt50].groupby('source_target').agg(newer_safe_LFC=('LFC','mean'),newer_safe_guides=('sgrna','nunique'))
    paired = cohort.merge(third,left_on='third_new_id',right_on='new_id',how='left',validate='one_to_one').merge(stats,left_on='source_ORF_ID',right_index=True,how='left',validate='one_to_one').merge(safe,left_on='source_ORF_ID',right_index=True,how='left',validate='one_to_one')
    paired['all_outcomes_finite'] = np.isfinite(paired[['third_A375_LFC_mean2','newer_A375_LFC_mean']]).all(axis=1)
    paired['in_exact2021_2026_inventory'] = paired.source_ORF_ID.isin(pd.read_csv(N / 'cycle_87_fresh_guide_cross_screen_nORF_validation/all127_exact_product_cohort_eligibility_and_clusters_before_outcome_join.csv').newer_source_ORF_ID)
    paired = paired.sort_values('source_ORF_ID').reset_index(drop=True)
    paired.to_csv(OUT / 'all196_metadata_eligible_product_outcomes.csv',index=False)
    primary = paired[paired.all_outcomes_finite].copy().reset_index(drop=True)
    assert len(primary)>=protocol['minimum_products']
    x = primary.newer_A375_LFC_mean.to_numpy()
    y = primary.third_A375_LFC_mean2.to_numpy()
    clusters = sorted(primary.dependence_cluster.unique())
    index_groups = [np.flatnonzero(primary.dependence_cluster.eq(c)) for c in clusters]
    rng = np.random.default_rng(20261021)
    draws = rng.integers(0,len(clusters),size=(2000,len(clusters)))
    np.savez_compressed(OUT / 'all2000_cluster_bootstrap_draws.npz',draws=draws,clusters=np.array(clusters))
    bootstrap = []
    for i,d in enumerate(draws):
        rows = np.concatenate([index_groups[j] for j in d])
        bootstrap.append({'draw':i,'products_in_draw':len(rows),'rho':rho(x[rows],y[rows])})
    boots = pd.DataFrame(bootstrap)
    boots.to_csv(OUT / 'all2000_cluster_bootstrap_results.csv',index=False)
    ci = np.quantile(boots.rho,[.025,.975]).tolist()
    rows = []
    def describe(name,frame,xcol='newer_A375_LFC_mean',ycol='third_A375_LFC_mean2'):
        good = frame[np.isfinite(frame[[xcol,ycol]]).all(axis=1)]
        rows.append({'analysis':name,'products':len(good),'clusters':good.dependence_cluster.nunique(),'rho':rho(good[xcol],good[ycol])})
    for c in cols[1:]:
        describe(c,primary,ycol=c)
    describe('within_third_study_replicates',primary,xcol=cols[1],ycol=cols[2])
    describe('source2026_median_guides',primary,xcol='newer_A375_LFC_median')
    describe('source2026_splice_safe_min4',primary[primary.newer_safe_guides.ge(4)],xcol='newer_safe_LFC')
    describe('one_product_per_dependence_cluster',primary.sort_values('source_ORF_ID').drop_duplicates('dependence_cluster'))
    describe('in_exact2021_2026_inventory',primary[primary.in_exact2021_2026_inventory])
    pd.DataFrame(rows).to_csv(OUT / 'all_declared_descriptive_sensitivities.csv',index=False)
    omissions = []
    for c in clusters:
        keep = primary.dependence_cluster.ne(c)
        omissions.append({'omitted_cluster':c,'remaining_products':int(keep.sum()),'rho':rho(x[keep],y[keep])})
    omit = pd.DataFrame(omissions)
    omit.to_csv(OUT / 'all_single_cluster_omissions.csv',index=False)
    k = math.ceil(.2*len(primary))
    topnew = set(primary.sort_values(['newer_A375_LFC_mean','source_ORF_ID']).head(k).source_ORF_ID)
    topthird = set(primary.sort_values(['third_A375_LFC_mean2','source_ORF_ID']).head(k).source_ORF_ID)
    primary['top20percent_source2026'] = primary.source_ORF_ID.isin(topnew)
    primary['top20percent_third'] = primary.source_ORF_ID.isin(topthird)
    primary.to_csv(OUT / 'all_primary_product_pairs_and_rank_overlap.csv',index=False)
    result = {'utc_started':started,'utc_completed':dt.datetime.now(dt.timezone.utc).isoformat(),
              'protocol_sha256':hashlib.sha256(protocol_path.read_bytes()).hexdigest(),
              'status':'completed_pending_separate_verification',
              'primary': {'products':len(primary),'clusters':len(clusters),'newer_guides':len(joined),
                          'rho':rho(x,y),'nominal95_cluster_CI':ci,'pearson':float(pearsonr(x,y).statistic)},
              'missing_products':len(paired)-len(primary),'source2025_nonfinite_cells_full_library':len(invalid),
              'secondary_descriptive':rows,
              'single_cluster_omission_rho_range':[float(omit.rho.min()),float(omit.rho.max())],
              'top20percent':{'selected_in_each':k,'overlap':len(topnew & topthird),'overlap_fraction':len(topnew & topthird)/k,'uniform_random_expectation':k*k/len(primary),'descriptive_only':True},
              'source_guides_independent':None,'guide_level_validation':False,'new_medical_target':False,
              'scope':'Third-study product-level A375 concordance; guide reuse unknown. No peptide-specific or disease-selectivity conclusion.',
              'limits':protocol['claim_limits']}
    (OUT / 'PRODUCT_CONCORDANCE_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__ == '__main__':
    main()
