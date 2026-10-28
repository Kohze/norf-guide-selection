"""Separate source calculations, custom rank statistic, and interval verification.

This is a computational check by the same assistant, not external peer review.
"""
from pathlib import Path
import csv
import datetime as dt
import hashlib
import json
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
N = ROOT / 'exploration_2026-10-11/novelty_discovery'

def average_ranks(values):
    # Independent implementation: ranks from unique-value counts.
    unique, inverse, counts = np.unique(values,return_inverse=True,return_counts=True)
    ends = np.cumsum(counts)
    middles = (ends-counts+1+ends)/2
    return middles[inverse]

def rank_correlation(a,b):
    a,b = average_ranks(a),average_ranks(b)
    a,b = a-a.mean(),b-b.mean()
    return float(np.dot(a,b)/np.sqrt(np.dot(a,a)*np.dot(b,b)))

def main():
    checks = []
    def check(name,condition):
        assert condition,name
        checks.append(name)
    result = json.loads((OUT / 'PRODUCT_CONCORDANCE_RESULTS.json').read_text())
    protocol_path = OUT / 'FROZEN_PRODUCT_CONCORDANCE_PROTOCOL.json'
    protocol = json.loads(protocol_path.read_text())
    check('protocol_hash_unchanged',hashlib.sha256(protocol_path.read_bytes()).hexdigest()==result['protocol_sha256'])
    check('protocol_before_outcome_analysis',protocol['utc']<result['utc_started'])
    p = pd.read_csv(OUT / 'all_primary_product_pairs_and_rank_overlap.csv')
    check('196_unique_products_192_clusters',len(p)==p.source_ORF_ID.nunique()==196 and p.dependence_cluster.nunique()==192)
    check('all_exact_AA_and_genomic_blocks',p.exact_blocks_if_third0_halfopen.all() and p.metadata_one_to_one.all())
    check('source25_no_CCDS_min4',p.third_any_ccds_overlap.eq(False).all() and p.third_num_guides.ge(4).all())
    check('source26_min4_guides',p.newer_guides_observed.ge(4).all())
    selected = set(p.third_new_id)
    records = {}
    with (OUT / 'sources/gene_summaries_final_v4.csv').open(encoding='latin1',newline='') as f:
        for row in csv.DictReader(f):
            if row['new_id'] in selected:
                check('source25_unique_id_'+row['new_id'],row['new_id'] not in records)
                records[row['new_id']] = (float(row['A375_r1_1_0_lfc'])+float(row['A375_r2_1_0_lfc']))/2
    expected_third = np.array([records[x] for x in p.third_new_id])
    check('source25_all196_raw_CSV_mean2',np.allclose(expected_third,p.third_A375_LFC_mean2,rtol=0,atol=1e-12))
    guides = pd.read_csv(OUT / 'all1147_qualified_source2026_A375_guide_outcomes.csv')
    raw = pd.read_csv(N / 'cycle_35_guide_replicates_and_medical_provenance/all_original_table13_B_guide_measurements.csv')
    raw = raw.set_index('sgrna').loc[guides.sgrna]
    cols_late = ['A375-311Cas9_Rep A','A375-311Cas9_Rep B','A375-311Cas9_Rep C']
    plasmid = raw[['pDNA_1','pDNA_2','pDNA_3']].to_numpy(float)
    late = raw[cols_late].to_numpy(float)
    calc = np.log2(late+1).mean(axis=1)-np.log2(plasmid+1).mean(axis=1)
    check('source26_all1147_direct_processed_count_reconstruction',np.max(np.abs(calc-guides.LFC))<1e-12)
    check('source26_all1147_plasmid_coverage',bool((plasmid>=10).all()))
    direct = {}
    for target,val in zip(guides.source_target,calc):
        direct.setdefault(target,[]).append(val)
    expected_new = np.array([np.mean(direct[x]) for x in p.source_ORF_ID])
    check('source26_all196_direct_product_means',np.max(np.abs(expected_new-p.newer_A375_LFC_mean))<1e-12)
    x,y = expected_new,expected_third
    check('independent_rank_primary',abs(rank_correlation(x,y)-result['primary']['rho'])<1e-12)
    npz = np.load(OUT / 'all2000_cluster_bootstrap_draws.npz')
    groups = [np.flatnonzero(p.dependence_cluster.eq(c)) for c in npz['clusters']]
    computed = []
    for draw in npz['draws']:
        idx = np.concatenate([groups[k] for k in draw])
        computed.append(rank_correlation(x[idx],y[idx]))
    saved = pd.read_csv(OUT / 'all2000_cluster_bootstrap_results.csv')
    check('independent_ranks_all2000_cluster_draws',np.max(np.abs(computed-saved.rho))<1e-12)
    check('nominal_interval_reconstructed',np.max(np.abs(np.quantile(computed,[.025,.975])-result['primary']['nominal95_cluster_CI']))<1e-12)
    k = result['top20percent']['selected_in_each']
    a = set(p.iloc[np.lexsort((p.source_ORF_ID.to_numpy(),x))[:k]].source_ORF_ID)
    b = set(p.iloc[np.lexsort((p.source_ORF_ID.to_numpy(),y))[:k]].source_ORF_ID)
    check('rank_overlap_reconstructed',len(a&b)==result['top20percent']['overlap']==17)
    baseline = json.loads((OUT / 'baseline/CURRENT_NORFS_MANUSCRIPT.json').read_text())
    check('v14_pdf_unchanged',hashlib.sha256(Path(baseline['pdf']).read_bytes()).hexdigest()==baseline['pdf_sha256'])
    check('no_independent_guide_or_medical_claim',result['source_guides_independent'] is None and result['guide_level_validation'] is False and result['new_medical_target'] is False)
    # Compact check names; retain the source uniqueness count without196 labels.
    individual = [s for s in checks if s.startswith('source25_unique_id_')]
    checks = [s for s in checks if not s.startswith('source25_unique_id_')]
    checks.append(f'source25_unique_ids_checked:{len(individual)}')
    sources = [OUT/'sources/gene_summaries_final_v4.csv',protocol_path,
               OUT/'all_metadata_eligibility_and_dependence_before_outcomes.csv',
               N/'cycle_35_guide_replicates_and_medical_provenance/all_original_table13_B_guide_measurements.csv',
               N/'cycle_83_human_nORF_guide_splice_liability/all_original2026_library_guide_splice_and_CDS_geometry_before_outcome_join.csv']
    report = {'utc':dt.datetime.now(dt.timezone.utc).isoformat(),'status':'pass','checks':checks,
              'separate_source_and_statistic_implementation':True,'independent_external_expert_review':False,
              'source_sha256':{str(s.relative_to(ROOT)):hashlib.sha256(s.read_bytes()).hexdigest() for s in sources},
              'primary_rho_verified':rank_correlation(x,y),
              'nominal95_CI_verified':np.quantile(computed,[.025,.975]).tolist()}
    (OUT / 'SEPARATE_COMPUTATIONAL_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'pass','checks':len(checks),'rho':report['primary_rho_verified'],'CI':report['nominal95_CI_verified']},indent=2))

if __name__ == '__main__':
    main()
