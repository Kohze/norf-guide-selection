"""Outcome-blind eligibility and dependence components for the third study."""
from pathlib import Path
import datetime as dt
import hashlib
import json
import re
import pandas as pd

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
N = ROOT / 'exploration_2026-10-11/novelty_discovery'

def main():
    matched = pd.read_csv(OUT / 'all_locus_compatible_exact_AA_metadata_before_outcomes.csv', dtype=str, keep_default_na=False)
    models = pd.read_csv(N / 'cycle_86_cross_screen_complete_product_and_reagent_identity/all2196_original2026_complete_sequence_models_before_cross_screen_comparison.csv', dtype=str, keep_default_na=False).set_index('source_ORF_ID')
    geom = pd.read_csv(N / 'cycle_83_human_nORF_guide_splice_liability/all_original2026_library_guide_splice_and_CDS_geometry_before_outcome_join.csv', dtype=str, keep_default_na=False)
    plasmid = pd.read_csv(N / 'cycle_35_guide_replicates_and_medical_provenance/all_original_table13_B_guide_measurements.csv', usecols=['sgrna','pDNA_1','pDNA_2','pDNA_3'])
    assert not plasmid.sgrna.duplicated().any()
    plasmid['coverage_pass'] = plasmid[['pDNA_1','pDNA_2','pDNA_3']].ge(10).all(axis=1)
    geom = geom.merge(plasmid[['sgrna','coverage_pass']], on='sgrna', how='left', validate='many_to_one')
    geom['qualified'] = geom.primary_geometry_eligible.eq('True') & geom.unique_primary_NGG_site.eq('True') & geom.source_distinct_target_assignments_for_guide.eq('1') & geom.coverage_pass.eq(True)
    guides = geom[geom.qualified].drop_duplicates(['source_target','sgrna'])
    guide_counts = guides.groupby('source_target').sgrna.nunique()
    matched['newer_qualified_guides'] = matched.source_ORF_ID.map(guide_counts).fillna(0).astype(int)
    # Test exact coding intervals allowing only a source convention of including
    # or omitting the terminal stop. Save both coordinate interpretations.
    exact0, exact1 = [], []
    for r in matched.to_dict('records'):
        m = models.loc[r['source_ORF_ID']]
        b = [list(map(int, x.split('-'))) for x in m.blocks_in_RNA_order0.split('|')]
        full = sorted(tuple(x) for x in b)
        if m.terminal_inframe_stop == 'True':
            remaining = 3
            while remaining:
                length = b[-1][1]-b[-1][0]
                remove = min(remaining,length)
                if m.strand == '+':
                    b[-1][1] -= remove
                else:
                    b[-1][0] += remove
                remaining -= remove
                if b[-1][0] == b[-1][1]:
                    b.pop()
        coding = sorted(tuple(x) for x in b)
        parsed = [re.fullmatch(r'(chr[^:]+):(\d+)-(\d+):([+-])', p) for p in r['third_loci'].split('|')]
        ok = bool(parsed and all(x and x[1]==m.chrom and x[4]==m.strand for x in parsed))
        third0 = sorted((int(x[2]),int(x[3])) for x in parsed) if ok else []
        third1 = sorted((int(x[2])-1,int(x[3])) for x in parsed) if ok else []
        exact0.append(ok and (third0==coding or third0==full))
        exact1.append(ok and (third1==coding or third1==full))
    matched['exact_blocks_if_third0_halfopen'] = exact0
    matched['exact_blocks_if_third1_inclusive'] = exact1
    # The exact coordinate requirement is additional specificity, never an AA rescue.
    matched['metadata_eligible'] = matched.metadata_one_to_one.eq('True') & matched.third_any_ccds_overlap.eq('FALSE') & pd.to_numeric(matched.third_num_guides).ge(4) & matched.newer_qualified_guides.ge(4) & (matched.exact_blocks_if_third0_halfopen | matched.exact_blocks_if_third1_inclusive)
    # Union products from the same stable gene, overlapping genomic coding bodies,
    # or shared third-study guide IDs. No outcome columns enter these components.
    ids = list(matched.loc[matched.metadata_eligible,'source_ORF_ID'])
    parent = {x:x for x in ids}
    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]
    def union(a,b):
        parent[find(a)] = find(b)
    eligible = matched[matched.metadata_eligible].to_dict('records')
    for i,a in enumerate(eligible):
        for b in eligible[:i]:
            samegene = bool(a['gene_2026'] and a['gene_2026']==b['gene_2026'])
            shared = bool(set(a['third_guide_ids'].split(';')) & set(b['third_guide_ids'].split(';')) - {''})
            overlap = a['chrom_2026']==b['chrom_2026'] and any(max(x,u)<min(y,v) for x,y in (map(int,p.split('-')) for p in a['blocks_2026_0'].split('|')) for u,v in (map(int,p.split('-')) for p in b['blocks_2026_0'].split('|')))
            if samegene or shared or overlap:
                union(a['source_ORF_ID'],b['source_ORF_ID'])
    groups = {}
    for x in ids:
        groups.setdefault(find(x),[]).append(x)
    component = {x: f'cluster{k:03d}' for k,g in enumerate(sorted(groups.values(),key=lambda g:min(g))) for x in g}
    matched['dependence_cluster'] = matched.source_ORF_ID.map(component).fillna('')
    matched.to_csv(OUT / 'all_metadata_eligibility_and_dependence_before_outcomes.csv', index=False)
    guides[guides.source_target.isin(ids)].to_csv(OUT / 'all_qualified_source2026_guides_before_outcomes.csv', index=False)
    report = {'utc': dt.datetime.now(dt.timezone.utc).isoformat(), 'stage': 'metadata_only',
              'metadata_outcomes_opened': False, 'exact_AA_locus_pairs': len(matched),
              'exact0_block_pairs': sum(exact0), 'exact1_block_pairs': sum(exact1),
              'one_to_one_no_CCDS': int((matched.metadata_one_to_one.eq('True') & matched.third_any_ccds_overlap.eq('FALSE')).sum()),
              'eligible_products': int(matched.metadata_eligible.sum()), 'dependence_clusters': len(groups),
              'component_sizes': sorted([len(g) for g in groups.values()],reverse=True),
              'eligible_newer_guides': int(guides.source_target.isin(ids).sum()),
              'guide_reuse_known': False,
              'interpretation': 'Cross-study product concordance only. Independent reagent validation remains unavailable.',
              'cohort_sha256': hashlib.sha256((OUT / 'all_metadata_eligibility_and_dependence_before_outcomes.csv').read_bytes()).hexdigest()}
    (OUT / 'COHORT_PREFLIGHT.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__ == '__main__':
    main()
