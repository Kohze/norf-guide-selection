"""Exact protein/locus matching without opening source fitness columns."""
from pathlib import Path
import datetime as dt
import hashlib
import json
import re
import pandas as pd

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
N = ROOT / 'exploration_2026-10-11/novelty_discovery'

def aa(s):
    # Formatting whitespace and one terminal stop marker only. Internal symbols stay.
    s = re.sub(r'\s+', '', s)
    return s[:-1] if s.endswith('*') else s

def loci(s):
    parts = s.split('|')
    result = []
    for part in parts:
        m = re.fullmatch(r'(chr[^:]+):(\d+)-(\d+):([+-])', part)
        if not m:
            return []
        result.append((m[1], int(m[2])-1, int(m[3]), m[4]))
    return result

def main():
    third = pd.read_csv(OUT / 'all_third_study_identity_and_reagent_metadata.csv', dtype=str, keep_default_na=False)
    new = pd.read_csv(N / 'cycle_86_cross_screen_complete_product_and_reagent_identity/all2196_original2026_complete_sequence_models_before_cross_screen_comparison.csv', dtype=str, keep_default_na=False)
    old = pd.read_csv(N / 'cycle_86_cross_screen_complete_product_and_reagent_identity/all553_original2021_complete_sequence_models_before_cross_screen_comparison.csv', dtype=str, keep_default_na=False)
    lookup = {}
    for r in third.to_dict('records'):
        seq = aa(r['Amino.acid.sequence'])
        if seq and re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+', seq):
            lookup.setdefault(seq, []).append(r)
    rows = []
    for r in new.to_dict('records'):
        seq = aa(r['literal_standard_AA'])
        for t in lookup.get(seq, []):
            parsed = loci(t['loci_all'])
            blocks = [tuple(map(int, x.split('-'))) for x in r['blocks_in_RNA_order0'].split('|')] if r['blocks_in_RNA_order0'] else []
            overlap = sum(max(0, min(b, y)-max(a, x)) for a,b in blocks for chrom,x,y,strand in parsed if chrom==r['chrom'] and strand==r['strand'])
            nt = sum(b-a for a,b in blocks)
            compatible = bool(blocks and parsed and overlap >= .9*min(nt, sum(y-x for chrom,x,y,strand in parsed)))
            rows.append({'source_ORF_ID': r['source_ORF_ID'], 'gene_2026': r['source_single_stable_gene_ID'],
                         'AA_sequence': seq, 'AA_length': len(seq), 'chrom_2026': r['chrom'], 'strand_2026': r['strand'],
                         'blocks_2026_0': r['blocks_in_RNA_order0'],
                         'third_new_id': t['new_id'], 'third_pool_id': t['id'],
                         'third_loci': t['loci_all'], 'third_gene_label': t['gene_id'],
                         'third_guide_ids': t['guides_list'], 'third_num_guides': t['num_guides'],
                         'third_any_ccds_overlap': t['any_ccds_overlap'],
                         'coding_overlap_nt_assuming_third1inclusive': overlap,
                         'locus_compatible': compatible,
                         'exact_old_AA_matches': '|'.join(old.loc[old.literal_standard_AA.map(aa).eq(seq), 'source_name'].astype(str))})
    matched = pd.DataFrame(rows)
    matched.to_csv(OUT / 'all_exact_AA_match_metadata_before_outcomes.csv', index=False)
    good = matched[matched.locus_compatible].copy()
    good['newer_match_count'] = good.groupby('third_new_id').source_ORF_ID.transform('nunique')
    good['third_match_count'] = good.groupby('source_ORF_ID').third_new_id.transform('nunique')
    good['metadata_one_to_one'] = good.newer_match_count.eq(1) & good.third_match_count.eq(1)
    good.to_csv(OUT / 'all_locus_compatible_exact_AA_metadata_before_outcomes.csv', index=False)
    report = {'utc': dt.datetime.now(dt.timezone.utc).isoformat(),
              'metadata_only': True, 'third_rows': len(third), 'newer_rows': len(new),
              'exact_AA_pairs': len(matched), 'locus_compatible_pairs': len(good),
              'one_to_one_pairs': int(good.metadata_one_to_one.sum()),
              'one_to_one_no_source_CCDS_overlap_pairs': int((good.metadata_one_to_one & good.third_any_ccds_overlap.eq('FALSE')).sum()),
              'third_guide_ids_are_sequences': False,
              'reagent_independence_established': False,
              'guide_level_files': 'Publisher direct resources returned403; stopped. No mirrors/retries.',
              'examples': good[['source_ORF_ID','third_new_id','third_pool_id','third_loci','third_any_ccds_overlap','metadata_one_to_one']].head(8).to_dict('records'),
              'metadata_sha256': hashlib.sha256((OUT / 'all_locus_compatible_exact_AA_metadata_before_outcomes.csv').read_bytes()).hexdigest()}
    (OUT / 'IDENTITY_PREFLIGHT.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__ == '__main__':
    main()
