from pathlib import Path
import datetime as dt
import hashlib
import json

OUT = Path(__file__).resolve().parent

if __name__ == '__main__':
    path = OUT / 'FROZEN_PRODUCT_CONCORDANCE_PROTOCOL.json'
    assert not path.exists(), 'Protocol immutable; do not overwrite'
    cohort = OUT / 'all_metadata_eligibility_and_dependence_before_outcomes.csv'
    protocol = {
        'utc': dt.datetime.now(dt.timezone.utc).isoformat(),
        'design': 'Retrospective extension frozen before first third-study numeric outcome comparison; not preregistration.',
        'known_results': 'C89-C93 and metadata-only C94 are known. Third-study fitness outcome columns have not been loaded.',
        'changed_scope_before_outcomes': 'The planned reagent-independent test cannot be done because S3/S4 downloads are blocked. Proceed only with product-level cross-study concordance; never call it independent-guide validation.',
        'identity': 'Exact complete AA sequence with only whitespace/terminal-stop formatting normalized, one-to-one source identities, same chromosome/strand and exact coding interval chain allowing only terminal-stop inclusion. All239 compatible pairs exactly match under third-study0basedhalfopen intervals; none under1basedinclusive.',
        'eligibility': 'Third-study any_ccds_overlap FALSE and >=4 designed guides; source2026 >=4 unique20nt guides with primary_geometry_eligible, unique primary NGG site, a single source target assignment, and all3pDNA>=10. No outcome filtering.',
        'metadata_products': 196,
        'metadata_clusters': 192,
        'cohort_sha256': hashlib.sha256(cohort.read_bytes()).hexdigest(),
        'third_source': 'Pinned author-repository gene_summaries_final_v4.csv, verified against git blob. Reversible latin1 decoding. Untreated A375 only.',
        'third_outcomes': ['A375_r1_1_0_lfc','A375_r2_1_0_lfc'],
        'third_aggregation': 'Arithmetic mean of the two reported product MAGeCK LFCs; these are source perproduct guide summaries, not individual-guide measurements.',
        'source2026_outcome': 'Arithmetic mean of qualified guide LFCs in A375; sourceC83 meanlog2(late processed abundance+1)-meanlog2(pDNA+1).',
        'primary_endpoint': 'Spearman correlation across all eligible products with finite measurements in both sources. Equal product weight.',
        'primary_uncertainty': '2000 nonparametric dependence-cluster bootstrap draws, seed20261021. All products in a sampled cluster included and repeated explicitly; ranks recomputed for every draw. Nominal percentile95CI; no confirmatory Pvalue.',
        'dependence': 'Union same stable gene, overlapping coding bodies regardless strand, or shared third-study guide IDs.192 metadata clusters.',
        'minimum_products': 20,
        'missing': 'No imputation; list every missing product and preserve all metadata-eligible rows.',
        'secondary_descriptive': [
            'Pearson correlation.',
            'Each third-study replicate compared to the same source2026 product mean.',
            'Within-third-study replicate concordance.',
            'Source2026 median-of-guides instead of mean.',
            'Source2026 splice-safe>50nt from CDS-bearing RNA boundaries, >=4 eligible guides.',
            'One lexicographically first product per metadata dependence cluster.',
            'Subset also in original127 exact2021-2026 product inventory.',
            'All single-cluster omissions.',
            'Top20percent depletion-ranked product overlap: k=ceil(.2*n), ties by source2026ID. Random-selection expectation k*k/n. No target classification or causal conclusion.'
        ],
        'no_new_model_inference': True,
        'claim_limits': [
            'No guide-sequence independence established.',
            'Source2025 D4-toD21 and source2026 plasmid-toendpoint differ in baselines, guide filtering and summarization.',
            'Study/cell conditions and source25CCDS22 vs currentsource26CDS annotations differ.',
            'Only A375 evaluated; no normal tissues or therapeutic window.',
            'Functional locus concordance does not prove that the protein product causes the phenotype.',
            'All intervals are exploratory nominal intervals in a long retrospective research program.'
        ],
        'promotion_rule': 'Report the estimate and uncertainty regardless direction. Do not promote medical leads or independent-guide confirmation from this analysis. Decide manuscript integration only after verification and assessment of material scientific value.'
    }
    path.write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'frozen':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'products':196,'clusters':192}))
