"""Bounded public-data intake. Never loads fitness outcomes at this stage."""
from pathlib import Path
import datetime as dt
import hashlib
import json
import shutil
import urllib.request

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
OLD = ROOT / 'exploration_2026-10-11/novelty_discovery'
PRIOR = OLD / 'cycle_91_independent_sORF_screen_guide_validation'
META = ['new_id', 'id', 'loci_all', 'num_guides', 'any_ccds_overlap',
        'ccds_overlap_per_exon', 'Amino.acid.sequence',
        'Amino.acid.sequence_iniMet', 'Microprotein.length', 'old_id',
        'guides_list', 'peptide_id_only', 'source_id', 'uniprot_id', 'gene_id']

def save(name, obj):
    (OUT / name).write_text(json.dumps(obj, indent=2) + '\n', encoding='utf-8')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    baseline = OUT / 'baseline'
    baseline.mkdir(exist_ok=True)
    for name in ['CURRENT_NORFS_MANUSCRIPT.json', 'SESSION_RESUME_STATE_2026-10-16.json']:
        dest = baseline / name
        if not dest.exists():
            shutil.copyfile(ROOT / name, dest)
    scope = OUT / 'INTAKE_SCOPE.json'
    if not scope.exists():
        save(scope.name, {
            'utc': dt.datetime.now(dt.timezone.utc).isoformat(),
            'baseline_manuscript': 14, 'cycle': 94,
            'question': 'Do exact nORF protein products show independent A375 fitness concordance when the third study uses different guide sequences?',
            'new_increment': 'Third-study product-level validation; not individual-guide validation.',
            'known_before_intake': 'C90-C93 outcomes and the C91 header/repository feasibility inspection are known. No full third-study table or cross-source outcomes were inspected.',
            'stage': 'metadata_only; freeze final eligibility and endpoint before opening outcome columns',
            'allowed_intake': 'Pinned author-repository processed gene-summary CSV, <=25 MB; no raw sequencing.',
            'metadata_columns': META,
            'candidate_identity': 'Exact complete AA sequence plus compatible locus; one-to-one identities primary. Initiator-M-only matches separate sensitivity, never silently substituted.',
            'guides': 'Inspect guide metadata and quantify exact20nt overlap before outcome comparisons. Missing guide sequences cannot establish reagent independence.',
            'common_context': 'A375 untreated; no 6TG or cross-cell substitution.',
            'planned_primary': 'Spearman correlation of independent-study mean untreated A375 replicate LFC against source2026 mean qualified guide A375 LFC; eligibility and dependence handling frozen after metadata feasibility.',
            'minimum_inferential_cohort': 20,
            'stop_rule': 'If exact identity/reagent independence is unavailable or fewer than20 eligible independent loci, save bounded feasibility/descriptive result without calling it external confirmation.',
            'interpretation_limits': ['DNA perturbation does not isolate peptide causality.', 'A375-only concordance does not establish disease selectivity or a therapeutic window.', 'No individual medical target declaration from correlation.'],
            'access_boundaries': 'Do not retry blocked author article or prior blocked sources; this fetch is independently purposed processed data.'
        })
    branch = json.loads((PRIOR / 'AUTHOR_DEFAULT_BRANCH_COMMIT.body').read_text())
    commit = branch['commit']['sha']
    tree = json.loads((PRIOR / 'AUTHOR_PUBLIC_REPOSITORY_TREE.body').read_text())
    entry = next(x for x in tree['tree'] if x['path'] == 'Figure1_S2/gene_summaries_final_v4.csv')
    url = f'https://raw.githubusercontent.com/elsasserlab/schlesinger_microprotein_screen/{commit}/{entry["path"]}'
    src = OUT / 'sources'
    src.mkdir(exist_ok=True)
    dest = src / 'gene_summaries_final_v4.csv'
    if not dest.exists():
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=35) as response:
            if response.status != 200:
                raise RuntimeError(f'Stop on HTTP {response.status}')
            data = response.read(25_000_001)
            if len(data) > 25_000_000:
                raise RuntimeError('Intake size bound exceeded')
            dest.write_bytes(data)
            save('SOURCE_RESPONSE.json', {
                'utc': dt.datetime.now(dt.timezone.utc).isoformat(),
                'url': url, 'status': response.status,
                'headers': dict(response.headers), 'bytes': len(data),
                'sha256': sha(dest), 'pinned_commit': commit,
                'git_blob_sha1': hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest(),
                'expected_git_blob_sha1': entry['sha'],
                'expected_bytes': entry['size']
            })
    receipt = json.loads((OUT / 'SOURCE_RESPONSE.json').read_text())
    assert dest.stat().st_size == entry['size']
    assert receipt['git_blob_sha1'] == entry['sha']
    assert sha(dest) == receipt['sha256']
    # Reversible mapping of the source's non-UTF8 bytes; no replacement loss.
    meta = pd.read_csv(dest, usecols=META, dtype=str, keep_default_na=False, encoding='latin1')
    meta.to_csv(OUT / 'all_third_study_identity_and_reagent_metadata.csv', index=False)
    save('METADATA_INTAKE_RESULT.json', {
        'utc': dt.datetime.now(dt.timezone.utc).isoformat(), 'rows': len(meta),
        'columns_loaded': list(meta), 'outcomes_loaded': False,
        'source_encoding': 'latin1 reversible decoding; raw bytes preserved',
        'source_verified': True, 'source_sha256': receipt['sha256'],
        'column_nonempty_counts': {c: int((meta[c] != '').sum()) for c in meta},
        'first5_metadata_rows': meta.head(5).to_dict('records')
    })
    print(json.dumps({'rows': len(meta), 'verified_source': True, 'outcomes_loaded': False,
                      'metadata_examples': meta[['new_id','loci_all','guides_list','gene_id','source_id','any_ccds_overlap']].head(4).to_dict('records')}, indent=2))

if __name__ == '__main__':
    main()
