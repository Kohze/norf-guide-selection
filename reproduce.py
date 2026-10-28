"""Portable, isolated reproduction of released nORF screen analyses."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parent
def io_path(p):
    p=Path(p).resolve()
    return Path('\\\\?\\'+str(p)) if os.name=='nt' and not str(p).startswith('\\\\?\\') else p
p=argparse.ArgumentParser();p.add_argument('--output',default='results/reproduction');p.add_argument('--full-resampling',action='store_true');p.add_argument('--skip-figures',action='store_true');args=p.parse_args()
OUT=io_path(Path(args.output));OUT.mkdir(parents=True,exist_ok=False)
manifest=json.loads((ROOT/'RELEASE_MANIFEST.json').read_text())
for row in manifest['files']:
    path=ROOT/row['path']
    assert hashlib.sha256(path.read_bytes()).hexdigest()==row['sha256'],row['path']
for folder in ['supplementary','figure_source_data','evidence','code','provenance']:
    shutil.copytree(ROOT/folder,OUT/folder)
for name in ['REPORTED_RESULTS.json','SUPPLEMENTARY_TABLE_MANIFEST.json']:shutil.copy2(ROOT/name,OUT/name)
(OUT/'figures').mkdir()
logs=[]
def run(relative,*arguments):
    log=OUT/(str(len(logs)+1).zfill(2)+'_'+Path(relative).stem+'.log')
    print('Running',relative,flush=True)
    with log.open('w',encoding='utf-8') as f:
        r=subprocess.run([sys.executable,str(OUT/relative),*arguments],cwd=OUT,stdout=f,stderr=subprocess.STDOUT)
    logs.append({'script':relative,'exit_code':r.returncode,'log':log.name})
    if r.returncode:
        print(log.read_text(encoding='utf-8',errors='replace')[-5000:]);raise RuntimeError(relative)
run('evidence/analysis_extensions/cycle_90_within_product_shared_guide_concordance/analyze_frozen_within_product_shared_guide_concordance.py')
run('evidence/analysis_extensions/cycle_93_transferred_guide_selection/analyze_guide_selection.py')
run('code/frozen_replay.py')
run('code/analyze_context_exclusion.py');run('code/verify_context_exclusion.py')
run('code/analyze_efficacy_control.py');run('code/verify_efficacy_control.py')
if args.full_resampling:
    run('evidence/analysis_extensions/cycle_88_balanced_guide_count_cluster_uncertainty/run_balanced_guide_count_cluster_uncertainty.py')
    fresh=json.loads((OUT/'evidence/analysis_extensions/cycle_88_balanced_guide_count_cluster_uncertainty/BALANCED_GUIDE_COUNT_CLUSTER_UNCERTAINTY_RESULTS.json').read_text())
    expected=json.loads((ROOT/'REPORTED_RESULTS.json').read_text());import numpy as np
    q=fresh['estimates']['reused_minus_fresh_rank_agreement']
    assert np.allclose(q['point'],expected['count_balanced_product_agreement_gap'],rtol=0,atol=2e-12)
    assert np.allclose(q['CI97_5'],expected['balanced_gap_CI97_5'],rtol=0,atol=2e-12)
if not args.skip_figures:run('code/make_journal_figures.py')
receipt={'utc':datetime.now(timezone.utc).isoformat(),'status':'pass','manifest_files_checked':len(manifest['files']),
         'source_data_unchanged':True,'output':str(OUT),'full_balanced_resampling_regenerated':args.full_resampling,
         'figure_masters_regenerated':not args.skip_figures,'steps':logs,
         'scope':'Frozen-input statistical replay and figure rebuilding; original-source fetch and genome geometry scans are separate.'}
(OUT/'REPRODUCTION_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
