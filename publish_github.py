"""Publish this reviewed repository after local GitHub CLI authentication."""
from pathlib import Path
import argparse,json,re,subprocess,sys
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--repo',default='Kohze/norf-guide-selection');args=p.parse_args()
assert re.fullmatch(r'[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+',args.repo)
def run(cmd,check=True):return subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True,check=check)
auth=run(['gh','api','user','--jq','.login'],False)
if auth.returncode:
    sys.exit('GitHub authentication is required. Run gh auth login -h github.com locally, then rerun this command. Do not put tokens in source files.')
if run(['git','status','--porcelain']).stdout.strip():sys.exit('Commit or review local changes before publishing.')
repo=run(['gh','repo','view',args.repo,'--json','url,isPrivate'],False)
if repo.returncode:
    result=run(['gh','repo','create',args.repo,'--public','--source',str(ROOT),'--remote','origin','--push'])
else:
    remote=run(['git','remote','get-url','origin'],False)
    wanted='https://github.com/'+args.repo+'.git'
    if remote.returncode:run(['git','remote','add','origin',wanted])
    elif remote.stdout.strip()!=wanted:sys.exit('Existing origin differs from the requested destination; no changes made.')
    # A normal push rejects unrelated or newer remote history; never force.
    result=run(['git','push','-u','origin','main'])
local=run(['git','rev-parse','HEAD']).stdout.strip()
remote=run(['gh','api','repos/'+args.repo+'/commits/main','--jq','.sha']).stdout.strip()
assert local==remote,'Remote commit verification failed.'
info=json.loads(run(['gh','repo','view',args.repo,'--json','url,isPrivate']).stdout)
assert not info['isPrivate'],'Expected the authorized public data/code release.'
receipt={'status':'published_and_verified','url':info['url'],'commit':local,'public':True,'force_push':False}
(ROOT.parent/(ROOT.name+'_publication_receipt.json')).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
print(json.dumps(receipt,indent=2))
