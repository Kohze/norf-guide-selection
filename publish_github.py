"""Publish this reviewed repository after local GitHub CLI authentication."""
from pathlib import Path
import argparse,json,os,re,shutil,subprocess,sys
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--repo',default='Kohze/norf-guide-selection')
p.add_argument('--check-only',action='store_true',help='Check local Git without uploading or accessing GitHub.')
args=p.parse_args()
assert re.fullmatch(r'[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+',args.repo)
installed=Path(os.environ.get('ProgramFiles','C:/Program Files'))/'Git/cmd/git.exe'
git=str(installed) if os.name=='nt' and installed.is_file() else shutil.which('git')
if git is None:sys.exit('Git was not found. Install Git for Windows and reopen PowerShell.')
env=os.environ.copy()
# The parent shell must not select a different repository or index.
for key in ['GIT_DIR','GIT_WORK_TREE','GIT_INDEX_FILE','GIT_COMMON_DIR',
            'GIT_OBJECT_DIRECTORY','GIT_ALTERNATE_OBJECT_DIRECTORIES','GIT_NAMESPACE','GIT_PREFIX']:
    env.pop(key,None)
# Trust this explicitly requested release only, for this process only.
count=int(env.get('GIT_CONFIG_COUNT','0'))
env['GIT_CONFIG_KEY_'+str(count)]='safe.directory'
env['GIT_CONFIG_VALUE_'+str(count)]=ROOT.as_posix()
env['GIT_CONFIG_COUNT']=str(count+1)
env['PATH']=str(Path(git).parent)+os.pathsep+env.get('PATH','')

def run(cmd,check=True):
    command=[git,'-C',str(ROOT),*cmd[1:]] if cmd[0]=='git' else cmd
    try:
        result=subprocess.run(command,cwd=ROOT,env=env,text=True,
            encoding='utf-8',errors='replace',capture_output=True)
    except OSError as exc:
        sys.exit('Cannot start '+cmd[0]+': '+str(exc))
    if check and result.returncode:
        detail=result.stderr.strip() or result.stdout.strip() or '(no error message returned)'
        sys.exit('Command failed (exit '+str(result.returncode)+'): '+subprocess.list2cmdline(command)+'\n'+detail)
    return result

print('Git executable: '+git,flush=True)
print('Release folder: '+str(ROOT),flush=True)
top=Path(run(['git','rev-parse','--show-toplevel']).stdout.strip()).resolve()
if top!=ROOT:sys.exit('Git selected a different repository; no publication attempted.')
if args.check_only:
    status=run(['git','status','--porcelain']).stdout.strip()
    if status:sys.exit('Review and commit local changes before publishing:\n'+status)
    print(json.dumps({'status':'local_checks_passed','commit':run(['git','rev-parse','HEAD']).stdout.strip(),
                     'network_accessed':False,'uploaded':False},indent=2))
    sys.exit(0)
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
