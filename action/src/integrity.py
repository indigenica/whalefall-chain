"""Read-only SHA-256 audit of originals and sealed previous campaign."""
from pathlib import Path
import hashlib,json,datetime,os
ROOT=Path(__file__).resolve().parents[1]
PREVIOUS=ROOT.parent/'cadena_cierre_20261004'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        while b:=f.read(4*1024*1024):h.update(b)
    return h.hexdigest()
def canon(p):return str(p.resolve()).removeprefix('\\\\?\\').casefold()
old=json.loads((PREVIOUS/'results/campaign_manifest.json').read_text())
previous_bad=[]
for rel,record in old['files'].items():
    p=PREVIOUS/rel
    if not p.is_file() or sha(p)!=record['sha256']:previous_bad.append(rel)
print('Previous campaign checked:',len(old['files']),'changes:',len(previous_bad),flush=True)
original=json.loads((PREVIOUS/'results/originals_manifest.json').read_text())
original_bad=[];total_bytes=0
for i,(name,record) in enumerate(original.items()):
    p=Path(name);total_bytes+=record['bytes']
    if not p.is_file() or sha(p)!=record['sha256']:original_bad.append(name)
    if i%500==0:print('Originals checked:',i,'/',len(original),flush=True)
known={canon(Path(k)) for k in original}
allowed={canon(ROOT),canon(PREVIOUS)};new_outside=[]
for base in [ROOT.parent,ROOT.parent.parent/'claude']:
    for folder,dirs,files in os.walk(base):
        # Match the baseline snapshot's exclusions. Existing bytecode caches
        # were not inventoried in that snapshot and cannot be called new files.
        dirs[:]=[d for d in dirs if d not in ('.git','__pycache__') and canon(Path(folder)/d) not in allowed]
        for name in files:
            p=Path(folder)/name
            if canon(p) not in known:new_outside.append(str(p))
sources=[(PREVIOUS/'results/thermal.json',ROOT/'vendor/thermal.json'),
    (PREVIOUS/'vendor/c32b_scars_kibble.py',ROOT/'vendor/c32b_scars_kibble.py'),
    (ROOT.parent.parent/'claude/tools/c33_q_decays_sympy.py',ROOT/'vendor/c33_q_decays_sympy.py'),
    (PREVIOUS/'src/common.py',ROOT/'src/common.py')]
copied=[dict(source=str(s),copy=str(p.relative_to(ROOT)),source_sha256=sha(s),copy_sha256=sha(p),identical=sha(s)==sha(p)) for s,p in sources]
files={str(p.relative_to(ROOT)):dict(bytes=p.stat().st_size,sha256=sha(p)) for top in ['src','vendor','logs','results']
    for p in (ROOT/top).rglob('*') if p.is_file() and p.name!='integrity.json'}
for p in ROOT.glob('*.md'):files[p.name]=dict(bytes=p.stat().st_size,sha256=sha(p))
result=dict(checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    originals_checked=len(original),original_bytes=total_bytes,originals_changed_or_removed=original_bad,
    previous_campaign_files_checked=len(old['files']),previous_campaign_changed_or_removed=previous_bad,
    new_files_outside_authorized_folders=new_outside,copied_sources=copied,files=files,
    scan_exclusions=['.git','__pycache__','the two authorized campaign folders'],
    passed=not(original_bad or previous_bad or new_outside) and all(c['identical'] for c in copied))
(ROOT/'results/integrity.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print('INTEGRITY',result['passed'],'Originals',len(original),'Previous',len(old['files']),'New artifacts',len(files),flush=True)
if not result['passed']:raise RuntimeError('Integrity audit failed; see results/integrity.json')
