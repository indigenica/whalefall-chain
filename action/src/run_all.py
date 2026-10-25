"""Run this isolated folder only. Fail on numerical contract violations."""
from pathlib import Path
import subprocess,sys,json,time,hashlib,datetime,os
ROOT=Path(__file__).resolve().parents[1]
# U26: the plateau (see common.py); pole -> results_pole/, logs_pole/
PLATEAU=os.environ.get('WHALEFALL_PLATEAU','pole')
RES=ROOT/('results' if PLATEAU=='sum' else 'results_pole'); LOGS=ROOT/('logs' if PLATEAU=='sum' else 'logs_pole')
for d in (RES,LOGS,ROOT/'runtime'): d.mkdir(exist_ok=True)
modules=['symbolic','portal','experiments','decays','late_fluid','plots']
record=dict(plateau=PLATEAU,start_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),python=sys.version,executable=sys.executable,modules=[])
for name in modules:
    start=time.monotonic();src=ROOT/'src'/f'{name}.py';log=LOGS/f'{name}.log'
    print('START',name,flush=True)
    with log.open('w',encoding='utf-8') as f:
        p=subprocess.run([sys.executable,'-B',str(src)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,
            env=dict(os.environ,PYTHONIOENCODING='utf-8',PYTHONDONTWRITEBYTECODE='1'))
    record['modules'].append(dict(module=name,returncode=p.returncode,seconds=time.monotonic()-start,
        code_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),log_sha256=hashlib.sha256(log.read_bytes()).hexdigest()))
    (RES/'execution.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    print('END',name,p.returncode,round(time.monotonic()-start,2),flush=True)
    if p.returncode:
        print(log.read_text(encoding='utf-8')[-2500:],flush=True);sys.exit(p.returncode)
checks={
    'symbolic_exact':json.loads((RES/'symbolic.json').read_text())['pass_exact'],
    'early_numerics':all(c['pass_numerics'] for c in json.loads((RES/'experiments.json').read_text())['cases']),
    'late_numerics':all(c['pass_numerics'] for c in json.loads((RES/'late_fluid.json').read_text())['cases']),
    'heavy_decay_channels':json.loads((RES/'decays.json').read_text())['all_modes_have_channels']}
record['checks']=checks;record['pass_execution']=all(checks.values())
record['end_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
(RES/'execution.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
if not record['pass_execution']:raise RuntimeError(checks)
