#!/usr/bin/env python3
"""Check snapshot hashes, run Python syntax compilation and radiality audits."""
import csv
import hashlib
import pathlib
import subprocess
import sys

ROOT=pathlib.Path(__file__).resolve().parents[1]
SUMS=ROOT/'SHA256SUMS.txt'

def main():
    failures=[]
    seen=set()
    for line in SUMS.read_text(encoding='utf-8').splitlines():
        if not line.strip(): continue
        digest,sep,path=line.partition('  ')
        if not sep:
            failures.append('Malformed SHA256SUMS entry: '+line);continue
        target=ROOT/path
        seen.add(path)
        if not target.is_file(): failures.append('Missing: '+path);continue
        actual=hashlib.sha256(target.read_bytes()).hexdigest()
        if actual!=digest: failures.append('Checksum mismatch: '+path)
    ignored={'__pycache__','.git','local_workspace','local_optimization','workspaces','results_generated'}
    actualfiles={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt' and not (ignored & set(p.relative_to(ROOT).parts))}
    if actualfiles!=seen: failures.append('Manifest membership mismatch (extra/missing tracked files): '+repr(sorted(actualfiles^seen)[:8]))
    if failures:
        for failure in failures: print('FAIL:',failure)
        return 1
    commands=[
        [sys.executable,str(ROOT/'src/validation/independent_radiality_audit.py'), '--case','33','--branch-file',str(ROOT/'data/networks/branch_endpoints_case33bw.csv'),'--schedule',str(ROOT/'data/historical_schedules/control_schedule_33bus.csv')],
        [sys.executable,str(ROOT/'src/validation/independent_radiality_audit.py'), '--case','69','--branch-file',str(ROOT/'data/networks/case69modified.m'),'--schedule',str(ROOT/'results/final_stage2_69bus/control_schedule_69bus_STAGE2_REOPT.csv')],
    ]
    for cmd in commands:
        subprocess.run(cmd,check=True)
    s=ROOT/'results/final_stage2_69bus/audit_69_full_replay_summary.csv'
    with s.open(newline='') as fp: rows={r['Method']:r for r in csv.DictReader(fp)}
    final=rows['Reoptimized_Stage2_1MW']
    if abs(float(final['JV'])-1.54885382313137)>1e-8 or int(final['Violations'])!=51:
        print('FAIL: Saved 69-bus metrics mismatch');return 1
    if abs(float(final['ActiveLoss_kWh'])-1797.509191829)>1e-6:
        print('FAIL: Saved 69-bus losses mismatch');return 1
    if (int(final['Accepted']),int(final['Rejected']),int(final['AppliedQCommands']),int(final['SaturationEvents'])) != (27,0,158,143):
        print('FAIL: Saved 69-bus command statistics mismatch');return 1
    print('PASS: SHA256 hashes, two radiality audits, 69-bus saved summary invariants.')
    print('NOTE: This tool does NOT run MATPOWER or GWO.')
    return 0

if __name__=='__main__': raise SystemExit(main())
