#!/usr/bin/env python3
"""Prepare a SEPARATE flat MATLAB/Octave workspace without rewriting original scripts."""
import argparse
import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

def copy(source, target):
    src = REPO / source
    if not src.is_file():
        raise FileNotFoundError(src)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, target)

def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--output', required=True, type=Path, help='A NEW (or empty) directory, outside tracked source folders')
    cli.add_argument('--mode', choices=['check', 'optimize'], default='check')
    args=cli.parse_args()
    output=args.output.expanduser().resolve()
    if output.exists() and any(output.iterdir()):
        raise SystemExit('Refusing to overwrite a non-empty directory: ' + str(output))
    if output == REPO or REPO in output.parents and output.name in {'src','data','results','tools'}:
        raise SystemExit('Refusing to write over a tracked repository directory')
    output.mkdir(parents=True,exist_ok=True)
    bundle = output / 'ALEXIS_v4_short_names_bundle'
    for name in ('myloads.mat','mypvs.mat','myloadspvs.mat'):
        copy(Path('data/profiles')/name, bundle/name)
    for name in ('control_schedule_33bus.csv','control_schedule_69bus.csv'):
        copy(Path('data/historical_schedules')/name, bundle/name)
    copy('data/networks/case69modified.m',bundle/'case69modified.m')
    for f in (REPO/'src/v4_replay').glob('*.m'):
        copy(f.relative_to(REPO),bundle/f.name)
    for f in (REPO/'src/final_stage2_69bus').glob('*.m'):
        copy(f.relative_to(REPO),output/f.name)
    if args.mode=='check':
        for f in (REPO/'results/final_stage2_69bus').iterdir():
            if f.is_file():
                copy(f.relative_to(REPO),output/'results_stage2_reopt_69_1MW'/f.name)
    print('Workspace prepared:', output)
    print('Mode:',args.mode)
    print('Next: export MATPOWER_DIR=/path/to/matpower')
    print('Then: cd ' + str(output))
    cmd = "stage2_reopt_69_octave('check')" if args.mode == 'check' else "stage2_reopt_69_octave('optimize')"
    print('octave --no-gui --quiet --eval "' + cmd + '"')
    if args.mode=='check':
        print('After the check, optional read-only: octave --no-gui --quiet --eval "audit_stage2_69_losses_octave()"')
        print('After the check, optional read-only: octave --no-gui --quiet --eval "audit_stage2_69_final_comparisons()"')

if __name__=='__main__': main()
