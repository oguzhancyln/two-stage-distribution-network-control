#!/usr/bin/env python3
"""Read-only schedule / radial topology checker. No numerical power flows or optimization.

Usage:
  python independent_radiality_audit.py --case 33 --branch-file cases/branch_endpoints_case33bw.csv --schedule schedules/control_schedule_33bus.csv
  python independent_radiality_audit.py --case 69 --branch-file cases/case69modified.m --schedule schedules/control_schedule_69bus_STAGE2_REOPT.csv
Input branches are 1-based and indexed by their row order in the MATPOWER case file.
"""
import argparse, csv, pathlib, re, sys

def parse_branches(path):
    s=path.read_text(encoding='utf-8')
    if path.suffix.lower() == '.csv':
        with path.open(newline='') as fp:
            rows=list(csv.DictReader(fp))
        if not rows or set(rows[0]) != {'index','from_bus','to_bus'}:
            raise ValueError('Expected index,from_bus,to_bus CSV')
        triples=[(int(r['index']),int(r['from_bus']),int(r['to_bus'])) for r in rows]
    else:
        m=re.search(r'\bmpc\.branch\s*=\s*\[(.*?)\];',s,re.S)
        if not m: raise ValueError('Cannot locate mpc.branch matrix in '+str(path))
        triples=[]
        for ln in m.group(1).splitlines():
            ln=ln.split('%')[0].strip().rstrip(';').strip()
            if not ln: continue
            vals=ln.split()
            if len(vals)<2: raise ValueError('Bad branch row '+repr(ln))
            triples.append((len(triples)+1,int(float(vals[0])),int(float(vals[1]))))
    if [x[0] for x in triples] != list(range(1,len(triples)+1)):
        raise ValueError('Branch indices must be consecutive, one-based and ordered')
    return triples

def parse_vect(txt):
    txt=txt.strip()
    if not (txt.startswith('[') and txt.endswith(']')):
        raise ValueError('Invalid vector: '+txt)
    st=txt[1:-1].strip()
    return [int(v) for v in st.replace(',',' ').split()] if st else []

def radiality(buses, branches, opens):
    active=[(a,b) for i,a,b in branches if i not in opens]
    if len(active)!=buses-1: return False
    parent=list(range(buses+1))
    def root(z):
        while parent[z]!=z:
            parent[z]=parent[parent[z]]
            z=parent[z]
        return z
    for a,b in active:
        if a not in range(1,buses+1) or b not in range(1,buses+1): raise ValueError('out-of-range bus in branch input')
        ra,rb=root(a),root(b)
        if ra==rb:return False
        parent[ra]=rb
    return len({root(i) for i in range(1,buses+1)})==1

def main():
    pa=argparse.ArgumentParser()
    pa.add_argument('--case',type=int,choices=[33,69],required=True)
    pa.add_argument('--branch-file',type=pathlib.Path,required=True)
    pa.add_argument('--schedule',type=pathlib.Path,required=True)
    a=pa.parse_args()
    branch=parse_branches(a.branch_file)
    with a.schedule.open(newline='') as fp: rows=list(csv.DictReader(fp))
    if len(rows)!=96: raise ValueError('Expected exactly 96 schedule rows, got '+str(len(rows)))
    if a.case==33:
        assert len(branch)==37,(len(branch),37)
        initial_open=set(range(33,38))
        max_caps=3
    else:
        assert len(branch)==73,(len(branch),73)
        initial_open=set(range(69,74))
        max_caps=5
    n_open=len(branch)-a.case+1
    oldopens=initial_open
    oldcaps=[0]*max_caps
    switches_init=switches_next=cap_init=cap_next=0
    prior=None
    abnormal=[]
    for j,r in enumerate(rows,1):
        if int(r['TimeStep']) != j: raise ValueError('TimeStep row mismatch '+str(j))
        op=set(parse_vect(r['OpenLines']))
        caps=parse_vect(r['CapSteps'])
        if len(op)!=n_open or min(op)<1 or max(op)>len(branch):
            raise ValueError('Wrong open-line count/range at '+str(j))
        if len(caps)!=max_caps or any(c<0 or c>5 for c in caps):
            raise ValueError('Wrong capacitor count/range at '+str(j))
        if not radiality(a.case,branch,op):
            raise ValueError('Non-radial network at '+str(j))
        if j==1:
            switches_init=len(op.symmetric_difference(initial_open))
            cap_init=sum(caps)
        else:
            switches_next+=len(op.symmetric_difference(oldopens))
            cap_next+=sum(abs(p-q) for p,q in zip(caps,oldcaps))
            if (op!=oldopens or caps!=oldcaps) and j not in [40,73]:
                abnormal.append(j)
        oldopens,oldcaps=op,caps
    data={'case':a.case,'radial_steps':len(rows),'line_initial':switches_init,'line_later':switches_next,'cap_initial':cap_init,'cap_later':cap_next,'unscheduled_changes':abnormal}
    for k,v in data.items(): print(f'{k}: {v}')
    if abnormal: raise ValueError('Unexpected discrete switching outside Stage-1 boundary')
    expected = (8,8,13,8) if a.case==33 else (10,16,14,13)
    actual=(switches_init,switches_next,cap_init,cap_next)
    if actual != expected: raise ValueError(f'Count does not match reported diagnostic: {actual} != {expected}')
    print('PASS: 96/96 radial endpoint checks and switching counts match manuscript')
if __name__=='__main__':main()
