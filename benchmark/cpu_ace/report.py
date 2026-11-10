#!/usr/bin/env python3
"""Summarize raw stage-one TOML; timings remain labelled by their boundary."""
import argparse, pathlib, statistics, tomllib
p=argparse.ArgumentParser()
p.add_argument('directory',type=pathlib.Path)
p.add_argument('--output',type=pathlib.Path)
a=p.parse_args()
reports=[]
for file in sorted(a.directory.glob('*.toml')):
    d=tomllib.loads(file.read_text())
    if d.get('schema')=='PFW.CPUACE.Stage1.v1': reports.append((file,d))
key=lambda c:(c['model'],c['atoms'],c['precision'],c['width'])
baselines={}
for _,d in reports:
    for c in d['cases']:
        if c['algorithm']=='historical': baselines[(d['threads'],key(c))]=c
lines=['# CPU ACE stage-one measured results','',
       'Warmed BenchmarkTools medians. Ratios compare identical numerical models and boundaries; raw samples and allocation observations remain in the TOML files.', '',
       '| Threads | Model | Atoms | Precision / width | Reference ms | Optimized ms | Incremental ratio | Historical ratio | Public ms | External public ms |',
       '|---:|---|---:|---|---:|---:|---:|---:|---:|---:|']
for file,d in reports:
    refs={key(c):c for c in d['cases'] if c['algorithm']=='reference'}
    for c in d['cases']:
        if c['algorithm']!='optimized': continue
        ref=refs[key(c)]; base=baselines.get((d['threads'],key(c)))
        t=c['prepared']['median_s']
        h=f"{base['prepared']['median_s']/t:.2f}x" if base else 'unmeasured'
        ext=statistics.median(c['external_raw_s'])*1e3 if 'external_raw_s' in c else None
        lines.append(f"| {d['threads']} | {c['model']} | {c['atoms']} | {c['precision']} / {c['width']} | {ref['prepared']['median_s']*1e3:.3f} | {t*1e3:.3f} | {ref['prepared']['median_s']/t:.2f}x | {h} | {c['complete_public']['median_s']*1e3:.3f} | {ext:.3f} |" if ext is not None else
                     f"| {d['threads']} | {c['model']} | {c['atoms']} | {c['precision']} / {c['width']} | {ref['prepared']['median_s']*1e3:.3f} | {t*1e3:.3f} | {ref['prepared']['median_s']/t:.2f}x | {h} | {c['complete_public']['median_s']*1e3:.3f} | unavailable |")
lines+=['','External times are the stored oracle thread count/boundary, not the PFW thread count. Modern references are combined Float64 public E/F/V; Float32 rows show a Float64 sensitivity comparison. Legacy references are explicitly three separate analytic ACE1 calls with JuLIP native neighbor caching at fixed geometry. Do not interpret the prepared/external ratio as a pure kernel speedup.','',
        '| Model / precision | Algorithm | First-site product µs | Density µs | Radial µs | Angular µs | Force back-contraction µs |',
        '|---|---|---:|---:|---:|---:|---:|']
for _,d in reports:
    for c in d['cases']:
        s=c.get('isolated_site_stages')
        if not s: continue
        vals=[s[k]['median_s']*1e6 for k in ['products_reverse','density','radial','angular','force_back_contraction']]
        lines.append(f"| {c['model']} / {c['precision']} / w{c['width']} | {c['algorithm']} | "+' | '.join(f'{v:.3f}' for v in vals)+' |')
lines+=['','Isolated site stages perturb state/cache behavior and cannot be added to reconstruct total execution time. Preparation observations retain species packing and colored scheduling separately.','', 'Under-sampled trials:','']
under=[]
for file,d in reports:
    for c in d['cases']:
        for boundary in ('prepared','complete_public','repeated_frames'):
            s=c[boundary]
            if s['undersampled']: under.append(f"- `{file.name}` {c['model']} / {c['atoms']} / {c['precision']} / w{c['width']} / {c['algorithm']} / {boundary}: {s['samples']} samples.")
lines+=under or ['None for the three reported whole-call boundaries.']
lines+=['','Omissions:','']
for file,d in reports:
    for o in d.get('omissions',[]): lines.append(f"- `{file.name}` {o['model']}: {o['status']}: {o['reason']}")
lines+=['','Source identities, package versions, scratch payloads, geometry/model hashes, absolute E/F/V errors and raw observations are in each TOML.']
text='\n'.join(lines)+'\n'
if a.output:a.output.write_text(text)
else:print(text)
