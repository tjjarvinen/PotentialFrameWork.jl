"""Render complete tables from raw records; no hard-coded performance results."""
import json, re, sys, tomllib
from pathlib import Path
out=Path(sys.argv[1]); meta=tomllib.loads((out/'metadata.toml').read_text()); records=tomllib.loads((out/'records.toml').read_text())['records']
keys=['polynomial_public','et_polynomial_public','spline_public','optimized_public','dag_public','optimized_prepared','dag_prepared']
labels=['ACE polynomial','ETACE polynomial','ACE spline','PFW prefix public','PFW DAG public','PFW prefix prepared','PFW DAG prepared']
summary={'fixtures':len(records),'trials':0,'samples':0,'undersampled':[],'maximum_errors':{},'tests':{}}
lines=['# Complete May30 observations','','All energies and virials are in eV; forces are in eV/Å. Energy error is per atom.','']
for r in records:
    lines += [f"## {r['fixture']}",'',f"Composition in Ti/Al/V order: {r['composition']}; directed neighbors: {r['directed_neighbors']}. Face heights: {r['face_heights_A']} Å; minimum distance: {r['minimum_distance_A']:.12g} Å. Verified center-1 V–V nearest image: atom {r['VV_neighbor']}, shift {r['VV_image']}, distance {r['VV_distance_A']:.12g} Å.",'',f"Cell (column-major): `{r['cell']}`. Geometry SHA256: `{r['geometry_sha256']}`.",'','| Comparison | Energy/atom | Force component | Virial component |','|---|---:|---:|---:|']
    for k,d in r['checks'].items():
        vals=[d[v] for v in ['energy_per_atom_eV','force_component_eV_per_A','virial_component_eV']]
        assert d['passed']; lines.append(f"| {k} | {' | '.join(f'{v:.5g}' for v in vals)} |")
        prev=summary['maximum_errors'].setdefault(k,[0.,0.,0.]); summary['maximum_errors'][k]=[max(x,y) for x,y in zip(prev,vals)]
    lines+=['','| Boundary | Median ms | Min–max ms | Bytes/call | Allocations | Samples |','|---|---:|---:|---:|---:|---:|']
    for k,label in zip(keys,labels):
        t=r['timings'][k]; summary['trials']+=1; summary['samples']+=t['samples']
        if t['undersampled']: summary['undersampled'].append([r['fixture'],k,t['samples']])
        lines.append(f"| {label} | {t['median_ns']/1e6:.4f} | {t['min_ns']/1e6:.4f}–{t['max_ns']/1e6:.4f} | {t['bytes']:,} | {t['allocations']:,} | {t['samples']} |")
    lines+=['',f"Setup elapsed seconds (excluded from warmed trials): `{r['setup']}`.",'']
initial=out/'initial-records.toml'
if initial.exists():
    old=tomllib.loads(initial.read_text())['records']
    missing=[(r['fixture'],k,t['samples']) for r in old for k,t in r['timings'].items() if t['undersampled']]
    summary['initial_undersampled']=missing
    summary['executed_trials']=summary['trials']+len(missing)
    summary['all_measured_samples']=summary['samples']+sum(n for _,_,n in missing)
for label,file in [('focused113','focused-113.log'),('focused110','focused-110.log'),('baseline113','baseline-tests-113.log'),('current113','tests-113.log'),('current110','tests-110.log')]:
    p=out/file
    if p.exists():
        text=p.read_text(); counts=re.findall(r'\|\s+(\d+)\s+(\d+)\s+[^\n]+',text)
        summary['tests'][label]={'passed':sum(int(a) for a,b in counts),'testsets':len(counts),'failure_markers':bool(re.search(r'Test Failed|Error During Test|ERROR:',text))}
(out/'OBSERVATIONS.md').write_text('\n'.join(lines))
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
