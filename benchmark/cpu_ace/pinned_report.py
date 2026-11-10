#!/usr/bin/env python3
"""Join controlled-affinity PFW and modern oracle observations without changing raw data."""
import json, pathlib, statistics, sys, tomllib
root=pathlib.Path(sys.argv[1])
rows=[]
lines=['# Controlled-affinity repeats','',
       'Float64, 1,024 atoms, ten warmed samples per boundary. Four threads use physical CPUs 0–3; sixteen use 0–15. BLAS uses one thread. External modern public calls use the same CPU masks and thread counts.','',
       '| Threads | Model | Historical prepared ms | Current reference ms | Optimized ms | Incremental ratio | Total historical ratio | PFW public ms | Modern public ms |',
       '|---:|---|---:|---:|---:|---:|---:|---:|---:|']
for threads in (4,16):
    current=tomllib.loads((root/'pinned'/f'current-t{threads}.toml').read_text())
    historical=tomllib.loads((root/'pinned'/f'historical-t{threads}.toml').read_text())
    for name in ('small-tial','medium-ti','large-tial'):
        select=lambda data,algorithm:next(c for c in data['cases'] if c['model']==name and c['algorithm']==algorithm)
        r=select(current,'reference');o=select(current,'optimized');h=select(historical,'historical')
        ext=tomllib.loads((root/f'oracles-pinned-t{threads}'/(name+'-oracle.toml')).read_text())
        raw=next(c for c in ext['records'] if c['g']==8)
        ext_ms=statistics.median(raw['raw_ns'])*1e-6
        hm=h['prepared']['median_s']*1e3;rm=r['prepared']['median_s']*1e3;om=o['prepared']['median_s']*1e3
        pm=o['complete_public']['median_s']*1e3
        lines.append(f'| {threads} | {name} | {hm:.3f} | {rm:.3f} | {om:.3f} | {rm/om:.2f}x | {hm/om:.2f}x | {pm:.3f} | {ext_ms:.3f} |')
        rows.append(dict(threads=threads,model=name,atoms=1024,precision='Float64',affinity=current['affinity'],
            historical_prepared_ms=hm,reference_prepared_ms=rm,optimized_prepared_ms=om,incremental_ratio=rm/om,
            historical_ratio=hm/om,pfw_public_ms=pm,modern_public_ms=ext_ms,
            oracle_threads=ext['threads'],samples_pfw=o['prepared']['samples'],samples_modern=raw['samples']))
lines+=['','These repeats confirm gains while showing sensitivity to placement and process/compiler state. Historical/source and current-reference ratios are distinct observations; all differences between their executors cannot be attributed to the product primitive. Raw samples/spreads and independent E/F/V errors remain in the TOML files.','',
        'Public comparisons include fresh neighbor preparation and output construction on both sides. Prepared PFW timings exclude those costs and are not compared as pure kernel speedups against external public calls.']
(root/'pinned'/'REPORT.md').write_text('\n'.join(lines)+'\n')
(root/'pinned'/'comparisons.json').write_text(json.dumps(dict(schema='PFW.CPUACE.PinnedComparison.v1',cases=rows),indent=2)+'\n')
