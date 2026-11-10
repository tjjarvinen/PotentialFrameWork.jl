"""Summarize preserved raw DAG observations; no benchmark execution."""
import json, pathlib, sys, tomllib
root=pathlib.Path(sys.argv[1])
def read(name): return tomllib.loads((root/name).read_text())
def ms(c,key): return c[key]['median_s']*1000
def group(cases):
    d={}
    for c in cases: d.setdefault((c['model'],c['atoms'],c['precision']),{})[c['algorithm']]=c
    return d
serial=read('serial.toml'); large=read('large.toml')
lines=['# DAG measured observations','', 'All timings below are warmed medians in milliseconds. Ratios are prefix / DAG. Raw TOML retains every sample, extrema, GC observations and allocation volumes.','',
       '| Model | Atoms | Precision | Polynomial prepared | Prefix prepared | DAG prepared | Gain | Prefix public | DAG public | Gain |',
       '|---|---:|---|---:|---:|---:|---:|---:|---:|---:|']
summary=[]
for key,g in group(serial['cases']+large['cases']).items():
    p,d=g['optimized'],g['dag']; r=g.get('reference')
    gain=ms(p,'prepared')/ms(d,'prepared'); public=ms(p,'complete_public')/ms(d,'complete_public')
    lines.append(f'| {key[0]} | {key[1]} | {key[2]} | {ms(r,"prepared"):.3f}' if r else f'| {key[0]} | {key[1]} | {key[2]} | —')
    lines[-1]+=f' | {ms(p,"prepared"):.3f} | {ms(d,"prepared"):.3f} | {gain:.2f}× | {ms(p,"complete_public"):.3f} | {ms(d,"complete_public"):.3f} | {public:.2f}× |'
    summary.append(dict(model=key[0],atoms=key[1],precision=key[2],prepared_gain=gain,public_gain=public,
                        fresh_gain=ms(p,'construct_and_public')/ms(d,'construct_and_public')))
lines+=['','## Fresh calculator plus complete public call','',
        '| Model | Atoms | Precision | Prefix ms | DAG ms | Gain |',
        '|---|---:|---|---:|---:|---:|']
for key,g in group(serial['cases']+large['cases']).items():
    p,d=g['optimized'],g['dag']
    lines.append(f'| {key[0]} | {key[1]} | {key[2]} | {ms(p,"construct_and_public"):.3f} | {ms(d,"construct_and_public"):.3f} | {ms(p,"construct_and_public")/ms(d,"construct_and_public"):.2f}× |')
lines+=['','## Product diagnostics at 128 atoms','',
        '| Model | Precision | Polynomial product/reverse | Prefix product/reverse | DAG forward | DAG reverse | DAG combined | Product gain | Readout |',
        '|---|---|---:|---:|---:|---:|---:|---:|---:|']
for key,g in group(serial['cases']).items():
    if key[1]!=128: continue
    r,p,d=g['reference']['stages'],g['optimized']['stages'],g['dag']['stages']
    lines.append(f'| {key[0]} | {key[2]} | {ms(r,"products_reverse"):.6f} | {ms(p,"products_reverse"):.6f} | {ms(d,"dag_forward"):.6f} | {ms(d,"dag_reverse"):.6f} | {ms(d,"products_reverse"):.6f} | {ms(p,"products_reverse")/ms(d,"products_reverse"):.2f}× | {ms(d,"readout"):.6f} |')
lines+=['','## Structural work (both precisions share indices)','',
        '| Model | Leaves | Original products | Orders | Unique terminals | Multiply nodes | Intermediate only | Depth | Sharing ratio | Prefix fwd/rev multiplies | DAG fwd/rev multiplies |',
        '|---|---:|---:|---|---:|---:|---:|---:|---:|---|---|']
for c in serial['cases']:
    if c['algorithm']!='dag' or c['atoms']!=128 or c['precision']!='Float64': continue
    g=c['graph']; orders=', '.join(f'{k}: {v}' for k,v in c['metadata']['product_orders'].items())
    lines.append(f'| {c["model"]} | {g["leaves"]} | {c["metadata"]["products"]} | {orders} | {g["unique_terminals"]} | {g["multiplication_nodes"]} | {g["intermediate_only_nodes"]} | {g["max_depth"]} | {g["sharing_ratio_vs_independent_nontrivial_products"]:.2f} | {g["prefix_forward_multiplications"]} / {g["prefix_reverse_multiplications"]} | {g["dag_forward_multiplications"]} / {g["dag_reverse_multiplications"]} |')
lines+=['','## Construction and memory at 128 atoms, one slot','',
        '| Model | Precision | Plan build ms | Plan build allocated bytes | Executor build ms | Plan array payload | DAG scratch | Retained polynomial | Total slot scratch | Fresh prefix public ms | Fresh DAG public ms |',
        '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
for c in serial['cases']:
    if c['algorithm']!='dag' or c['atoms']!=128: continue
    g=c['graph']; prefix=group(serial['cases'])[(c['model'],128,c['precision'])]['optimized']
    lines.append(f'| {c["model"]} | {c["precision"]} | {ms(c,"plan_build"):.4f} | {c["plan_build"]["host_bytes"]} | {ms(c,"executor_build"):.4f} | {g["plan_array_payload_bytes"]} | {g["dag_per_slot_bytes"]} | {g["retained_polynomial_per_slot_bytes"]} | {g["total_per_slot_bytes"]} | {ms(prefix,"construct_and_public"):.3f} | {ms(c,"construct_and_public"):.3f} |')
lines+=['','## Numerical disagreement, max across measured sizes','',
        '| Model | Precision | DAG vs polynomial E eV | F eV/Å | W eV | DAG vs independent Float64 E eV | F eV/Å | W eV |',
        '|---|---|---:|---:|---:|---:|---:|---:|']
for model in dict.fromkeys(c['model'] for c in serial['cases']):
    for precision in ('Float64','Float32'):
        cs=[c for c in serial['cases']+large['cases'] if c['algorithm']=='dag' and c['model']==model and c['precision']==precision]
        vals=[]
        for category in ('reference_errors','external_errors'):
            for key in ('energy_abs_eV','forces_max_eV_per_A','virial_max_eV'):
                vals.append(max(c[category][key] for c in cs if category in c))
        lines.append(f'| {model} | {precision} | '+' | '.join(f'{v:.5g}' for v in vals)+' |')
lines+=['','## Controlled physical-core scaling, 1,024 atoms, Float64','',
        '| Cores | Model | Prefix prepared ms | DAG prepared ms | Gain | Prefix public ms | DAG public ms | All-slot scratch bytes | Shared force bytes | DAG prepared host allocation bytes |',
        '|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|']
for threads in (1,4,16):
    for key,g in group(read(f'scaling-t{threads}.toml')['cases']).items():
        p,d=g['optimized'],g['dag']; mem=d['graph']
        lines.append(f'| {threads} | {key[0]} | {ms(p,"prepared"):.3f} | {ms(d,"prepared"):.3f} | {ms(p,"prepared")/ms(d,"prepared"):.2f}× | {ms(p,"complete_public"):.3f} | {ms(d,"complete_public"):.3f} | {mem["all_slot_bytes"]} | {mem["single_shared_force_bytes"]} | {d["prepared"]["host_bytes"]} |')
lines+=['','## Immutable stage-one versus same-checkout prefix','',
        '| Model | Atoms | Immutable prepared ms | Integrated prefix ms | Immutable public ms | Integrated prefix public ms | Immutable fresh public ms | Integrated fresh prefix ms |',
        '|---|---:|---:|---:|---:|---:|---:|---:|']
sg=group(serial['cases'])
for c in read('immutable-stageone.toml')['cases']:
    p=sg[(c['model'],c['atoms'],c['precision'])]['optimized']
    lines.append(f'| {c["model"]} | {c["atoms"]} | {ms(c,"prepared"):.3f} | {ms(p,"prepared"):.3f} | {ms(c,"complete_public"):.3f} | {ms(p,"complete_public"):.3f} | {ms(c,"construct_and_public"):.3f} | {ms(p,"construct_and_public"):.3f} |')
# Recursively audit every BenchmarkTools record, including isolated stages.
trials=[]
def visit(value,path):
    if isinstance(value,dict):
        if 'raw_s' in value: trials.append(dict(path=path,samples=value['samples'],undersampled=value['undersampled'],min_s=value['min_s'],max_s=value['max_s'],host_bytes=value['host_bytes']))
        else:
            for k,v in value.items(): visit(v,path+'/'+k)
    elif isinstance(value,list):
        for i,v in enumerate(value): visit(v,path+f'/{i}')
for path in root.glob('*.toml'): visit(read(path.name),path.name)
lines+=['',f'{len(trials)} trials audited; {sum(t["undersampled"] for t in trials)} under-sampled. Counts and full extrema are in summary.json and source TOML.','']
(root/'OBSERVATIONS.md').write_text('\n'.join(lines))
(root/'summary.json').write_text(json.dumps(dict(comparisons=summary,trials=trials),indent=2)+'\n')
