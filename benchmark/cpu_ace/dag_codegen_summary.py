"""Inspect the saved concrete callee LLVM; static counts are not runtime shares."""
import json, pathlib, re, sys
root=pathlib.Path(sys.argv[1]); records=[]
for p in sorted((root/'callee-codegen').glob('*.llvm.txt')):
    s=p.read_text()
    row=dict(file=p.name,
        float_ops=len(re.findall(r'f(?:mul|add) (?:float|<\d+ x float>)',s)),
        double_ops=len(re.findall(r'f(?:mul|add) (?:double|<\d+ x double>)',s)),
        float_to_double=len(re.findall(r'fpext[^\n]*float[^\n]*double',s)),
        direct_gc_allocator_calls=len(re.findall(r'call[^\n]*(?:gc_pool_alloc|gc_alloc|julia.gc_alloc)',s)),
        bounds_error_calls=len(re.findall(r'call[^\n]*throw_boundserror',s)))
    assert not row['float_to_double']
    assert not row['direct_gc_allocator_calls']
    if 'Float32' in p.name: assert not row['double_ops']
    records.append(row)
(root/'codegen-summary.json').write_text(json.dumps(records,indent=2)+'\n')
