#!/usr/bin/env python3
"""Summarise benchmark CSVs using only Python's standard library."""
import csv
from pathlib import Path
import sys
import re

here = Path(__file__).resolve().parent
source = Path(sys.argv[1]) if len(sys.argv) > 1 else here / "results"
target = Path(sys.argv[2]) if len(sys.argv) > 2 else here / "RESULTS.md"
rows = []
for path in sorted(source.glob("threads-*.csv")):
    with path.open() as handle:
        for row in csv.DictReader(handle):
            for key in ("threads", "N", "width", "phase_count"):
                row[key] = int(row[key])
            row["parallel"] = row["parallel"] == "true"
            rows.append(row)
if not rows:
    raise SystemExit("No benchmark CSVs found")
index = {}
for r in rows:
    key = (r["N"], r["threads"], r["implementation"], r["regime"], r["width"], r["parallel"])
    if key in index:
        raise SystemExit(f"Duplicate case {key}")
    index[key] = r
sizes = sorted({r["N"] for r in rows})
threads = sorted({r["threads"] for r in rows})
def get(n, t, impl, regime, width, parallel=True):
    return index[n, t, impl, regime, width, parallel]
def ms(r):
    return float(r["median_ms"])
def span(r):
    return f'{ms(r):.3f} [{float(r["min_ms"]):.3f}, {float(r["max_ms"]):.3f}]'

metadata_files=sorted(source.glob("metadata-threads-*.txt"))
metadata=metadata_files[0].read_text().splitlines()[0] if metadata_files else "Hardware metadata unavailable"
warmups=re.search(r"warmups=(\d+)",metadata)
warmup_text=warmups[1] if warmups else "recorded"
sample_text=", ".join(sorted({r["samples"] for r in rows}))
lines = ["# Recorded Argon Lennard-Jones comparison", "",
         f"Recorded environment: `{metadata}`", "",
         "Full calls use `AtomsCalculators.energy_forces_virial` with the matched inverse-distance cutoff tilt.",
         "Reference: EmpiricalPotentials 0.2.6; systems: AtomsSystems 0.1.2 perturbed periodic FCC Ar.",
         f"Worker processes ran sequentially, with {warmup_text} warmups and {sample_text} measured samples per case.",
         "CPU affinity was uncontrolled. Exact versions/hardware are in the metadata files.",
         "All values below are milliseconds; brackets show minimum/maximum around the median.",
         "`EP/PFW4` compares implementations; `W4/W8` compares vector widths (>1 favours W=8).",
         "Ratios use measured medians, not universal speedup claims or confidence intervals.", ""]
for regime, heading in [("rebuild", "Complete public calls, rebuilding topology each call"),
                         ("fixed", "Unchanged-geometry diagnostics, setup excluded")]:
    lines += [f"## {heading}", "",
              "| N | Threads | EP threaded median [min, max] | PFW W=4 phased median [min, max] | EP/PFW4 | PFW W=8 phased median [min, max] | W4/W8 |",
              "|---:|---:|---:|---:|---:|---:|---:|"]
    for n in sizes:
        for t in threads:
            ep = get(n, t, "EmpiricalPotentials", regime, 0)
            pf = get(n, t, "PotentialFrameWork", regime, 4)
            p8 = get(n, t, "PotentialFrameWork", regime, 8)
            lines.append(f"| {n} | {t} | {span(ep)} | {span(pf)} | {ms(ep)/ms(pf):.2f} | {span(p8)} | {ms(pf)/ms(p8):.3f} |")
    lines += [""]

base_thread = min(threads)
lines += ["## Serial vector-width scaling", "",
          f"These rows use the {base_thread}-thread process; the executor is serial.",
          "W=1 uses scalar pair arithmetic; W=4/8 use explicit SIMD vectors (Float64).",
          "W1/W4 and W1/W8 compare against the scalar path; values >1 favour SIMD.", "",
          "| N | Regime | W=1 median | W=4 median | W=8 median | W1/W4 | W1/W8 | W4/W8 |",
          "|---:|---|---:|---:|---:|---:|---:|---:|"]
for n in sizes:
    for regime in ("rebuild","fixed"):
        p1,p4,p8=[get(n,base_thread,"PotentialFrameWork",regime,w,False) for w in (1,4,8)]
        lines.append(f"| {n} | {regime} | {ms(p1):.3f} | {ms(p4):.3f} | {ms(p8):.3f} | {ms(p1)/ms(p4):.3f} | {ms(p1)/ms(p8):.3f} | {ms(p4)/ms(p8):.3f} |")

n = max(sizes)
lines += ["", f"## Thread scaling at N={n}", "",
          f"Each ratio compares against the same implementation/regime at {base_thread} thread(s).", "",
          "| Threads | EP rebuild | PFW W=4 rebuild | PFW W=8 rebuild | EP fixed | PFW W=4 fixed | PFW W=8 fixed |",
          "|---:|---:|---:|---:|---:|---:|---:|"]
for t in threads:
    ratios = []
    for regime in ("rebuild", "fixed"):
        for impl,w in [("EmpiricalPotentials",0),("PotentialFrameWork",4),("PotentialFrameWork",8)]:
            ratios.append(ms(get(n,base_thread,impl,regime,w))/ms(get(n,t,impl,regime,w)))
    lines.append(f'| {t} | ' + ' | '.join(f'{x:.2f}' for x in ratios) + ' |')

lines += ["", f"## Allocated bytes per complete call at N={n}", "",
          "| Threads | EP threaded allocated bytes | PFW W=4 allocated bytes | PFW sorted force bytes | PFW phases |",
          "|---:|---:|---:|---:|---:|"]
for t in threads:
    ep=get(n,t,"EmpiricalPotentials","rebuild",0)
    pf=get(n,t,"PotentialFrameWork","rebuild",4)
    lines.append(f'| {t} | {float(ep["allocated_bytes_median"]):,.0f} | {float(pf["allocated_bytes_median"]):,.0f} | {int(pf["sorted_force_bytes"]):,} | {pf["phase_count"]} |')

lines += ["", "## Numerical agreement", "",
          f"All {len(rows)} measured cases passed energy, force and virial checks before/after timing.",
          "Largest recorded absolute differences across cases:", ""]
for name,column,unit in [("energy","energy_abs_error_eV","eV"),
                         ("force component","force_max_abs_error_eV_per_A","eV/Å"),
                         ("virial component","virial_max_abs_error_eV","eV")]:
    lines.append(f'- {name}: {max(float(r[column]) for r in rows):.3e} {unit}')
lines += ["", "This compares complete implementations: a single-species analytic callback vs a general LJ implementation,",
          "directed vs unordered pair ownership, automatic vs analytic",
          "derivatives/precomputed cutoff coefficients, list construction, allocation, scheduling and reduction.",
          "The fixed-topology table excludes setup, requires unchanged coordinates/cell/species, and is not",
          "a public Verlet reuse claim. Small systems and phase occupancy can limit thread scaling.",
          "Allocated bytes measure Julia allocations, not peak RSS; reference site assembly uses full force",
          "partials per task, whereas PFW's sorted force workspace is independent of thread count.", "",
          "See [methodology](README.md), [CSV records](results/), [time scaling](scaling.svg),",
          "[vector-width scaling](vector-width-scaling.svg), [native-code inspection](codegen/README.md) and per-thread metadata."]
target.write_text("\n".join(lines)+"\n")
print(target)
