#!/usr/bin/env python3
"""Optional standalone plots; requires matplotlib, outside Julia timing runs."""
import csv
from pathlib import Path
import sys
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here=Path(__file__).resolve().parent
source=Path(sys.argv[1]) if len(sys.argv)>1 else here/"results"
rows=[]
for path in source.glob("threads-*.csv"):
    with path.open() as f:
        rows.extend(csv.DictReader(f))
sizes=sorted({int(r["N"]) for r in rows})
threads=sorted({int(r["threads"]) for r in rows})
if len(sizes)!=4:
    raise SystemExit("This four-panel plot expects four system sizes; tables support arbitrary sizes")
fig,axes=plt.subplots(2,2,figsize=(10,7),sharex=True)
for ax,n in zip(axes.flat,sizes):
    for impl,width,color in [("EmpiricalPotentials",0,"#bd5d19"),("PotentialFrameWork",4,"#1464a2"),("PotentialFrameWork",8,"#287e40")]:
        for regime,style in [("rebuild","-"),("fixed","--")]:
            points=sorted([r for r in rows if int(r["N"])==n and r["implementation"]==impl and
                           int(r["width"])==width and r["regime"]==regime and r["parallel"]=="true"],
                          key=lambda r:int(r["threads"]))
            xs=[int(r["threads"]) for r in points]
            ys=[float(r["median_ms"]) for r in points]
            label=("EP" if width==0 else f"PFW W={width}")+" "+regime
            ax.plot(xs,ys,marker="s" if width==8 else "o",linestyle=style,color=color,label=label,linewidth=1.8,markersize=4)
            if regime=="rebuild":
                ax.fill_between(xs,[float(r["min_ms"]) for r in points],
                                [float(r["max_ms"]) for r in points],color=color,alpha=0.12)
    ax.set_title(f"{n:,} Argon atoms")
    ax.set_yscale("log")
    ax.set_xticks(threads)
    ax.set_ylabel("Milliseconds (log scale)")
    ax.set_xlabel("Julia threads")
    ax.grid(True,which="both",alpha=0.2)
fig.suptitle("Argon Lennard-Jones: complete calls and fixed-topology diagnostics",fontsize=13)
handles,labels=axes[0,0].get_legend_handles_labels()
fig.legend(handles,labels,loc="lower center",ncol=3,bbox_to_anchor=(0.5,0.025),frameon=False)
meta_files=sorted(source.glob("metadata-threads-*.txt"))
meta=meta_files[0].read_text().splitlines()[0] if meta_files else ""
cpu=re.search(r"CPU=(.*?) compiler_target=",meta)
julia=re.search(r"Julia=([^ ]+)",meta)
samples=", ".join(sorted({r["samples"] for r in rows}))
caption=f"Median of {samples} samples; shading=min/max. {cpu[1] if cpu else 'See metadata'}. Julia {julia[1] if julia else '?'}. Fixed excludes setup."
fig.text(0.5,0.005,caption,ha="center",fontsize=7)
fig.tight_layout(rect=(0,0.09,1,0.94))
for suffix in ("svg","png"):
    fig.savefig(here/f"scaling.{suffix}",dpi=180)
print(here/"scaling.svg")

# Isolate W=4 versus W=8 so nearly overlapping absolute-time curves remain clear.
fig2,axes2=plt.subplots(2,2,figsize=(10,7),sharex=True,sharey=True)
for ax,n in zip(axes2.flat,sizes):
    for regime,color,marker in [("rebuild","#bd5d19","o"),("fixed","#1464a2","s")]:
        ratios=[]
        for t in threads:
            def measured(width):
                matches=[r for r in rows if int(r["N"])==n and int(r["threads"])==t and
                         r["implementation"]=="PotentialFrameWork" and int(r["width"])==width and
                         r["regime"]==regime and r["parallel"]=="true"]
                if len(matches)!=1:
                    raise ValueError(f"Expected one W={width} case at N={n}, T={t}, {regime}")
                return float(matches[0]["median_ms"])
            ratios.append(measured(4)/measured(8))
        ax.plot(threads,ratios,label=regime,color=color,marker=marker,linewidth=1.8)
    ax.axhline(1.0,color="0.4",linestyle="--",linewidth=1)
    ax.set_title(f"{n:,} Argon atoms")
    ax.set_xticks(threads)
    ax.set_xlabel("Julia threads")
    ax.set_ylabel("W=4 time / W=8 time")
    ax.grid(True,alpha=0.2)
fig2.suptitle("Vector-width scaling: above 1 favours W=8; below 1 favours W=4",fontsize=12)
h,l=axes2[0,0].get_legend_handles_labels()
fig2.legend(h,l,loc="lower center",ncol=2,bbox_to_anchor=(0.5,0.025),frameon=False)
ratio_caption=f"Ratios of medians ({samples} samples); no confidence intervals. {cpu[1] if cpu else 'See metadata'}. Julia {julia[1] if julia else '?'}."
fig2.text(0.5,0.005,ratio_caption,ha="center",fontsize=7)
fig2.tight_layout(rect=(0,0.09,1,0.94))
for suffix in ("svg","png"):
    fig2.savefig(here/f"vector-width-scaling.{suffix}",dpi=180)
print(here/"vector-width-scaling.svg")
