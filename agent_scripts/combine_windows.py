"""Combine windows.json from several window_figures.py runs (different accumulations) per trigger
(inverse-variance weights) and plot net heralded spot hits per trigger vs delay (agent, session 4).

    python agent_scripts/combine_windows.py <out folder> <exp folder 1> <exp folder 2> ...
"""
import sys, json, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out = pathlib.Path(sys.argv[1]); rows = {}
for f in sys.argv[2:]:
    for o in json.loads((pathlib.Path(f) / "figures" / "windows.json").read_text()):
        rows.setdefault(o["delay"], []).append(o)
res = []
for d in sorted(rows):
    v = np.array([o["net_spot_per_readout"] / o["acc"] for o in rows[d]]); e = np.array([o["err"] / o["acc"] for o in rows[d]])
    w = 1 / e ** 2; m = (w * v).sum() / w.sum(); s = 1 / np.sqrt(w.sum())
    res.append(dict(delay=d, width=rows[d][0]["width"], per_trigger=m, err=s, z=m / s, limit_2sigma=max(m, 0) + 2 * s,
                    readouts_on=sum(o["readouts_on"] for o in rows[d]), triggers_on=sum(o["readouts_on"] * o["acc"] for o in rows[d])))
    print(f"{d:5d} ns: {m:+.2e} +- {s:.2e} spot hits/trigger, z {m / s:+.2f}, 2-sigma limit {res[-1]['limit_2sigma']:.1e} "
          f"({res[-1]['readouts_on']} ON readouts, {res[-1]['triggers_on'] / 1e6:.1f} M triggers)")
fig, ax = plt.subplots(figsize=(8, 4.5))
x = np.array([r["delay"] + r["width"] / 2 for r in res])
ax.errorbar(x, [r["per_trigger"] * 1e6 for r in res], yerr=[r["err"] * 1e6 for r in res], xerr=[r["width"] / 2 for r in res], fmt="o", capsize=4)
for r, xx in zip(res, x):
    ax.annotate(f"z {r['z']:+.1f}", (xx, (r["per_trigger"] + r["err"]) * 1e6), ha="center", va="bottom", fontsize=9)
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("gate window after the herald trigger (ns)"); ax.set_ylabel("net spot hits per 10^6 triggers (ON - OFF)")
ax.set_title(f"All passes combined, {res[0]['width']:g} ns gates, gain 9000, external trigger")
fig.tight_layout(); fig.savefig(out / "figures" / "summary_vs_delay_combined.png", dpi=110)
(out / "figures" / "windows_combined.json").write_text(json.dumps(res, indent=1))
