"""Where (and in which gate) are the heralded photons? For disc_ab.py ladder experiments over a large ROI: per arm,
the mean photon-hit map minus its own ROI mean (removes the readout-time background), minus the median of the other
arms' maps. Reports the brightest smoothed excess per arm (position, significance vs the map's own noise) and the
net hits in a disc there, and saves ladder_map.png.

    python agent_scripts/ladder_map.py EXP [EXP ...] [--thr 30] [--smooth 4] [--r 12]
"""
import argparse
import json
import pathlib
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht

p = argparse.ArgumentParser()
p.add_argument("exps", nargs="+")
p.add_argument("--thr", type=float, default=30)
p.add_argument("--smooth", type=float, default=4)
p.add_argument("--r", type=float, default=12)
a = p.parse_args()

arms, roi = {}, None
for f in a.exps:
    for r in ht.load_runs(f):
        roi = r["meta"]["roi"]
        n = r["meta"]["settings"]["label"].split("_")[-1]
        act = r["meta"]["actual"]
        hit = ht.above_pedestal(r) >= a.thr
        v = arms.setdefault(n, dict(sum=0, n=0, delay=act["gate_delay_ns"], width=act["gate_width_ns"]))
        v["sum"] = v["sum"] + hit.sum(0); v["n"] += len(hit)
names = sorted(arms, key=lambda n: arms[n]["delay"])
maps = {n: arms[n]["sum"] / arms[n]["n"] for n in names}
maps = {n: m - m.mean() for n, m in maps.items()}
H, W = maps[names[0]].shape
yy, xx = np.mgrid[roi[1]:roi[1] + H, roi[0]:roi[0] + W]
edge = 3 * a.smooth
res = {}
fig, ax = plt.subplots(1, len(names), figsize=(5 * len(names), 5), squeeze=False); ax = ax[0]
for i, n in enumerate(names):
    others = [maps[o] for o in names if o != n]
    diff = maps[n] - (np.median(others, axis=0) if len(others) > 1 else others[0])
    g = gaussian_filter(diff, a.smooth)
    inner = g[int(edge):-int(edge), int(edge):-int(edge)]
    sd = inner.std()
    iy, ix = np.unravel_index(np.argmax(inner), inner.shape); iy += int(edge); ix += int(edge)
    px, py = int(xx[iy, ix]), int(yy[iy, ix])
    d = np.hypot(xx - px, yy - py)
    net = float(diff[d <= a.r].sum() - diff[(d > 3 * a.r) & (d < 5 * a.r)].mean() * (d <= a.r).sum())
    res[n] = dict(gate=f"{arms[n]['delay']:.0f}+{arms[n]['width']:.0f}", peak_xy=[px, py], peak_sigma=float(inner.max() / sd),
                  net_hits_per_readout=net, readouts=arms[n]["n"])
    print(f"{n:>6s} {res[n]['gate']:>12s} ns: brightest excess at ({px}, {py}), {inner.max() / sd:.1f} sd of the map; "
          f"net {net:+.2f} hits/readout in r<={a.r:g} | {arms[n]['n']} readouts", flush=True)
    v = np.abs(g).max()
    h = ax[i].imshow(g, cmap="RdBu_r", vmin=-v, vmax=v, extent=[roi[0], roi[0] + W, roi[1] + H, roi[1]])
    ax[i].add_patch(plt.Circle((px, py), 3 * a.r, fill=False, color="k", lw=0.8))
    ax[i].set_title(f"{n} {res[n]['gate']} ns minus others\npeak ({px},{py}) {inner.max() / sd:.1f} sd")
    plt.colorbar(h, ax=ax[i])
fig.tight_layout()
out = pathlib.Path(a.exps[0]) / "figures"; out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "ladder_map.png", dpi=80); (out / "ladder_map.json").write_text(json.dumps(res, indent=1))
print(out / "ladder_map.png")
