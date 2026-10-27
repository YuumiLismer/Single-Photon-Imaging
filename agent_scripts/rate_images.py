"""Side-by-side heralded images (ON - OFF photon hits per readout) for rate-test experiments:
one row per experiment (or merged group), columns 4 ns and 31 ns arms, shared colour scale.
    python agent_scripts/rate_images.py OUT.png "label=EXP[,EXP...]" ["label=EXP..."]"""
import sys, numpy as np
sys.path.insert(0, r"C:/Users/ReschLab/Single-Photon-Imaging")
import hoi_troubleshoot as ht
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
out, rows = sys.argv[1], [a.split("=", 1) for a in sys.argv[2:]]
maps = []
for label, folders in rows:
    rs = [r for f in folders.split(",") for r in ht.load_runs(f)]
    acc = {}
    for r in rs:
        lab = r["meta"]["settings"]["label"]; arm = "off" if "off" in lab else ("w4" if "w4" in lab else "w31")
        hit = (ht.above_pedestal(r) >= 200)
        s, n, gps = acc.get(arm, (0, 0, []))
        acc[arm] = (s + hit.sum(0), n + len(hit), gps + [r["meta"]["actual"]["accumulations"] / ht.readout_interval_s(r)])
    off = acc["off"][0] / acc["off"][1]
    roi = rs[0]["meta"]["roi"]
    maps.append((label, roi, np.median(acc["off"][2]),
                 [(k, gaussian_filter(acc[k][0] / acc[k][1] - off, 2), acc[k][1]) for k in ("w4", "w31")]))
v = max(np.abs(m[75 - 20:75 + 20, 75 - 20:75 + 20]).max() for _, _, _, arms in maps for _, m, _ in arms)
fig, ax = plt.subplots(len(maps), 2, figsize=(10, 4.4 * len(maps)), squeeze=False)
for i, (label, roi, gps, arms) in enumerate(maps):
    ext = [roi[0], roi[0] + 150, roi[1] + 150, roi[1]]
    for j, (k, m, n) in enumerate(arms):
        h = ax[i, j].imshow(m, cmap="RdBu_r", vmin=-v, vmax=v, extent=ext)
        ax[i, j].set_title(f"{label} (~{gps/1e3:.0f} k gates/s): {'547 + 4 ns' if k == 'w4' else '546 + 31 ns'}\nON - OFF, {n} readouts, 2 px smooth")
        ax[i, j].set_xlabel("x (px)"); ax[i, j].set_ylabel("y (px)")
        plt.colorbar(h, ax=ax[i, j], label="heralded hits / px / readout")
fig.suptitle("Heralded spot vs herald rate: narrow (4 ns) and wide (31 ns) gate; OFF = 1546 + 4 ns; same colour scale")
fig.tight_layout(); fig.savefig(out, dpi=110); print(out)
