"""Final heralded single-photon spot image from disc_ab.py ON/OFF experiments: pooled ON - OFF photon-hit map
(ROI and zoom), radial profile against the diode spot, and second-moment widths (roundness).

    python agent_scripts/spot_image.py EXP [EXP ...] --spot 513 454 --diode DIODE_EXP [--thr 30]
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
p.add_argument("--spot", type=int, nargs=2, required=True)
p.add_argument("--diode", default=None)
p.add_argument("--thr", type=float, default=30)
p.add_argument("--on", default="on31"); p.add_argument("--off", default="off31")
p.add_argument("--zoom", type=int, default=30)
a = p.parse_args()

hon = hoff = 0; non = noff = 0; roi = None
for f in a.exps:
    for r in ht.load_runs(f):
        arm = r["meta"]["settings"]["label"].split("_")[-1]
        if arm not in (a.on, a.off):
            continue
        roi = r["meta"]["roi"]
        hit = ht.above_pedestal(r) >= a.thr
        if arm == a.on:
            hon = hon + hit.sum(0); non += len(hit)
        else:
            hoff = hoff + hit.sum(0); noff += len(hit)
diff = hon / non - hoff / noff
H, W = diff.shape
yy, xx = np.mgrid[roi[1]:roi[1] + H, roi[0]:roi[0] + W]
d = np.hypot(xx - a.spot[0], yy - a.spot[1])
diff = diff - diff[(d > 35) & (d < 65)].mean()                 # remove any residual uniform offset
# second moments of the positive excess within 20 px
m = np.clip(gaussian_filter(diff, 1), 0, None) * (d <= 20)
cx, cy = (m * xx).sum() / m.sum(), (m * yy).sum() / m.sum()
sx = np.sqrt((m * (xx - cx) ** 2).sum() / m.sum()); sy = np.sqrt((m * (yy - cy) ** 2).sum() / m.sum())
sxy = (m * (xx - cx) * (yy - cy)).sum() / m.sum()
ev = np.linalg.eigvalsh([[sx ** 2, sxy], [sxy, sy ** 2]])
ell = float(np.sqrt(ev[0] / ev[1]))
rb = np.arange(0, 41, 2)
prof = [diff[(d >= r0) & (d < r0 + 2)].mean() for r0 in rb[:-1]]
perr = [diff[(d >= r0) & (d < r0 + 2)].std() / np.sqrt(((d >= r0) & (d < r0 + 2)).sum()) for r0 in rb[:-1]]
res = dict(readouts_on=non, readouts_off=noff, thr=a.thr, centroid=[round(float(cx), 1), round(float(cy), 1)],
           sigma_xy_px=[round(float(sx), 2), round(float(sy), 2)], axis_ratio=round(ell, 2),
           net_hits_r10_per_readout=float(diff[d <= 10].sum()))
fig = plt.figure(figsize=(17, 5.6))
ax1 = fig.add_subplot(1, 3, 1); ax2 = fig.add_subplot(1, 3, 2); ax3 = fig.add_subplot(1, 3, 3)
z = a.zoom
sl = (slice(a.spot[1] - roi[1] - z, a.spot[1] - roi[1] + z), slice(a.spot[0] - roi[0] - z, a.spot[0] - roi[0] + z))
ext = [a.spot[0] - z, a.spot[0] + z, a.spot[1] + z, a.spot[1] - z]
for ax, s, t in ((ax1, 1.0, "1 px smooth"), (ax2, 2.0, "2 px smooth")):
    im = gaussian_filter(diff, s)[sl]
    h = ax.imshow(im, cmap="inferno", vmin=0, vmax=im.max(), extent=ext, interpolation="nearest")
    ax.set_title(f"heralded photons (ON - OFF), {t}"); ax.set_xlabel("x (px)"); ax.set_ylabel("y (px)")
    plt.colorbar(h, ax=ax, label="photon hits / px / readout")
ax3.errorbar(rb[:-1] + 1, prof, perr, fmt="o-", color="C3", label="heralded photons (ON - OFF)")
if a.diode:
    fr = np.load(pathlib.Path(a.diode) / "figures" / "diode_frame.npy").astype(float)
    fr = fr - np.median(fr)
    Y, X = np.mgrid[:fr.shape[0], :fr.shape[1]]
    dd = np.hypot(X - a.spot[0], Y - a.spot[1])
    dp = np.array([fr[(dd >= r0) & (dd < r0 + 2)].mean() for r0 in rb[:-1]])
    ax3.plot(rb[:-1] + 1, dp / dp.max() * max(prof), "k--", label="diode spot (scaled)")
ax3.axhline(0, color="gray", lw=0.5); ax3.set_xlabel("radius from spot centre (px)")
ax3.set_ylabel("hits / px / readout"); ax3.legend(); ax3.set_title(f"radial profile; sigma x/y {sx:.1f}/{sy:.1f} px, axis ratio {ell:.2f}")
fig.suptitle(f"Heralded single-photon spot, PI-MAX4, gate 546 + 31 ns, threshold {a.thr:g} ADU: {non} ON + {noff} OFF readouts, "
             f"centroid ({cx:.1f}, {cy:.1f})")
fig.tight_layout()
out = pathlib.Path(a.exps[0]) / "figures"; out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "spot_image.png", dpi=120)
(out / "spot_image.json").write_text(json.dumps(res, indent=1))
print(json.dumps(res)); print(out / "spot_image.png")
