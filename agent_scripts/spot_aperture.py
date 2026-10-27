"""Heralded ON-OFF for a compact spot from disc_ab.py experiments (one cycle per folder): photon hits within
--r px of --spot minus the per-pixel level in a ring (--ring r1 r2) times the aperture area, per readout;
ON-OFF paired per cycle and combined with inverse-variance weights; for several hit thresholds. Saves the
ON-OFF hit image (first threshold) in the first folder's figures.

    python agent_scripts/spot_aperture.py --spot 513 454 EXP [EXP ...] [--thr 30 200] [--r 10]
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
p.add_argument("--thr", type=float, nargs="+", default=[30, 200])
p.add_argument("--r", type=float, default=10)
p.add_argument("--ring", type=float, nargs=2, default=(35, 65))
p.add_argument("--on", default="on31"); p.add_argument("--off", default="off31")
a = p.parse_args()

cycles = [ht.load_runs(f) for f in a.exps]
out, img = {}, None
for thr in a.thr:
    per, hon, hoff, non, noff, acc, gps = [], 0, 0, 0, 0, None, []
    for rs in cycles:
        v = {}
        for r in rs:
            roi = r["meta"]["roi"]; H, W = r["roi_stack"].shape[1:]
            yy, xx = np.mgrid[roi[1]:roi[1] + H, roi[0]:roi[0] + W]
            d = np.hypot(xx - a.spot[0], yy - a.spot[1]); ap = d <= a.r; ring = (d > a.ring[0]) & (d < a.ring[1])
            hit = ht.above_pedestal(r) >= thr
            arm = r["meta"]["settings"]["label"].split("_")[-1]
            v[arm] = hit[:, ap].sum(1) - hit[:, ring].mean(1) * ap.sum()
            if arm == a.on:
                hon = hon + hit.sum(0); non += len(hit); acc = r["meta"]["actual"]["accumulations"]
                gps.append(acc / ht.readout_interval_s(r))
            elif arm == a.off:
                hoff = hoff + hit.sum(0); noff += len(hit)
        if a.on in v and a.off in v:
            S, B = v[a.on], v[a.off]
            per.append((S.mean() - B.mean(), np.hypot(S.std(ddof=1) / np.sqrt(len(S)), B.std(ddof=1) / np.sqrt(len(B)))))
    m = np.array([x[0] for x in per]); e = np.array([x[1] for x in per]); w = 1 / e ** 2
    net, err = float((m * w).sum() / w.sum()), float(1 / np.sqrt(w.sum()))
    out[thr] = dict(net_per_readout=net, err=err, z=net / err, per_trigger=net / acc, cycles=len(per), readouts_on=non,
                    gates_per_s=float(np.median(gps)), per_cycle=[round(x, 2) for x in m])
    print(f"thr {thr:4g} ADU: net {net:+.2f} +- {err:.2f} hits/readout in r<={a.r:g} px, z {net / err:+.1f}, "
          f"{net / acc:.1e}/trigger | {len(per)} cycles, {non} ON readouts, ~{np.median(gps):,.0f} gates/s | per cycle "
          + " ".join(f"{x:+.1f}" for x in m), flush=True)
    if img is None:
        img = (hon / non - hoff / noff, thr, roi)
fold = pathlib.Path(a.exps[0]) / "figures"; fold.mkdir(parents=True, exist_ok=True)
(fold / "spot_aperture.json").write_text(json.dumps({str(k): v for k, v in out.items()}, indent=1))
diff, thr, roi = img
H, W = diff.shape
fig, ax = plt.subplots(1, 2, figsize=(12, 5.4))
for a_, s in ((ax[0], 1), (ax[1], 3)):
    m_ = gaussian_filter(diff, s); v_ = np.abs(m_).max()
    h = a_.imshow(m_, cmap="RdBu_r", vmin=-v_, vmax=v_, extent=[roi[0], roi[0] + W, roi[1] + H, roi[1]])
    a_.add_patch(plt.Circle(tuple(a.spot), a.r, fill=False, color="k", lw=0.8))
    a_.set_title(f"ON - OFF photon hits/px/readout, thr {thr:g} ADU, {s} px smooth")
    plt.colorbar(h, ax=a_)
o = out[thr]
fig.suptitle(f"Heralded single photons at {tuple(a.spot)}: net {o['net_per_readout']:+.2f} +- {o['err']:.2f}/readout, "
             f"z {o['z']:+.1f} ({o['cycles']} cycles, {o['readouts_on']} ON readouts)")
fig.tight_layout(); fig.savefig(fold / "spot_aperture.png", dpi=110)
print(fold / "spot_aperture.png")
