"""Plain aperture estimate for the wide spot: hits inside the diode half-max outline minus outside hits scaled by
area, per readout; ON-OFF paired per cycle (same rate), at several thresholds."""
import sys, glob, numpy as np
sys.path.insert(0, r"C:/Users/ReschLab/Single-Photon-Imaging")
import hoi_troubleshoot as ht
from scipy.ndimage import gaussian_filter
D = r"C:/Users/ReschLab/Desktop/HOI_freshstaret/data/troubleshoot/2026-10-27"
fr = np.load(D + "/1446_agent_diode_locate/figures/diode_frame.npy").astype(float)
tf = gaussian_filter(fr - np.median(fr), 4)
cyc = sorted(glob.glob(D + "/*_agent_disc_c[1-8]"))
for thr in (200, 100, 30):
    per = {"on31": [], "onwide": []}; noise = []
    for f in cyc:
        vals = {}
        for r in ht.load_runs(f):
            roi = r["meta"]["roi"]; x0, y0 = roi[0], roi[1]; H, W = r["roi_stack"].shape[1:]
            t = tf[y0:y0 + H, x0:x0 + W]; inside = t > 0.5 * t.max(); outside = t < 0.05 * t.max()
            hit = ht.above_pedestal(r) >= thr
            y = hit[:, inside].sum(1) - hit[:, outside].mean(1) * inside.sum()
            vals[r["meta"]["settings"]["label"].split("_")[-1]] = y
        for on, off in (("on31", "off31"), ("onwide", "offwide")):
            per[on].append((vals[on].mean() - vals[off].mean(), np.hypot(vals[on].std(ddof=1), vals[off].std(ddof=1)) / np.sqrt(len(vals[on]))))
        noise.append(vals["off31"].std(ddof=1))
    for on in per:
        m = np.array([p[0] for p in per[on]]); e = np.array([p[1] for p in per[on]])
        w = 1 / e ** 2; net = (m * w).sum() / w.sum(); err = 1 / np.sqrt(w.sum())
        print(f"thr {thr:3d} {on:7s}: net {net:+7.2f} +- {err:5.2f} photons/readout (z {net/err:+.1f}) | per-readout noise ~{np.mean(noise):.1f} | per-cycle: " + " ".join(f"{v:+.0f}" for v in m))
