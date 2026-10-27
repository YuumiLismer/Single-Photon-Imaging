"""Heralded ON/OFF test for a wide, unfocused spot (no lens; 27 Oct 2026). Saves a large ROI (--roi px square)
around --spot each readout and scores photon hits with the diode-spot image as the matched filter, minus the
same weights applied to the ROI's outer edge region. One cycle of each arm per call (re-lock the herald rate
between calls); --analyse pools several experiment folders.

    python agent_scripts/disc_ab.py --spot 298 277 --acc 20000 --block 25 --tag disc_c1 \
        --arm on31:546:31 --arm off31:1546:31 --arm onwide:450:300 --arm offwide:1450:300
    python agent_scripts/disc_ab.py --analyse EXP [EXP ...] --diode DIODE_EXP --spot 298 277
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
p.add_argument("--spot", type=int, nargs=2, required=True)
p.add_argument("--roi", type=int, default=400)
p.add_argument("--arm", action="append", default=[], help="name:delay_ns:width_ns")
p.add_argument("--acc", type=int, default=20000)
p.add_argument("--block", type=int, default=25)
p.add_argument("--thr", type=float, default=200)
p.add_argument("--tag", default="disc")
p.add_argument("--analyse", nargs="+", default=None)
p.add_argument("--diode", default=None, help="diode_locate experiment: its frame is the matched filter")
a = p.parse_args()

if not a.analyse:
    exp = ht.new_experiment(f"agent_{a.tag}")
    base = dict(gain=9000, trigger_source="External", trigger_threshold_v=0.5, roi_center=tuple(a.spot),
                roi_size=a.roi, kind="signal", accumulations=a.acc, n_readouts=a.block)
    cam, _ = ht.open_camera()
    try:
        ht.cool_down(cam, setpoint_c=-20.0, poll_s=30)
        for s in a.arm:
            n, d, w = s.split(":")
            if float(w) * 1e-9 * a.acc > ht.MAX_GATE_OPEN_S:
                sys.exit(f"arm {n}: gate-open time per readout above the safety limit")
            r = ht.take_run(cam, **{**base, "label": f"{a.tag}_{n}", "gate_delay_ns": float(d), "gate_width_ns": float(w)})
            if r is None:
                raise RuntimeError("no readouts")
            if float((r["first"] >= ht.SATURATION_ADU).mean()) > 1e-3:
                raise RuntimeError("saturation > 0.1%: stopping")
            print(f"  {n} {d}+{w} ns: {ht.readout_interval_s(r):.2f} s/readout", flush=True)
    finally:
        ht.shutdown(cam)
    print("EXP", exp)
    sys.exit()

# ---------------- analysis ----------------
runs = [r for f in a.analyse for r in ht.load_runs(f)]
roi = runs[0]["meta"]["roi"]; x0, y0 = roi[0], roi[1]
H, W = runs[0]["roi_stack"].shape[1:]
yy, xx = np.mgrid[y0:y0 + H, x0:x0 + W]
if a.diode:
    fr = np.load(pathlib.Path(a.diode) / "figures" / "diode_frame.npy").astype(float)
    t = gaussian_filter(fr - np.median(fr), 4)[y0:y0 + H, x0:x0 + W]
    t = np.clip(t, 0, None); t[t < 0.05 * t.max()] = 0
else:
    t = np.exp(-((xx - a.spot[0]) ** 2 + (yy - a.spot[1]) ** 2) / (2 * 60.0 ** 2))
edge = (t == 0) & ((np.abs(xx - a.spot[0]) > 0.4 * W) | (np.abs(yy - a.spot[1]) > 0.4 * H))
arms = {}
for r in runs:
    n = r["meta"]["settings"]["label"].split("_")[-1]
    st = ht.above_pedestal(r); hit = st >= a.thr
    ring = hit[:, edge].mean(1)                                                 # hit density far from the spot
    # template-weighted hits above the edge level, in units of spot photons (exact for a spot of shape t)
    score = ((hit * t).sum((1, 2)) - ring * t.sum()) * t.sum() / (t * t).sum()
    acc = r["meta"]["actual"]["accumulations"]
    d = arms.setdefault(n, dict(score=[], hits=0, n=0, gps=[], delay=r["meta"]["actual"]["gate_delay_ns"],
                               width=r["meta"]["actual"]["gate_width_ns"]))
    d["score"].append(score); d["hits"] = d["hits"] + hit.sum(0); d["n"] += len(hit)
    d["gps"].append(acc / ht.readout_interval_s(r))
res = {}
for on, off in (("on31", "off31"), ("onwide", "offwide")):
    if on not in arms or off not in arms:
        continue
    S, B = np.concatenate(arms[on]["score"]), np.concatenate(arms[off]["score"])
    net = S.mean() - B.mean(); err = np.hypot(S.std(ddof=1) / np.sqrt(len(S)), B.std(ddof=1) / np.sqrt(len(B)))
    acc = runs[0]["meta"]["actual"]["accumulations"]
    res[on] = dict(gate=f"{arms[on]['delay']:.0f}+{arms[on]['width']:.0f}", net_photons_per_readout=net, err=err, z=net / err,
                   per_trigger=net / acc, readouts=[len(S), len(B)], gates_per_s=float(np.median(arms[on]["gps"])))
    print(f"{on:7s} {res[on]['gate']:>8s} ns: net {net:+.2f} +- {err:.2f} spot photons/readout, z {net / err:+.1f}, "
          f"{net / acc:.1e}/trigger | {len(S)} ON + {len(B)} OFF readouts, {res[on]['gates_per_s']:,.0f} gates/s")
out = pathlib.Path(a.analyse[0]) / "figures"; out.mkdir(parents=True, exist_ok=True)
(out / "disc_ab.json").write_text(json.dumps(res, indent=1))
pairs = [(on, off) for on, off in (("on31", "off31"), ("onwide", "offwide")) if on in arms and off in arms]
fig, ax = plt.subplots(1, len(pairs) + 1, figsize=(6 * (len(pairs) + 1), 5.4), squeeze=False); ax = ax[0]
h = ax[0].imshow(t, cmap="magma", extent=[x0, x0 + W, y0 + H, y0]); ax[0].set_title("matched filter (diode spot)")
plt.colorbar(h, ax=ax[0])
for i, (on, off) in enumerate(pairs, 1):
    m = gaussian_filter(arms[on]["hits"] / arms[on]["n"] - arms[off]["hits"] / arms[off]["n"], 12)
    v = np.abs(m).max()
    h = ax[i].imshow(m, cmap="RdBu_r", vmin=-v, vmax=v, extent=[x0, x0 + W, y0 + H, y0])
    ax[i].contour(xx, yy, t, [0.5 * t.max()], colors="k", linewidths=0.6)
    r_ = res[on]
    ax[i].set_title(f"{on} - {off} ({r_['gate']} ns), 12 px smooth\nnet {r_['net_photons_per_readout']:+.2f} +- {r_['err']:.2f}/readout, z {r_['z']:+.1f}")
    plt.colorbar(h, ax=ax[i], label="hits/px/readout")
fig.suptitle(f"Heralded photons in the unfocused spot at {tuple(a.spot)}; black contour = diode half-max")
fig.tight_layout(); fig.savefig(out / "disc_ab.png", dpi=100)
print(out / "disc_ab.png")
