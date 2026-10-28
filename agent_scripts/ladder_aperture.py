"""Heralded-photon arrival search over a ladder of gates (disc_ab.py arms, any names) at a compact spot: photon hits
within --r px of --spot minus the ring level times the aperture area, per readout, for every arm; each arm compared
with the median of all arms (the photons sit in one gate, every other gate is a control). Saves ladder.png/json.

    python agent_scripts/ladder_aperture.py EXP [EXP ...] --spot 513 454 [--thr 30]
"""
import argparse
import json
import pathlib
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht

p = argparse.ArgumentParser()
p.add_argument("exps", nargs="+")
p.add_argument("--spot", type=int, nargs=2, required=True)
p.add_argument("--thr", type=float, default=30)
p.add_argument("--r", type=float, default=10)
p.add_argument("--ring", type=float, nargs=2, default=(35, 65))
a = p.parse_args()

arms = {}
for f in a.exps:
    for r in ht.load_runs(f):
        roi = r["meta"]["roi"]; H, W = r["roi_stack"].shape[1:]
        yy, xx = np.mgrid[roi[1]:roi[1] + H, roi[0]:roi[0] + W]
        d = np.hypot(xx - a.spot[0], yy - a.spot[1]); ap = d <= a.r; ring = (d > a.ring[0]) & (d < a.ring[1])
        hit = ht.above_pedestal(r) >= a.thr
        n = r["meta"]["settings"]["label"].split("_")[-1]
        act = r["meta"]["actual"]
        v = arms.setdefault(n, dict(s=[], ring=[], delay=act["gate_delay_ns"], width=act["gate_width_ns"],
                                    acc=act["accumulations"], dt=[]))
        v["s"].append(hit[:, ap].sum(1) - hit[:, ring].mean(1) * ap.sum())
        v["ring"].append(hit[:, ring].mean(1)); v["dt"].append(ht.readout_interval_s(r))
names = sorted(arms, key=lambda n: arms[n]["delay"])
m = {n: np.concatenate(arms[n]["s"]) for n in names}
mean = np.array([m[n].mean() for n in names]); err = np.array([m[n].std(ddof=1) / np.sqrt(len(m[n])) for n in names])
base = np.median(mean)
res = {}
for n, mu, e in zip(names, mean, err):
    v = arms[n]
    res[n] = dict(gate=f"{v['delay']:.0f}+{v['width']:.0f}", spot=float(mu), err=float(e), z_vs_median=float((mu - base) / e),
                  readouts=len(m[n]), ring_per_px=float(np.concatenate(v["ring"]).mean()), s_per_readout=float(np.median(v["dt"])))
    print(f"{n:>8s} {res[n]['gate']:>12s} ns: spot {mu:+7.2f} +- {e:.2f} hits/readout (z {res[n]['z_vs_median']:+5.1f} vs median) "
          f"| ring {res[n]['ring_per_px']:.2e}/px | {len(m[n])} readouts, {res[n]['s_per_readout']:.2f} s", flush=True)
out = pathlib.Path(a.exps[0]) / "figures"; out.mkdir(parents=True, exist_ok=True)
(out / "ladder.json").write_text(json.dumps(res, indent=1))
fig, ax = plt.subplots(figsize=(9, 4.5))
x = [arms[n]["delay"] + arms[n]["width"] / 2 for n in names]; xe = [arms[n]["width"] / 2 for n in names]
ax.errorbar(x, mean, err, xe, fmt="o", color="C3"); ax.axhline(base, color="gray", lw=0.6)
ax.set_xlabel("gate (ns after herald trigger)"); ax.set_ylabel(f"spot hits/readout (r<={a.r:g} px, minus ring)")
ax.set_title(f"heralded-photon arrival search at {tuple(a.spot)}, thr {a.thr:g} ADU")
fig.tight_layout(); fig.savefig(out / "ladder.png", dpi=110); print(out / "ladder.png")
