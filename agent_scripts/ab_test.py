"""On-delay vs off-delay A/B test (agent, 20 Oct): the decisive '>= 5 sigma above an off-delay
background' check for one candidate window.

    python agent_scripts/ab_test.py --delay 327 --width 500 --off 10327 --acc 50000 --readouts 12 --cycles 25

Alternates on/off runs (same width, accumulations and gain) so drifts hit both alike, and scores
each readout as spot hits minus area-scaled hits at r > 50 px (see agent_log E3), plus raw hits.
Saves the summed on/off images and the difference image.
"""
import argparse, json, sys, pathlib, datetime
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht

p = argparse.ArgumentParser()
p.add_argument("--delay", type=float, required=True); p.add_argument("--width", type=float, required=True)
p.add_argument("--off", type=float, required=True); p.add_argument("--acc", type=int, default=50_000)
p.add_argument("--readouts", type=int, default=12); p.add_argument("--cycles", type=int, default=25)
p.add_argument("--max-minutes", type=float, default=28)
p.add_argument("--spot", type=int, nargs=2, default=(501, 454)); p.add_argument("--thr", type=float, default=200)
p.add_argument("--tag", default="ab")
a = p.parse_args()
spot = tuple(a.spot)
folder = ht.new_experiment(f"agent_{a.tag}") / "figures"   # runs go to the same experiment folder
base = dict(gain=9000, gate_width_ns=a.width, accumulations=a.acc, n_readouts=a.readouts, trigger_source="External",
            trigger_threshold_v=0.5, roi_center=(500, 455), notes=f"agent {a.tag}")
acc = {"on": dict(net=[], raw=[], img=0, n=0, t=0.0), "off": dict(net=[], raw=[], img=0, n=0, t=0.0)}


def add(key, r):
    x0, y0, n = r["meta"]["roi"]
    yy, xx = np.indices((n, n)); rr = np.hypot(xx - (spot[0] - x0), yy - (spot[1] - y0))
    ap, far = rr <= 10, rr > 50
    hit = ht.above_pedestal(r) >= a.thr
    A = acc[key]
    A["net"] += list(hit[:, ap].sum(1) - ap.sum() * hit[:, far].mean(1)); A["raw"] += list(hit[:, ap].sum(1))
    A["img"] = A["img"] + hit.sum(0); A["n"] += len(hit); A["t"] += ht.readout_interval_s(r) * len(hit)


def report():
    on, off = acc["on"], acc["off"]
    z = ht._diff_z(on["net"], off["net"]); zr = ht._diff_z(on["raw"], off["raw"])
    return dict(n_on=on["n"], n_off=off["n"], net_on=np.mean(on["net"]), net_off=np.mean(off["net"]), z_net=z,
                raw_on=np.mean(on["raw"]), raw_off=np.mean(off["raw"]), z_raw=zr,
                s_per_readout_on=on["t"] / on["n"], s_per_readout_off=off["t"] / off["n"])


cam, limits = ht.open_camera()
t0 = datetime.datetime.now()
try:
    ht.cool_down(cam, setpoint_c=-20.0, poll_s=30)
    for c in range(a.cycles):
        for key, d in (("on", a.delay), ("off", a.off)) if c % 2 == 0 else (("off", a.off), ("on", a.delay)):
            r = ht.take_run(cam, **{**base, "label": f"{a.tag}_{key}_{d:g}ns_c{c + 1}", "gate_delay_ns": d})
            if r is None:
                raise RuntimeError("no readouts (triggers stopped?)")
            if float((r["first"] >= ht.SATURATION_ADU).mean()) > 1e-3:
                raise RuntimeError("saturation > 0.1%: stopping")
            add(key, r)
        rep = report()
        print(f"CYCLE {c + 1}: " + json.dumps({k: round(float(v), 3) for k, v in rep.items()}), flush=True)
        (folder / "ab.json").write_text(json.dumps(ht._jsonable(dict(args=vars(a), cycles=c + 1, **rep)), indent=1))
        if (datetime.datetime.now() - t0).total_seconds() > a.max_minutes * 60:
            print("time limit reached"); break
finally:
    ht.shutdown(cam)
    if acc["on"]["n"] and acc["off"]["n"]:
        on, off = acc["on"]["img"] / acc["on"]["n"], acc["off"]["img"] / acc["off"]["n"]
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
        for a_, im, t in ((ax[0], on, f"on {a.delay:g} ns"), (ax[1], off, f"off {a.off:g} ns"),
                          (ax[2], gaussian_filter(on - off, 3), "on - off (3 px smoothing)")):
            h = a_.imshow(im, cmap="RdBu_r" if "-" in t else "inferno"); a_.set_title(t + ", hits per px per readout"); fig.colorbar(h, ax=a_)
        rep = report()
        fig.suptitle(f"{a.tag}: {ht.fmt_ns(a.width)} gate, {a.acc:,} acc, {rep['n_on']}/{rep['n_off']} readouts; "
                     f"net spot hits on {rep['net_on']:.2f} vs off {rep['net_off']:.2f} per readout, z = {rep['z_net']:+.1f}")
        fig.tight_layout(); fig.savefig(folder / "ab_test.png", dpi=110)
print("FIGS", folder)
