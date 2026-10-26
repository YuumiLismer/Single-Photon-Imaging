"""Interleaved multi-arm on/off test (agent, 20 Oct session 2).

    python agent_scripts/multi_ab.py --arm on:27:5000:12 --arm off:15027:5000:12 \
        --arm short:15027:50:6 --acc 50000 --cycles 20 --max-minutes 38 --tag wide_ab

Each --arm is name:delay_ns:width_ns:readouts_per_block[:opt...]. Options (session 6, to diagnose the External
gate-end quirk): I = Internal trigger at --internal-hz, F = falling edge, AC = AC coupling, HZ = high-impedance
termination, t<volts> = trigger threshold (e.g. t0.3), r<hz> = Internal trigger at this rate (e.g. r400000), g<gain> = emICCD gain for this arm (<= 9000),
NC = 'Clean Until Trigger' off, FF = 'Readout Control Mode' Full Frame (default Frame Transfer) for this arm
(set, read back, and restored to the default after the arm). Every cycle takes one block of each arm,
order rotated per cycle so drifts hit all arms alike. Per readout it records photon-hit pixels
(>= --thr ADU above the readout median) in apertures r<=6/10/15 px around --spot, the mean hit
density at r > 50 px (the far field, for the drift-free 'net' score of agent_log E3), the
connected hit clusters ('events') in r<=10 and the readout interval. Records go to
figures/records.json so the test can be re-analysed (e.g. regression on readout interval,
because the ungated background scales with wall-clock time and the arms have different dead time).
Optionally starts with a trigger-rate measurement (--rate).
"""
import argparse, json, sys, pathlib, datetime
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, label
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht

p = argparse.ArgumentParser()
p.add_argument("--arm", action="append", required=True)
p.add_argument("--acc", type=int, default=50_000); p.add_argument("--cycles", type=int, default=20)
p.add_argument("--max-minutes", type=float, default=38); p.add_argument("--gain", type=int, default=9000)
p.add_argument("--spot", type=int, nargs=2, default=(501, 454)); p.add_argument("--thr", type=float, default=200)
p.add_argument("--rate", action="store_true", help="measure the trigger rate first")
p.add_argument("--tag", default="multi_ab")
p.add_argument("--internal-hz", type=float, default=100_000, help="trigger rate of ':I' arms")
p.add_argument("--max-readout-s", type=float, default=1.05,
               help="stop if a block's median readout interval exceeds this (Tristan: readouts <= 1 s wall-clock)")
a = p.parse_args()
arms = []
for s in a.arm:
    n, d, w, k, *opts = s.split(":")
    ov = {}
    for o in opts:
        u = o.upper()
        if u == "I": ov["trigger_source"] = "Internal"
        elif u == "F": ov["trigger_edge"] = "Falling Edge"
        elif u == "AC": ov["trigger_coupling"] = "AC"
        elif u == "HZ": ov["trigger_termination"] = "High Impedance"
        elif u.startswith("T"): ov["trigger_threshold_v"] = float(o[1:])
        elif u.startswith("G"):
            ov["gain"] = int(o[1:])
            if ov["gain"] > ht.MAX_GAIN: sys.exit(f"arm {n}: gain above {ht.MAX_GAIN}")
        elif u == "NC": ov["_attr:Clean Until Trigger"] = False
        elif u == "FF": ov["_attr:Readout Control Mode"] = "Full Frame"
        elif u.startswith("R"): ov["trigger_source"] = "Internal"; ov["internal_trigger_hz"] = float(o[1:])
        else: sys.exit(f"arm {n}: unknown option {o}")
    arms.append(dict(name=n, delay=float(d), width=float(w), readouts=int(k), overrides=ov,
                     source=ov.get("trigger_source", "External")))
    if float(w) * 1e-9 * a.acc > ht.MAX_GATE_OPEN_S:
        sys.exit(f"arm {n}: gate-open time per readout above the safety limit")
spot = tuple(a.spot)
exp = ht.new_experiment(f"agent_{a.tag}"); folder = exp / "figures"
base = dict(gain=a.gain, accumulations=a.acc, trigger_source="External", trigger_threshold_v=0.5,
            roi_center=spot, notes=f"agent {a.tag}")
rec = []; imgs = {m["name"]: [0, 0] for m in arms}


def score(arm, r, cycle):
    x0, y0, n = r["meta"]["roi"]
    yy, xx = np.indices((n, n)); rr = np.hypot(xx - (spot[0] - x0), yy - (spot[1] - y0))
    far = rr > 50
    hit = ht.above_pedestal(r) >= a.thr
    dt = np.diff(np.asarray(r["t"], float))
    dt = np.concatenate([[np.median(dt) if len(dt) else np.nan], dt])   # first readout: use the median
    for i, h in enumerate(hit):
        lab, _ = label(h)
        ev10 = len(np.unique(lab[(rr <= 10) & h]))
        evfar = len(np.unique(lab[far & h]))
        rec.append(dict(arm=arm, cycle=cycle, dt=float(dt[i]), far=float(h[far].mean()), evfar=evfar / far.sum(),
                        h6=int(h[rr <= 6].sum()), h10=int(h[rr <= 10].sum()), h15=int(h[rr <= 15].sum()),
                        ev10=ev10, a6=int((rr <= 6).sum()), a10=int((rr <= 10).sum()), a15=int((rr <= 15).sum())))
    imgs[arm][0] = imgs[arm][0] + hit.sum(0); imgs[arm][1] += len(hit)


def summary():
    out = {}
    for m in arms:
        R = [x for x in rec if x["arm"] == m["name"]]
        if len(R) < 2:
            continue
        net = np.array([x["h10"] - x["a10"] * x["far"] for x in R])
        out[m["name"]] = dict(n=len(R), net10=net.mean(), err=net.std(ddof=1) / np.sqrt(len(R)),
                              raw10=np.mean([x["h10"] for x in R]), dt=np.mean([x["dt"] for x in R]))
    names = [m["name"] for m in arms]
    if len(names) >= 2 and all(k in out for k in names[:2]):
        o, f = out[names[0]], out[names[1]]
        out["z_net_01"] = (o["net10"] - f["net10"]) / np.hypot(o["err"], f["err"])
    return out


ATTR_DEFAULTS = {"Clean Until Trigger": True, "Readout Control Mode": "Frame Transfer"}


def set_attr(name, value):
    """Set a camera attribute (bool or enum label), commit, and check the read-back."""
    import camera_control as cc
    if isinstance(value, bool):
        cc.safe_set_attribute_value(cam, name, value)
    else:
        ok, err = ht._set_enum(cam, cam.get_all_attributes(), name, value)
        if not ok:
            raise RuntimeError(f"could not set {name} = {value}: {err}")
    cc.safe_commit_parameters(cam)
    got = cam.get_attribute_value(name)
    if str(got) != str(value):
        raise RuntimeError(f"{name} read-back {got!r} != {value!r}")
    print(f"  {name} -> {got} (read back)")


cam, limits = ht.open_camera()
t0 = datetime.datetime.now()
try:
    ht.cool_down(cam, setpoint_c=-20.0, poll_s=30)
    if a.rate:
        rate = ht.measure_trigger_rate(cam, base)
        (folder / "trigger_rate.json").write_text(json.dumps(dict(rate_hz=rate, when=str(datetime.datetime.now()))))
    for c in range(a.cycles):
        k = c % len(arms)
        for m in arms[k:] + arms[:k]:
            ov = {k: v for k, v in m["overrides"].items() if not k.startswith("_")}
            extra = {k[6:]: v for k, v in m["overrides"].items() if k.startswith("_attr:")}
            for k, v in extra.items():
                set_attr(k, v)
            try:
                r = ht.take_run(cam, **{**base, "label": f"{a.tag}_{m['name']}_{m['delay']:g}ns_c{c + 1}",
                                        "gate_delay_ns": m["delay"], "gate_width_ns": m["width"], "n_readouts": m["readouts"],
                                        "notes": base["notes"] + "".join(f" | {k} = {v}" for k, v in extra.items()),
                                        **{"internal_trigger_hz": a.internal_hz, **ov}})
            finally:
                for k in extra:
                    set_attr(k, ATTR_DEFAULTS[k])
            if r is None:
                raise RuntimeError("no readouts (triggers stopped?)")
            if float((r["first"] >= ht.SATURATION_ADU).mean()) > 1e-3:
                raise RuntimeError("saturation > 0.1%: stopping")
            score(m["name"], r, c + 1)
            dt_blk = ht.readout_interval_s(r)
            if dt_blk > a.max_readout_s:
                raise RuntimeError(f"readout interval {dt_blk:.2f} s > {a.max_readout_s} s (trigger rate dropped?): "
                                   "stopping; rerun with fewer --acc")
        s = summary()
        print(f"CYCLE {c + 1}: " + json.dumps(ht._jsonable({k: ({kk: round(float(vv), 3) for kk, vv in v.items()}
                                                                 if isinstance(v, dict) else round(float(v), 2)) for k, v in s.items()})), flush=True)
        (folder / "records.json").write_text(json.dumps(dict(args=vars(a), arms=arms, records=rec)))
        if (datetime.datetime.now() - t0).total_seconds() > a.max_minutes * 60:
            print("time limit reached"); break
finally:
    for k in {k[6:] for m in arms for k in m["overrides"] if k.startswith("_attr:")}:
        try:
            set_attr(k, ATTR_DEFAULTS[k])
        except Exception as e:
            print(f"WARNING: could not restore {k}:", e)
    ht.shutdown(cam)
    (folder / "records.json").write_text(json.dumps(dict(args=vars(a), arms=arms, records=rec)))
    done = [m["name"] for m in arms if imgs[m["name"]][1]]
    if len(done) >= 2:
        ims = {k: imgs[k][0] / imgs[k][1] for k in done}
        fig, ax = plt.subplots(1, len(done) + 1, figsize=(5 * (len(done) + 1), 4.6))
        for a_, k in zip(ax, done):
            h = a_.imshow(gaussian_filter(ims[k], 1), cmap="inferno"); a_.set_title(f"{k}: hits/px/readout (1 px smooth)"); fig.colorbar(h, ax=a_)
        h = ax[-1].imshow(gaussian_filter(ims[done[0]] - ims[done[1]], 3), cmap="RdBu_r")
        ax[-1].set_title(f"{done[0]} - {done[1]} (3 px smooth)"); fig.colorbar(h, ax=ax[-1])
        fig.suptitle(f"{a.tag}: " + json.dumps(ht._jsonable({k: (round(float(v['net10']), 2) if isinstance(v, dict) else round(float(v), 2))
                                                              for k, v in summary().items()})))
        fig.tight_layout(); fig.savefig(folder / "multi_ab.png", dpi=110)
print("EXP", exp)
