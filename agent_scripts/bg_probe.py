"""Where does the wall-clock-proportional background come from? (agent, session 2)

    python agent_scripts/bg_probe.py --cycles 6 --tag bg_probe

Interleaves short readout blocks of different camera states, all with an off-delay gate (no
heralded light), and reports photon-hit pixels (>= 200 ADU above the readout median) over the
150x150 ROI per readout and per second of wall time. States (see agent_log session 2):
  ext      External trigger (idler), gate 15027 ns + 50 ns, 50k acc (the normal hunting state)
  int      Internal trigger 1 kHz, same gate, 1k acc  -> same ~1 s wall time, 50x fewer gates
  bracket  as 'ext' but with MCP bracket gating on (MCP pulsed with the gate), if the camera has it
  int_bracket, int_open, int_open_bracket (session 3): internal-trigger states without/with bracket gating,
           closed (50 ns at 15 us) and open (10 us x 10k = 0.1 s at 10 kHz)
If the hits/s are equal for ext and int, the background does not come from the triggers/gates. If
bracket gating removes most of it, the background is MCP dark events or light passing the gated-off
photocathode, and bracket gating should be used for the hunt.
"""
import argparse, json, sys, pathlib, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht
import camera_control as cc

p = argparse.ArgumentParser()
p.add_argument("--cycles", type=int, default=6); p.add_argument("--readouts", type=int, default=8)
p.add_argument("--states", default="ext,int,bracket"); p.add_argument("--tag", default="bg_probe")
p.add_argument("--spot", type=int, nargs=2, default=(501, 454))
a = p.parse_args()
exp = ht.new_experiment(f"agent_{a.tag}"); folder = exp / "figures"
base = dict(gain=9000, gate_delay_ns=15027.0, gate_width_ns=50.0, accumulations=50_000, n_readouts=a.readouts,
            trigger_source="External", trigger_threshold_v=0.5, roi_center=(500, 455), notes=f"agent {a.tag}")
INT = dict(trigger_source="Internal", internal_trigger_hz=1000, accumulations=1000)
STATES = dict(ext=dict(), int=INT, bracket=dict(), int_bracket=INT,
              # session 3: internal 10 kHz, 10 us gate x 10k (0.1 s gate-open, 1 s wall) with/without bracket
              int_open=dict(INT, internal_trigger_hz=10000, accumulations=10000, gate_delay_ns=27.0, gate_width_ns=10000.0),
              int_open_bracket=dict(INT, internal_trigger_hz=10000, accumulations=10000, gate_delay_ns=27.0,
                                    gate_width_ns=10000.0))
states = a.states.split(",")
rec = []


def set_bracket(cam, on):
    cc.safe_set_attribute_value(cam, "Bracket Gating", bool(on)); cc.safe_commit_parameters(cam)
    got = cam.get_attribute_value("Bracket Gating")
    if bool(got) != bool(on):
        raise RuntimeError(f"Bracket Gating read back {got!r}, wanted {on}")


cam, limits = ht.open_camera()
try:
    attrs = cam.get_all_attributes()
    info = {}
    for n in ("Intensifier Options", "Bracket Gating", "Gating Speed", "Phosphor Type", "Intensifier Diameter",
              "Enable Modulation", "Photocathode Sensitivity", "Phosphor Decay Delay", "Phosphor Decay Delay Resolution"):
        if n in attrs:
            try:
                info[n] = ht._jsonable(cam.get_attribute_value(n))
            except Exception as e:
                info[n] = f"<{e}>"
            try:
                info[n + " (allowed)"] = ht._jsonable(getattr(attrs[n], "values", None) or getattr(attrs[n], "labels", None))
            except Exception:
                pass
    print("INTENSIFIER", json.dumps(info))
    (folder / "intensifier_info.json").write_text(json.dumps(info, indent=1))
    if "Bracket Gating" not in attrs:
        print("no Bracket Gating attribute: skipping bracket states")
        states = [st for st in states if "bracket" not in st]
    t = cam.get_attribute_value("Sensor Temperature Status")   # read only; never touch the set point
    print("sensor", cam.get_attribute_value("Sensor Temperature Reading"), t)
    for c in range(a.cycles):
        k = c % len(states)
        for st in states[k:] + states[:k]:
            if "Bracket Gating" in attrs:
                set_bracket(cam, "bracket" in st)
            r = ht.take_run(cam, **{**base, **STATES[st], "label": f"{a.tag}_{st}_c{c + 1}"})
            if r is None:
                raise RuntimeError("no readouts")
            if float((r["first"] >= ht.SATURATION_ADU).mean()) > 1e-3:
                raise RuntimeError("saturation > 0.1%: stopping")
            x0, y0, n = r["meta"]["roi"]
            yy, xx = np.indices((n, n)); rr = np.hypot(xx - (a.spot[0] - x0), yy - (a.spot[1] - y0))
            hit = ht.above_pedestal(r) >= 200
            dt = ht.readout_interval_s(r)
            rec.append(dict(state=st, cycle=c + 1, dt=dt, roi_hits=float(hit.sum((1, 2)).mean()),
                            spot10=float(hit[:, rr <= 10].sum(1).mean()),
                            ped=float(np.median(r["roi_stack"].astype(float))), n=len(hit)))
            print(f"REC {json.dumps(rec[-1])}", flush=True)
        (folder / "records.json").write_text(json.dumps(dict(info=info, records=rec), indent=1))
finally:
    try:
        if "Bracket Gating" in cam.get_all_attributes():
            set_bracket(cam, False); print("  Bracket Gating -> False")
    except Exception as e:
        print("  bracket reset:", e)
    ht.shutdown(cam)
for st in states:
    R = [x for x in rec if x["state"] == st]
    if R:
        h = np.array([x["roi_hits"] for x in R]); t = np.array([x["dt"] for x in R])
        print(f"SUMMARY {st}: {len(R)} blocks, {h.mean():.1f} ROI hits/readout, {t.mean():.2f} s/readout, "
              f"{(h / t).mean():.1f} +- {(h / t).std(ddof=1) / np.sqrt(len(R)) if len(R) > 1 else 0:.1f} hits/s")
print("EXP", exp)
