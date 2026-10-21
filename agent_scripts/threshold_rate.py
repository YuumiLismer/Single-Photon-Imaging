"""External trigger rate vs trigger threshold and edge (agent, session 6).

    python agent_scripts/threshold_rate.py [--thr 0.2 0.3 ...]

Per setting: 3 readouts of 20k accumulations x 3 ns gate (gives the rate if < ~88 kHz); if that is limited by the
0.226 s readout time, 3 more readouts with accumulations sized for ~0.6 s per readout. Rate = accumulations /
readout interval. A plateau vs threshold means real pulses; a rate that keeps rising as the threshold drops means
noise/ringing. Writes figures/threshold_rate.json/.png.
"""
import argparse, sys, json, pathlib
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht

p = argparse.ArgumentParser()
p.add_argument("--thr", type=float, nargs="*", default=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.3])
p.add_argument("--edges", nargs="*", default=["Rising Edge"])
a = p.parse_args()
exp = ht.new_experiment("agent_thr_rate")
base = dict(gain=9000, trigger_source="External", roi_center=(500, 455), gate_delay_ns=27.0, gate_width_ns=3.0,
            kind="signal", notes="agent session 6 threshold vs rate", keep_roi_frames=0)
out = []
cam, lim = ht.open_camera()
try:
    ht.cool_down(cam, setpoint_c=-20.0, poll_s=30)
    for edge in a.edges:
        for v in a.thr:
            s = {**base, "trigger_threshold_v": v, "trigger_edge": edge}
            r = ht.take_run(cam, save=False, **{**s, "label": f"thr{v}", "accumulations": 20000, "n_readouts": 3, "max_seconds": 8})
            if r is None or ht.n_rd(r) < 2:
                out.append(dict(edge=edge, thr=v, rate_hz=None, note="no readouts")); continue
            dt = ht.readout_interval_s(r); rate = 20000 / dt; limited = dt < 0.3
            if limited:
                acc = 100000
                r = ht.take_run(cam, save=False, **{**s, "label": f"thr{v}b", "accumulations": acc, "n_readouts": 3, "max_seconds": 8})
                dt = ht.readout_interval_s(r); rate = acc / dt; limited = dt < 0.3
            out.append(dict(edge=edge, thr=v, rate_hz=rate, dt=dt, lower_bound=limited))
            print(f"RATE {edge} {v:.2f} V: {rate / 1e3:.1f} kHz{' (lower bound)' if limited else ''} dt {dt:.3f}", flush=True)
finally:
    ht.shutdown(cam)
(exp / "figures" / "threshold_rate.json").write_text(json.dumps(ht._jsonable(out), indent=1))
fig, ax = plt.subplots(figsize=(6, 4))
for edge in a.edges:
    o = [x for x in out if x["edge"] == edge and x["rate_hz"]]
    ax.plot([x["thr"] for x in o], [x["rate_hz"] / 1e3 for x in o], "o-", label=edge)
ax.set_xlabel("trigger threshold (V, 50 ohm)"); ax.set_ylabel("camera trigger rate (k/s)"); ax.set_yscale("log")
ax.grid(alpha=.3); ax.legend(); ax.set_title("External trigger rate vs threshold")
fig.tight_layout(); fig.savefig(exp / "figures" / "threshold_rate.png", dpi=110)
print("EXP", exp)
