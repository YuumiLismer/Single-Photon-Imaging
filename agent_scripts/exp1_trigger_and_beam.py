"""Agent experiment 1 (20 Oct, remote): trigger-rate vs threshold/edge, then the wide-gate beam check.

Why: the trigger rate is 45 kHz idler-open / 5 kHz blocked. A plateau in rate vs threshold means
real pulses; a rate that climbs as the threshold drops means noise. Then confirm the beam spot
at (501, 454) with more readouts than the last hunt (which used 20 and found contrast ~1.4 there).
"""
import sys, json, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")   # plt.show() must not block an unattended script
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht

SPOT = (501, 454)
ht.new_experiment("agent_exp1_wide_gate")      # runs and figures go in one dated folder
out = {}
cam, limits = ht.open_camera()
try:
    ht.cool_down(cam, setpoint_c=-20.0, poll_s=30)
    base = dict(gain=1, gate_delay_ns=309, gate_width_ns=3, accumulations=100_000, n_readouts=4,
                trigger_source="External", trigger_edge="Rising Edge", roi_center=(500, 455),
                frame_timeout_s=12, notes="agent exp1 trigger sweep")
    rates = []
    for edge, thr in [("Rising Edge", v) for v in (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0)] + [("Falling Edge", 0.5)]:
        r = ht.take_run(cam, save=False, **{**base, "label": f"trig_{edge[:4]}_{thr}", "trigger_edge": edge,
                                             "trigger_threshold_v": thr})
        if r is None or ht.n_rd(r) < 3:
            rate = 0.0
        else:
            rate = int(r["meta"]["actual"]["accumulations"]) / ht.readout_interval_s(r)
        rates.append((edge, thr, rate))
        print(f"RATE {edge} {thr} V: {rate / 1e3:.2f} kHz", flush=True)
    out["rates"] = rates

    # wide-gate beam checks (gain 9000; gate x acc = 0.5 s <= 1 s)
    wide = dict(gain=9000, gate_delay_ns=309, trigger_source="External", trigger_threshold_v=0.5,
                trigger_edge="Rising Edge", roi_center=(500, 455), notes="agent exp1 wide gate")
    for lab, w, acc, n in (("agent_wide100us", 100_000, 5000, 40), ("agent_wide20us", 20_000, 25_000, 40)):
        r = ht.take_run(cam, **{**wide, "label": lab, "gate_width_ns": w, "accumulations": acc, "n_readouts": n})
        sat = float((r["first"] >= ht.SATURATION_ADU).mean())
        s = ht.spot_series(r, SPOT)
        ap = ht.spot_masks(r, SPOT)[0]
        h = ht.hit_series(r, 200, ap)
        found = ht.find_spot(r)
        res = dict(sat=sat, spot_adu=s.mean(), spot_err=s.std(ddof=1) / np.sqrt(len(s)),
                   z=s.mean() / (s.std(ddof=1) / np.sqrt(len(s))), hits=h.mean(),
                   contrast=ht.spot_contrast(r, SPOT), found=found, contrast_found=ht.spot_contrast(r, found),
                   dt=ht.readout_interval_s(r))
        out[lab] = res
        print("WIDE", lab, json.dumps(ht._jsonable(res)), flush=True)
        ht.plot_run(r, None, threshold=200, spot=SPOT, save_as=f"agent_exp1_{lab}")
        if sat > 1e-3:
            print("SATURATION > 0.1%: stopping"); break
finally:
    ht.shutdown(cam)
print("RESULT", json.dumps(ht._jsonable(out)))
