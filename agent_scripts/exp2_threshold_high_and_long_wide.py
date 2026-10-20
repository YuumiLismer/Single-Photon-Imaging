"""Agent experiment 2: (a) trigger rate at thresholds above 1 V, to find the herald pulse height
(exp1 showed a plateau of 38-43 kHz from 0.2 to 1.0 V); (b) a long wide-gate run at the light-safety
limit (20 us gate x 50,000 acc = 1 s gate-open per readout) to measure how much signal-arm light
reaches (501, 454) now (exp1's 100 us x 5000 x 40 run showed nothing there: -1.3 sigma)."""
import sys, json, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht

SPOT = (501, 454)
N_WIDE = int(sys.argv[1]) if len(sys.argv) > 1 else 100
cam, limits = ht.open_camera()
try:
    ht.cool_down(cam, setpoint_c=-20.0, poll_s=30)
    base = dict(gain=1, gate_delay_ns=309, gate_width_ns=3, accumulations=100_000, n_readouts=4,
                trigger_source="External", trigger_edge="Rising Edge", roi_center=(500, 455),
                frame_timeout_s=12, notes="agent exp2 trigger sweep")
    for thr in (1.5, 2.0, 2.5, 3.0, 0.5):
        r = ht.take_run(cam, save=False, **{**base, "label": f"trig_{thr}", "trigger_threshold_v": thr})
        rate = 0.0 if (r is None or ht.n_rd(r) < 3) else int(r["meta"]["actual"]["accumulations"]) / ht.readout_interval_s(r)
        print(f"RATE Rising {thr} V (read back {None if r is None else r['meta']['actual']['trigger_threshold_v']}): {rate / 1e3:.2f} kHz", flush=True)
    r = ht.take_run(cam, label="agent_wide20us_1s", gain=9000, gate_delay_ns=309, gate_width_ns=20_000,
                    accumulations=50_000, n_readouts=N_WIDE, trigger_source="External", trigger_threshold_v=0.5,
                    roi_center=(500, 455), notes="agent exp2 long wide gate")
    sat = float((r["first"] >= ht.SATURATION_ADU).mean())
    s = ht.spot_series(r, SPOT)
    h = ht.hit_series(r, 200, ht.spot_masks(r, SPOT)[0])
    res = dict(sat=sat, spot_adu=s.mean(), spot_err=s.std(ddof=1) / np.sqrt(len(s)), hits=h.mean(),
               hits_err=h.std(ddof=1) / np.sqrt(len(h)), contrast=ht.spot_contrast(r, SPOT),
               found=ht.find_spot(r), dt=ht.readout_interval_s(r))
    print("WIDE", json.dumps(ht._jsonable(res)), flush=True)
    ht.plot_run(r, None, threshold=200, spot=SPOT, save_as="agent_exp2_wide20us_1s")
finally:
    ht.shutdown(cam)
print("DONE")
