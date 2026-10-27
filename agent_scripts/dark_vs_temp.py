"""Does the time-dependent single-pixel background drop with sensor temperature? Alternates set points
(default -20, -25, -20, -25 C), takes one OFF-gate block at each (herald-triggered, gate far from the photons) and
reports hits per second of readout wall time over the ROI, split into single pixels and clusters.

    python agent_scripts/dark_vs_temp.py [--temps -20 -25 -20 -25] [--readouts 30]
"""
import argparse
import pathlib
import sys

import numpy as np
from scipy.ndimage import label

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht

p = argparse.ArgumentParser()
p.add_argument("--temps", type=float, nargs="+", default=[-20, -25, -20, -25])
p.add_argument("--readouts", type=int, default=30)
p.add_argument("--acc", type=int, default=20000)
p.add_argument("--spot", type=int, nargs=2, default=(298, 277))
p.add_argument("--thr", type=float, default=200)
a = p.parse_args()

exp = ht.new_experiment("agent_dark_vs_temp")
base = dict(gain=9000, trigger_source="External", trigger_threshold_v=0.5, roi_center=tuple(a.spot), roi_size=400,
            kind="dark_offdelay", accumulations=a.acc, n_readouts=a.readouts, gate_delay_ns=1546.0, gate_width_ns=31.0)
rows = []
cam, _ = ht.open_camera()
try:
    for i, sp in enumerate(a.temps):
        locked = ht.cool_down(cam, setpoint_c=sp, timeout_min=6, poll_s=15)
        temp = float(cam.get_attribute_value("Sensor Temperature Reading"))
        r = ht.take_run(cam, **{**base, "label": f"dark_{i}_{sp:g}C"})
        st = ht.above_pedestal(r)
        dt = ht.readout_interval_s(r)
        n1 = n2 = 0
        for fr in st:
            lab, n = label(fr >= a.thr)
            sizes = np.bincount(lab.ravel())[1:]
            n1 += (sizes == 1).sum(); n2 += (sizes >= 2).sum()
        k = len(st)
        rows.append((sp, temp, locked, dt, n1 / k / dt, n2 / k / dt))
        print(f"set {sp:g} C (read {temp:g} C, {'locked' if locked else 'NOT locked'}): {dt:.2f} s/readout | "
              f"single-px {n1 / k / dt:7.1f}/s  clusters {n2 / k / dt:5.1f}/s", flush=True)
    if a.temps[-1] != -20:
        ht.cool_down(cam, setpoint_c=-20.0, timeout_min=6, poll_s=15)
finally:
    ht.shutdown(cam)
for sp in sorted(set(a.temps)):
    v = [r[4] for r in rows if r[0] == sp]
    print(f"mean single-px rate at {sp:g} C: {np.mean(v):.1f}/s ({len(v)} blocks)")
print("EXP", exp)
