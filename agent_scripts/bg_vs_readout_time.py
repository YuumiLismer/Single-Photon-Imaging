"""Per run: median readout interval vs background hits per readout (ROI minus spot)."""
import sys, numpy as np
sys.path.insert(0, r"C:/Users/ReschLab/Single-Photon-Imaging")
import hoi_troubleshoot as ht
yy, xx = np.mgrid[:150, :150]; far = np.hypot(xx - 75, yy - 75) > 15
rows = []
for exp in sys.argv[1:]:
    for r in ht.load_runs(exp):
        hit = ht.above_pedestal(r) >= 200
        bg = hit[:, far].sum(1).mean()
        dt = ht.readout_interval_s(r)
        m = r["meta"]; a = m["actual"]
        rows.append((m["settings"].get("label", "?"), a["gate_delay_ns"], a["gate_width_ns"], dt, bg))
rows.sort(key=lambda x: x[0])
for l, d, w, dt, bg in rows:
    print(f"{l:34s} {d:6.0f}+{w:<3.0f} {dt:5.2f} s/readout  bg {bg:6.1f}  bg per s {bg/dt:6.1f}")
dts = np.array([r[3] for r in rows]); bgs = np.array([r[4] for r in rows])
print("corr(bg, dt) =", np.corrcoef(dts, bgs)[0, 1])
