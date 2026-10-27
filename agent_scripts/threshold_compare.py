"""Lowest thresholds and plain ADU summing (no threshold) on the 1533 focused heralded run."""
import sys, numpy as np
sys.path.insert(0, r"C:/Users/ReschLab/Single-Photon-Imaging")
import hoi_troubleshoot as ht
rs = ht.load_runs(r"C:/Users/ReschLab/Desktop/HOI_freshstaret/data/troubleshoot/2026-10-26/1533_agent_target31_newspot")
yy, xx = np.mgrid[:150, :150]
for mode in ("thr200", "thr30", "adu"):
    out = {}
    for r in rs:
        roi = r["meta"]["roi"]; near = np.hypot(xx - (462 - roi[0]), yy - (460 - roi[1])) <= 10
        ring = (np.hypot(xx - (462 - roi[0]), yy - (460 - roi[1])) > 40)
        st = ht.above_pedestal(r).astype(float)
        v = st if mode == "adu" else (st >= float(mode[3:])).astype(float)
        sig = v[:, near].sum(1) - v[:, ring].mean(1) * near.sum()       # local ring subtraction per readout
        out.setdefault(r["meta"]["settings"]["kind"], []).append(sig)
    a, b = np.concatenate(out["signal"]), np.concatenate(out["dark_offdelay"])
    net = a.mean() - b.mean(); err = np.hypot(a.std(ddof=1) / np.sqrt(len(a)), b.std(ddof=1) / np.sqrt(len(b)))
    print(f"{mode:6s}: net {net:+10.2f} +- {err:8.2f}  z {net/err:+.1f}")
