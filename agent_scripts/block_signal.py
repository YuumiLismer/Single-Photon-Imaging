"""Per-block (50 readouts) heralded spot score for ON runs, OFF mean as reference."""
import sys, numpy as np
sys.path.insert(0, r"C:/Users/ReschLab/Single-Photon-Imaging"); sys.path.insert(0, r"C:/Users/ReschLab/Single-Photon-Imaging/agent_scripts")
import hoi_troubleshoot as ht
from analyze_ab import beam_template
exp, TEMPLATE = sys.argv[1], sys.argv[2]
rs = ht.load_runs(exp)
w, _ = beam_template(rs[0]["meta"]["roi"], (462, 460), template_exp=TEMPLATE)
eff = (w * w / w.sum()).sum()
ctrl = [np.roll(np.roll(w, dy, 0), dx, 1) for dx, dy in ((45, 45), (-45, 45), (45, -45), (-45, -45))]
out = []
for r in rs:
    hit = ht.above_pedestal(r) >= 200
    M = ((hit * w).sum((1, 2)) - np.mean([(hit * c).sum((1, 2)) for c in ctrl], 0)) / eff
    s = r["meta"]["settings"]; a = r["meta"]["actual"]
    out.append((r["meta"]["created"] if "created" in r["meta"] else s.get("label"), s.get("label"), a["gate_delay_ns"], a["gate_width_ns"], M.mean(), M.std(ddof=1) / np.sqrt(len(M)), ht.readout_interval_s(r)))
offm = np.mean([o[4] for o in out if "off" in o[1]])
print(f"OFF mean spot score {offm:+.2f}")
for o in sorted(out, key=lambda o: str(o[0])):
    if "off" in o[1]: continue
    print(f"{o[1]:28s} {o[2]:.0f}+{o[3]:.0f}  net {o[4]-offm:+.2f} +- {o[5]:.2f}  dt {o[6]:.2f}")
