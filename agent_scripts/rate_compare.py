"""Rate test summary: per arm mean net spot score (vs OFF), ratio 4 ns / 31 ns, block chi2, gates/s.
    python agent_scripts/rate_compare.py EXP [EXP ...]          # each experiment on its own
    python agent_scripts/rate_compare.py --merge EXP [EXP ...]  # all runs pooled (one cycle per experiment, rate re-locked between)"""
import sys, numpy as np
sys.path.insert(0, r"C:/Users/ReschLab/Single-Photon-Imaging"); sys.path.insert(0, r"C:/Users/ReschLab/Single-Photon-Imaging/agent_scripts")
import hoi_troubleshoot as ht
from analyze_ab import beam_template
from scipy.stats import chi2
TEMPLATE = r"C:/Users/ReschLab/Desktop/HOI_freshstaret/data/troubleshoot/2026-10-26/1525_agent_beam_centered_min"
args = sys.argv[1:]
merge = args[:1] == ["--merge"]
exps = args[1:] if merge else args
groups = [(f"{len(exps)} merged", exps)] if merge else [(e, [e]) for e in exps]
for exp, folders in groups:
    rs = [r for f in folders for r in ht.load_runs(f)]
    w, _ = beam_template(rs[0]["meta"]["roi"], (462, 460), template_exp=TEMPLATE)
    eff = (w * w / w.sum()).sum()
    ctrl = [np.roll(np.roll(w, dy, 0), dx, 1) for dx, dy in ((45, 45), (-45, 45), (45, -45), (-45, -45))]
    arms = {}
    for r in rs:
        hit = ht.above_pedestal(r) >= 200
        M = ((hit * w).sum((1, 2)) - np.mean([(hit * c).sum((1, 2)) for c in ctrl], 0)) / eff
        lab = r["meta"]["settings"]["label"]; arm = "off" if "off" in lab else ("w4" if "w4" in lab else "w31")
        acc = r["meta"]["actual"]["accumulations"]
        arms.setdefault(arm, []).append((M, acc / ht.readout_interval_s(r)))
    off = np.concatenate([m for m, _ in arms["off"]])
    out = {}
    for k in ("w4", "w31"):
        allM = np.concatenate([m for m, _ in arms[k]])
        net = allM.mean() - off.mean(); err = np.hypot(allM.std(ddof=1) / np.sqrt(len(allM)), off.std(ddof=1) / np.sqrt(len(off)))
        bm = np.array([m.mean() for m, _ in arms[k]]); be = np.array([m.std(ddof=1) / np.sqrt(len(m)) for m, _ in arms[k]])
        c2 = (((bm - bm.mean()) / be) ** 2).sum(); dof = len(bm) - 1
        out[k] = (net, err)
        print(f"{exp.split('_agent_')[-1]:14s} {k:4s} net {net:+.3f} +- {err:.3f} (z {net/err:+.1f})  blocks {len(bm)}  chi2/dof {c2:.1f}/{dof} p {chi2.sf(c2, dof):.3f}")
    g = np.median([gps for v in arms.values() for _, gps in v])
    r4, r31 = out["w4"], out["w31"]; ratio = r4[0] / r31[0]
    rerr = abs(ratio) * np.hypot(r4[1] / r4[0], r31[1] / r31[0])
    print(f"{'':14s} ratio 4ns/31ns {ratio:.2f} +- {rerr:.2f} | median gates/s {g:,.0f}\n")
