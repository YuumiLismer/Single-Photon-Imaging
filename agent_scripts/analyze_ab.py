"""Offline re-analysis of an interleaved on/off experiment (multi_ab.py / ab_test.py) (agent, session 2).

    python agent_scripts/analyze_ab.py <experiment folder> --on on --off off [--thr 200]

Improvements over the online score (see agent_log session 2):
 * per-pixel baseline: each readout minus its median (as before) minus the per-pixel temporal median
   over the whole experiment. The beam spot sits ~12 ADU below its surroundings in long runs
   (a local 'hole' in the ungated background), which biased fixed-threshold hit counts;
 * regression instead of a fixed far-field subtraction: y = mu_arm + k * far_hits (fit on all arms),
   so readout-interval (ungated background) and in-gate stray light differences are absorbed by k;
 * three spot statistics: hits in r<=6, r<=10, and a matched filter weighted by the 12:21
   high-pump beam profile (signal minus signal-blocked, independent data; peak-normalised, r<=15).
Prints mu_on - mu_off (per readout and per trigger) with its z for each statistic.
"""
import argparse, sys, pathlib, json, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht

TEMPLATE_RUNS = ("2026-10-20/1221_notebook_session/runs/20261020T122130_A_wide_gate_signal.npz",
                 "2026-10-20/1221_notebook_session/runs/20261020T122510_A_wide_gate_blocked_dark.npz")


def beam_template(roi, spot=(501, 454), rmax=15):
    s, b = (ht.load_run(ht.DATA_ROOT / p) for p in TEMPLATE_RUNS)
    assert s["meta"]["roi"] == list(roi), "template ROI differs"
    d = s["roi_stack"].astype(float).mean(0) - b["roi_stack"].astype(float).mean(0)
    d -= np.median(d)
    x0, y0, n = roi
    yy, xx = np.indices((n, n)); rr = np.hypot(xx - (spot[0] - x0), yy - (spot[1] - y0))
    w = np.clip(d, 0, None) * (rr <= rmax)
    return w / w.max(), rr


def arm_of(run):
    m = re.search(r"_(on|off|short|[a-z0-9]+)_[\d.]+ns_c\d+", run["meta"]["settings"]["label"])
    return m.group(1) if m else run["meta"]["settings"]["label"]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("folder"); p.add_argument("--on", default="on"); p.add_argument("--off", default="off")
    p.add_argument("--thr", type=float, default=200); p.add_argument("--spot", type=int, nargs=2, default=(501, 454))
    a = p.parse_args()
    runs = ht.load_runs(a.folder)
    runs = [r for r in runs if "_c" in r["meta"]["settings"]["label"] and r["roi_stack"].size]
    roi = runs[0]["meta"]["roi"]
    w, rr = beam_template(roi, tuple(a.spot))
    prof = w / w.sum()            # beam profile (fraction of the spot's hits per pixel), r <= 15
    eff = dict(h6=prof[rr <= 6].sum(), h10=prof[rr <= 10].sum(), mf=(w * prof).sum())
    print("fraction of a beam-shaped spot counted by each statistic:", {k: round(float(v), 3) for k, v in eff.items()})
    stacks = [ht.above_pedestal(r) for r in runs]
    base = np.median(np.concatenate(stacks), axis=0)
    far = rr > 50
    rows = []
    for r, st in zip(runs, stacks):
        hit = (st - base) >= a.thr
        acc = int(r["meta"]["actual"]["accumulations"])
        for h in hit:
            rows.append(dict(arm=arm_of(r), far=h[far].mean() * far.sum(), h6=h[rr <= 6].sum(), h10=h[rr <= 10].sum(),
                             mf=(h * w).sum(), acc=acc, dt=ht.readout_interval_s(r)))
    arms = sorted({x["arm"] for x in rows})
    print("arms:", {k: sum(x["arm"] == k for x in rows) for k in arms})
    X = np.array([[x["arm"] == k for k in arms] + [x["far"]] for x in rows], float)
    out = {}
    for stat in ("h6", "h10", "mf"):
        y = np.array([x[stat] for x in rows], float)
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        res = y - X @ beta
        cov = np.linalg.inv(X.T @ X) * res.var(ddof=X.shape[1])
        i, j = arms.index(a.on), arms.index(a.off)
        d = beta[i] - beta[j]; e = np.sqrt(cov[i, i] + cov[j, j] - 2 * cov[i, j])
        # simple version too: per-arm mean of y - k*far with k from the fit
        out[stat] = dict(diff=d, err=e, z=d / e, k=beta[-1], mu={k: beta[n] for n, k in enumerate(arms)},
                         eff=eff[stat], spot_hits_per_trigger=d / eff[stat] / rows[0]["acc"],
                         err_per_trigger=e / eff[stat] / rows[0]["acc"],
                         per_trigger_2sigma_limit=(max(d, 0) + 2 * e) / eff[stat] / rows[0]["acc"])
        print(f"{stat:4s}: on-off = {d:+.3f} +- {e:.3f} per readout (z {d / e:+.2f}); k_far = {beta[-1]:.4f}; "
              f"means " + ", ".join(f"{k} {beta[n]:.2f}" for n, k in enumerate(arms)) +
              f"; spot hits/trigger {out[stat]['spot_hits_per_trigger']:+.2e} +- {out[stat]['err_per_trigger']:.2e}"
              f" (2-sigma limit {out[stat]['per_trigger_2sigma_limit']:.2e})")
    # beam level today (short-gate arm vs off arm), local-annulus estimate per s of gate-open
    print("dt per arm:", {k: round(np.mean([x["dt"] for x in rows if x["arm"] == k]), 3) for k in arms})
    fig_dir = pathlib.Path(a.folder) / "figures"
    (fig_dir / "analysis.json").write_text(json.dumps(ht._jsonable(out), indent=1))
    return out


if __name__ == "__main__":
    main()
