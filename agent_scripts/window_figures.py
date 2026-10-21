"""Per-window ON / OFF / ON-minus-OFF figures and a summary plot for an interleaved delay scan
(multi_ab.py with arms named on<delay>/off<delay>) (agent, session 4).

    python agent_scripts/window_figures.py <experiment folder> [--thr 200]

For each window: hit images (per-run per-pixel baseline as in analyze_ab.py) of the ON gate
(delay d) and the OFF gate (same width, d + 5 us), and ON minus OFF, with the spot (501, 454) marked.
Score per window = regression of the beam-profile matched filter (and r<=10 hits) on arm + readout
interval over ALL runs of the experiment (shared hits-per-second term), ON minus OFF of that window.
Net spot hits per readout = matched-filter difference / fraction of a beam-shaped spot it counts.
Writes figures/window_<d>ns.png, figures/summary_vs_delay.png and figures/windows.json.
"""
import argparse, sys, pathlib, json, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import hoi_troubleshoot as ht
from analyze_ab import beam_template, arm_of

p = argparse.ArgumentParser()
p.add_argument("folder"); p.add_argument("--thr", type=float, default=200)
p.add_argument("--spot", type=int, nargs=2, default=(501, 454))
p.add_argument("--max-dt", type=float, default=None, help="drop runs with readout interval above this (s)")
p.add_argument("--suffix", default="", help="appended to figure/json names")
p.add_argument("--method", choices=("control", "regress"), default="control",
               help="control (default since session 5): per-readout matched filter minus the same template at 4 control "
                    "positions (+-45 px), ON minus OFF with Welch errors. regress: old arm + readout-interval regression, "
                    "needed only when readouts were long (>1 s) and the ungated background depended on wall time")
p.add_argument("--skip", nargs="*", default=[], help="arms to leave out of the fit, e.g. on27 off27 (session 5: the 27 ns gate "
               "did not open, so its short readouts distorted the shared readout-interval term)")
p.add_argument("--shared-off", default=None, help="OFF arm used for every ON arm that has no off<d> of its own "
               "(session 6 gate-start scans with one off-delay reference)")
p.add_argument("--baseline", choices=("readout", "pixel"), default="readout",
               help="readout (default since session 6): hit = ADU above each readout's own pedestal. pixel: also subtract "
                    "each pixel's median over the run's readouts (old default); with a compact spot and strong heralded "
                    "signal (session 6: ~12 hits/readout in r<=10) the spot pixels' medians rise to ~200 ADU and about "
                    "half the photons are lost")
p.add_argument("--rate-band", type=float, default=None,
               help="readout post-selection: keep only readouts whose interval (accumulations / interval = trigger rate) "
                    "is within this fraction of the arm's median interval, e.g. 0.15 drops readouts taken while the "
                    "herald rate spiked or dropped by more than ~15%%. Not useful when readouts are readout-time limited "
                    "(0.226 s), where the interval no longer tracks the rate")
a = p.parse_args()
folder = pathlib.Path(a.folder); figs = folder / "figures"
runs = [r for r in ht.load_runs(str(folder)) if "_c" in r["meta"]["settings"]["label"] and r["roi_stack"].size]
if a.max_dt:
    runs = [r for r in runs if ht.readout_interval_s(r) <= a.max_dt]
    print(len(runs), "runs kept with readout interval <=", a.max_dt, "s")
roi = runs[0]["meta"]["roi"]; x0, y0, n = roi
w, rr = beam_template(roi, tuple(a.spot)); prof = w / w.sum(); eff_mf = (w * prof).sum()
ctrl = [np.roll(np.roll(w, dy, 0), dx, 1) for dx, dy in ((45, 45), (-45, 45), (45, -45), (-45, -45))]
def readout_dts(r):
    t = np.asarray(r["t"], float); dt = np.diff(t)
    return np.concatenate([[np.median(dt) if len(dt) else np.nan], dt])


med_dt = {}
if a.rate_band:
    for r in runs:
        med_dt.setdefault(arm_of(r), []).extend(readout_dts(r)[1:])
    med_dt = {k: float(np.median(v)) for k, v in med_dt.items()}
kept = {}
rows, img, per = [], {}, {}
for r in runs:
    st = ht.above_pedestal(r)
    hit = ((st - np.median(st, 0)) if a.baseline == "pixel" else st) >= a.thr
    arm = arm_of(r); s = r["meta"]["settings"]
    if a.rate_band and arm in med_dt:
        keep = np.abs(readout_dts(r)[:len(hit)] / med_dt[arm] - 1) <= a.rate_band
        k0, k1 = kept.get(arm, (0, 0)); kept[arm] = (k0 + int(keep.sum()), k1 + len(keep))
        hit = hit[keep]
        if len(hit) == 0:
            continue
    if arm in a.skip:
        continue
    rows.append(dict(arm=arm, n=len(hit), dt=ht.readout_interval_s(r), h10=hit[:, rr <= 10].sum(1).mean(),
                     mf=(hit * w).sum((1, 2)).mean(), acc=int(r["meta"]["actual"]["accumulations"]),
                     delay=float(s["gate_delay_ns"]), width=float(s["gate_width_ns"]), gain=int(s["gain"])))
    sm, k = img.get(arm, (0, 0)); img[arm] = (sm + hit.sum(0), k + len(hit))
    c = np.mean([(hit * cw).sum((1, 2)) for cw in ctrl], 0)
    per.setdefault(arm, {"mf": [], "h10": []})
    per[arm]["mf"].append((hit * w).sum((1, 2)) - c); per[arm]["h10"].append(hit[:, rr <= 10].sum(1))
arms = sorted({x["arm"] for x in rows})
X = np.array([[x["arm"] == k for k in arms] + [x["dt"]] for x in rows], float); W = np.sqrt([x["n"] for x in rows])


def fit(stat):
    y = np.array([x[stat] for x in rows], float)
    beta, *_ = np.linalg.lstsq(X * W[:, None], y * W, rcond=None)
    res = (y - X @ beta) * W
    cov = np.linalg.inv((X * W[:, None]).T @ (X * W[:, None])) * res.var(ddof=X.shape[1])
    return beta, cov


fits = {s: fit(s) if a.method == "regress" else (None, None) for s in ("mf", "h10")}
if kept:
    print("rate-band post-selection, readouts kept per arm:", {k: f"{v[0]}/{v[1]}" for k, v in kept.items()})
wins = sorted((int(re.sub(r"\D", "", k)), k) for k in arms if k.startswith("on") and re.sub(r"\D", "", k))
out = []
for d, on in wins:
    off = f"off{d}"
    tag = f"{d}ns" if on == f"on{d}" else on      # figure name; e.g. 'onc577' (session 6 extra arms vs a shared OFF)
    if off not in arms and a.shared_off:
        off = a.shared_off
    if on not in arms or off not in arms:
        continue
    i, j = arms.index(on), arms.index(off); res = {}
    for s, (beta, cov) in fits.items():
        if a.method == "regress":
            diff = beta[i] - beta[j]; err = np.sqrt(cov[i, i] + cov[j, j] - 2 * cov[i, j])
        else:
            x, y = np.concatenate(per[on][s]), np.concatenate(per[off][s])
            diff = x.mean() - y.mean(); err = np.hypot(x.std(ddof=1) / np.sqrt(len(x)), y.std(ddof=1) / np.sqrt(len(y)))
        res[s] = (diff, err)
    R = [x for x in rows if x["arm"] == on]; Roff = [x for x in rows if x["arm"] == off]
    acc = R[0]["acc"]; nro = sum(x["n"] for x in R); width = R[0]["width"]; gain = R[0]["gain"]
    net, err = res["mf"][0] / eff_mf, res["mf"][1] / eff_mf
    rec = dict(arm=on, delay=R[0]["delay"], width=width, gain=gain, acc=acc, readouts_on=nro, readouts_off=sum(x["n"] for x in Roff),
               off_delay=Roff[0]["delay"], net_spot_per_readout=net, err=err, z=net / err,
               r10_diff=res["h10"][0], r10_err=res["h10"][1], z_r10=res["h10"][0] / res["h10"][1],
               per_trigger=net / acc, limit_2sigma_per_trigger=(max(net, 0) + 2 * err) / acc,
               s_per_readout_on=float(np.mean([x["dt"] for x in R])))
    out.append(rec)
    ion, ioff = img[on][0] / img[on][1], img[off][0] / img[off][1]
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    vmax = np.percentile(gaussian_filter(ion, 1), 99.5)
    for a_, im, t in ((ax[0], ion, f"ON: gate {d}-{d + width:g} ns"), (ax[1], ioff, f"OFF: gate {rec['off_delay']:g}-{rec['off_delay'] + width:g} ns")):
        h = a_.imshow(gaussian_filter(im, 1), cmap="inferno", vmin=0, vmax=vmax, extent=(x0 - .5, x0 + n - .5, y0 + n - .5, y0 - .5))
        a_.set_title(t + "\nphoton hits per pixel per readout"); fig.colorbar(h, ax=a_, fraction=0.046)
    dimg = gaussian_filter(ion - ioff, 3); lim = np.abs(dimg).max()
    h = ax[2].imshow(dimg, cmap="RdBu_r", vmin=-lim, vmax=lim, extent=(x0 - .5, x0 + n - .5, y0 + n - .5, y0 - .5))
    ax[2].set_title("ON minus OFF = net heralded image\n(3 px smoothing; red = more in ON)"); fig.colorbar(h, ax=ax[2], fraction=0.046)
    for a_ in ax:
        a_.add_patch(plt.Circle(tuple(a.spot), 10, fill=False, color="cyan", lw=1.2))
        a_.plot(*a.spot, "c+", ms=8); a_.set_xlabel("x (px)"); a_.set_ylabel("y (px)")
    fig.suptitle(f"Gate {rec['delay']:g} ns, width {width:g} ns, gain {gain}, {acc:,} accumulations/readout, "
                 f"{nro} ON + {rec['readouts_off']} OFF readouts, External trigger\n"
                 f"net spot hits/readout {net:+.2f} +- {err:.2f}  ->  z = {net / err:+.1f}   "
                 f"(cyan circle = spot (501, 454), r = 10 px; >= 4 = candidate, >= 5 = detection)")
    fig.tight_layout(); fig.savefig(figs / f"window_{tag}{a.suffix}.png", dpi=100); plt.close(fig)
    print(f"{on:>7s} {R[0]['delay']:g}+{width:g} ns: net {net:+.3f} +- {err:.3f} spot hits/readout, z {net / err:+.2f} | r<=10 {res['h10'][0]:+.3f} +- {res['h10'][1]:.3f}"
          f" | {nro} readouts, {rec['s_per_readout_on']:.2f} s/readout | 2-sigma limit {rec['limit_2sigma_per_trigger']:.1e}/trigger")
if out:
    fig, ax = plt.subplots(figsize=(8, 4.5))
    D = np.array([o["delay"] + o["width"] / 2 for o in out])
    ax.errorbar(D, [o["net_spot_per_readout"] for o in out], yerr=[o["err"] for o in out], xerr=[o["width"] / 2 for o in out], fmt="o", capsize=4)
    ax.axhline(0, color="k", lw=0.8)
    for o, x in zip(out, D):
        ax.annotate(f"z {o['z']:+.1f}", (x, o["net_spot_per_readout"] + o["err"]), ha="center", va="bottom", fontsize=9)
    ax.set_xlabel("gate window after the herald trigger (ns)"); ax.set_ylabel("net spot hits per readout (ON - OFF)")
    ax.set_title(f"Heralded-photon search, {out[0]['width']:g} ns gates, {out[0]['acc']:,} triggers/readout, gain {out[0]['gain']}")
    fig.tight_layout(); fig.savefig(figs / f"summary_vs_delay{a.suffix}.png", dpi=110); plt.close(fig)
(figs / f"windows{a.suffix}.json").write_text(json.dumps(ht._jsonable(out), indent=1))
