"""Arrival-time profile of the heralded photons: net spot hits per trigger vs gate delay, from the
windows.json of one or more window_figures.py runs (agent session 5; Tristan: "see the decay").

    python agent_scripts/arrival_profile.py <out folder> <fit exp folder> [<other exp folders> ...]

Each window is drawn as a horizontal bar over its gate (delay .. delay + width) at its net heralded
spot hits per trigger, with 1 sigma error bars. The first experiment (the short-gate scan) is fitted
with a Gaussian, convolved with its gate (box of the gate width), plus a constant: centre and FWHM
are quoted. If every window of the first experiment has the same END (a gate-start scan, used in
session 5 because gates that close before ~560-577 ns after the trigger do not open), the signal is
the arrival distribution integrated from the gate start onwards, so a complementary error function
(Gaussian arrival) is fitted vs gate start instead; such sets are drawn at their gate start.
No camera time. Writes arrival_profile.png and arrival_profile.json.
"""
import sys, json, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import erf

out = pathlib.Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
sets = []
for f in sys.argv[2:]:
    W = json.loads((pathlib.Path(f) / "figures" / "windows.json").read_text())
    sets.append((pathlib.Path(f).name, W))


def model(d, amp, mu, sig, c, width):
    # Gaussian arrival (area amp per trigger) integrated over the gate [d, d + width]
    s = max(abs(sig), 0.3)
    return amp * 0.5 * (erf((d + width - mu) / (np.sqrt(2) * s)) - erf((d - mu) / (np.sqrt(2) * s))) + c


name, W = sets[0]
d = np.array([w["delay"] for w in W], float); y = np.array([w["per_trigger"] for w in W]); e = np.array([w["err"] / w["acc"] for w in W])
width = W[0]["width"]; fit = {}
ends = {round(w["delay"] + w["width"]) for w in W}
start_scan = len(ends) == 1 and len(W) > 2
if start_scan:
    end = ends.pop()
    def model(d, amp, mu, sig, c, width):   # noqa: F811  (gate [d, end]: photons arriving after d)
        s = max(abs(sig), 0.3)
        return amp * 0.5 * (erf((end - mu) / (np.sqrt(2) * s)) - erf((d - mu) / (np.sqrt(2) * s))) + c
try:
    i = int(np.argmax(y))
    p0 = (max(y.max(), 1e-6), d[i] + width / 2, width / 2, 0.0) if not start_scan else         (max(y.max(), 1e-6), float(d[np.argmin(np.abs(y - y.max() / 2))]), 2.0, 0.0)
    f = lambda x, amp, mu, sig, c: model(x, amp, mu, sig, c, width)
    popt, pcov = curve_fit(f, d, y, p0=p0, sigma=e, absolute_sigma=True, maxfev=20000)
    perr = np.sqrt(np.diag(pcov))
    fit = dict(amp_per_trigger=popt[0], amp_err=perr[0], centre_ns=popt[1], centre_err=perr[1],
               fwhm_ns=2.3548 * abs(popt[2]), fwhm_err=2.3548 * perr[2], const=popt[3], gate_ns=width,
               chi2=float((((y - f(d, *popt)) / e) ** 2).sum()), dof=len(d) - 4)
except Exception as ex:
    fit = dict(error=str(ex))
fig, ax = plt.subplots(figsize=(10, 5.2))
for k, (nm, Wk) in enumerate(sets):
    for j, w in enumerate(Wk):
        pt, er = w["per_trigger"], w["err"] / w["acc"]
        ss = len({round(v["delay"] + v["width"]) for v in Wk}) == 1 and len(Wk) > 2
        ax.errorbar(w["delay"] if ss else w["delay"] + w["width"] / 2, pt, yerr=er, xerr=0 if ss else w["width"] / 2, fmt="o" if k == 0 else "s",
                    color=f"C{k}", ms=5 if k == 0 else 4, capsize=2, alpha=1 if k == 0 else 0.6,
                    label=(f"{nm}: gate [start, {w['delay'] + w['width']:g} ns], plotted at start, {w['acc']:,} acc" if ss else
                           f"{nm} ({w['width']:g} ns gates, {w['acc']:,} acc, {w['readouts_on']} ON readouts)") if j == 0 else None)
if "centre_ns" in fit:
    xx = np.linspace(d.min() - 20, d.max() + 20, 600)
    if np.isfinite(fit["fwhm_err"]):
        lab = (f"fit: centre {fit['centre_ns']:.1f} +- {fit['centre_err']:.1f} ns, "
               f"FWHM {fit['fwhm_ns']:.1f} +- {fit['fwhm_err']:.1f} ns (gate-deconvolved)")
    else:   # edge sharper than the scan step: quote the bracketing points instead
        lo = float(d[d <= fit["centre_ns"]].max()) if (d <= fit["centre_ns"]).any() else np.nan; hi = float(d[d >= fit["centre_ns"]].min()) if (d >= fit["centre_ns"]).any() else np.nan
        step = float(np.min(np.diff(np.unique(d))))
        fit.update(bracket_ns=[lo, hi], fwhm_upper_ns=float(hi - lo) if np.isfinite(hi - lo) else step)
        lab = (f"fit: arrival between {lo:g} and {hi:g} ns (centre ~{(lo + hi) / 2:.1f} ns); "
               f"width unresolved, FWHM < ~{fit['fwhm_upper_ns']:g} ns (step-limited)")
    ax.plot(xx if start_scan else xx + width / 2, model(xx, fit["amp_per_trigger"], fit["centre_ns"], fit["fwhm_ns"] / 2.3548, fit["const"], width),
            "C0-", lw=1.2, label=lab)
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("time after the herald trigger (ns)  [bar = gate window; start-scan points at the gate start]")
ax.set_ylabel("net heralded spot hits per trigger (ON - OFF)")
ax.set_title("Arrival-time profile of heralded photons at the camera (External trigger, gain 9000, spot (501, 454))")
ax.legend(fontsize=8); ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig(out / "arrival_profile.png", dpi=120)
(out / "arrival_profile.json").write_text(json.dumps(dict(fit=fit, sets=[n for n, _ in sets]), indent=1, default=float))
print(json.dumps(fit, default=float))
