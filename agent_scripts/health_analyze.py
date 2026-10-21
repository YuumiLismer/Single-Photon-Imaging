"""Offline analysis of the internal-trigger camera health check (agent_scripts/health_check.py).

    python agent_scripts/health_analyze.py <experiment folders...> [--out <figures folder>]

Figures (in --out): timing.png, gain.png, gate.png, single_photons.png, full_frame.png, events.png;
numbers in health.json. Hits = pixels >= 200 ADU above the readout's median (ROI 150x150 at (500, 455)).
'Closed' reference = the same settings with a 3 ns gate (gate-open time negligible, same wall time).
"""
import argparse, sys, pathlib, json, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht

SPOT = (501, 454)
REF_1221 = "2026-10-20/1221_notebook_session/runs/20261020T122130_A_wide_gate_signal.npz"
REF_1225 = "2026-10-20/1221_notebook_session/runs/20261020T122510_A_wide_gate_blocked_dark.npz"


def A(r):
    return r["meta"]["actual"]


def lab(r):
    return r["meta"]["settings"]["label"]


def open_s(r):
    return float(A(r)["gate_width_ns"]) * 1e-9 * int(A(r)["accumulations"])


def above(r):
    st = r["roi_stack"].astype(float)
    return st - np.median(st.reshape(len(st), -1), 1)[:, None, None]


def radii(r, spot=SPOT):
    x0, y0, n = r["meta"]["roi"]
    yy, xx = np.indices((n, n))
    return np.hypot(xx - (spot[0] - x0), yy - (spot[1] - y0))


def core_minus_annulus(img, rr):
    """Sum over r<6 of (pixel - mean of the 10-16 px annulus): compact-spot signal (E2 definition)."""
    ann = img[(rr >= 10) & (rr <= 16)].mean()
    return float((img[rr < 6] - ann).sum())


def clusters(frame, thr=200, lo=60):
    """Events: connected regions >= lo ADU that contain a pixel >= thr. List of (y, x, peak, sum, npix)."""
    lb, n = ndi.label(frame >= lo)
    if n == 0:
        return []
    idx = np.arange(1, n + 1)
    peak = ndi.maximum(frame, lb, idx)
    tot = ndi.sum(frame, lb, idx)
    npx = ndi.sum(np.ones_like(frame), lb, idx)
    com = ndi.center_of_mass(frame, lb, idx)
    return [(c[0], c[1], p, t, k) for c, p, t, k in zip(com, peak, tot, npx) if p >= thr]


def timing(runs, out, res):
    T = [r for r in runs if re.match(r"[SA]_", lab(r))]
    if not T:
        return
    fig, ax = plt.subplots(figsize=(5.5, 4))
    rows = []
    for r in T:
        dt = ht.readout_interval_s(r)
        src = str(A(r)["trigger_source"])
        exp_dt = int(A(r)["accumulations"]) / float(A(r)["internal_trigger_hz"]) if src == "Internal" else None
        rows.append(dict(label=lab(r), source=src, freq=A(r)["internal_trigger_hz"], acc=A(r)["accumulations"],
                         dt=dt, expected=exp_dt))
        if src == "Internal":
            ax.plot(exp_dt, dt, "o", color="C0")
        else:
            ax.axhline(dt, color="C3", ls=":", lw=1)
    x = np.array([0.05, 2.5])
    ax.plot(x, x, "k--", lw=.8, label="interval = accumulations / internal rate")
    ax.plot([], [], ":", color="C3", label="External (APD) runs in between (20k acc)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("accumulations / internal trigger rate (s)"); ax.set_ylabel("measured readout interval (s)")
    ax.set_title("Internal trigger takes: timing follows the clock", fontsize=10)
    ax.legend(fontsize=8); ax.grid(alpha=.3, which="both")
    fig.tight_layout(); fig.savefig(out / "timing.png", dpi=130); plt.close(fig)
    res["timing"] = rows


def gain_series(runs, out, res):
    G1 = [r for r in runs if A(r)["gain"] == 1 and A(r)["gate_width_ns"] < 5 and len(r["roi_stack"]) > 2]
    if G1:
        st = G1[0]["roi_stack"].astype(float)
        d = (st[1:] - st[:-1]).ravel() / np.sqrt(2)
        res["pedestal_adu_gain1"] = float(np.median(st))
        res["read_noise_adu_gain1"] = float(1.4826 * np.median(np.abs(d - np.median(d))))
    B = [r for r in runs if lab(r).startswith("B_")]
    if not B:
        return
    rows = []
    for g in sorted({A(r)["gain"] for r in B}):
        o = [r for r in B if A(r)["gain"] == g and A(r)["gate_width_ns"] > 5][0]
        c = [r for r in B if A(r)["gain"] == g and A(r)["gate_width_ns"] < 5][0]
        so, sc = o["roi_stack"].astype(float), c["roi_stack"].astype(float)
        dc = (sc[1:] - sc[:-1]).ravel() / np.sqrt(2)
        rows.append(dict(gain=g, roi_mean_open=float(so.mean()), roi_mean_closed=float(sc.mean()),
                         light_adu_px=float(so.mean() - sc.mean()),
                         light_err=float(np.hypot(so.mean((1, 2)).std(ddof=1) / np.sqrt(len(so)),
                                                  sc.mean((1, 2)).std(ddof=1) / np.sqrt(len(sc)))),
                         noise_closed_adu=float(1.4826 * np.median(np.abs(dc - np.median(dc)))),
                         hits200_open=float((above(o) >= 200).sum((1, 2)).mean()),
                         hits200_closed=float((above(c) >= 200).sum((1, 2)).mean()),
                         sat_open=float((o["first"] >= ht.SATURATION_ADU).mean())))
    res["gain"] = rows
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    g = np.array([x["gain"] for x in rows]); L = np.array([x["light_adu_px"] for x in rows])
    E = np.array([x["light_err"] for x in rows])
    ok = L > 0
    ax[0].errorbar(g[ok], L[ok], E[ok], fmt="o-", label="ROI mean: 10 us gate minus 3 ns gate")
    ax[0].plot(g, L[-1] * g / g[-1], "k--", lw=.8, label="proportional to gain")
    ax[0].set_xscale("log"); ax[0].set_yscale("log")
    ax[0].set_xlabel("emICCD gain"); ax[0].set_ylabel("in-gate light (ADU / pixel / readout)")
    ax[0].set_title("Intensifier gain response (internal 10 kHz, 0.1 s open)", fontsize=10)
    ax[0].legend(fontsize=8); ax[0].grid(alpha=.3, which="both")
    ax[1].plot(g, [x["hits200_open"] for x in rows], "o-", label="10 us gate (0.1 s open)")
    ax[1].plot(g, [x["hits200_closed"] for x in rows], "s-", label="3 ns gate (closed)")
    ax[1].set_xscale("log"); ax[1].set_xlabel("emICCD gain")
    ax[1].set_ylabel("pixels >= 200 ADU above pedestal / readout (ROI)")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=.3); ax[1].set_title("Photon hits vs gain", fontsize=10)
    fig.tight_layout(); fig.savefig(out / "gain.png", dpi=130); plt.close(fig)


def gate_series(runs, out, res):
    C = [r for r in runs if lab(r).startswith("C_")]
    if not C:
        return
    rows = []
    for r in C:
        h = above(r) >= 200
        rr = radii(r)
        rows.append(dict(label=lab(r), hz=float(A(r)["internal_trigger_hz"]), acc=int(A(r)["accumulations"]),
                         width_ns=float(A(r)["gate_width_ns"]), open_s=open_s(r), wall_s=ht.readout_interval_s(r),
                         hits_px=float(h.mean()),
                         hits_px_err=float(h.mean((1, 2)).std(ddof=1) / np.sqrt(len(h))),
                         core_hits=float(h[:, rr <= 6].sum(1).mean()),
                         ann_hits_per_px=float(h[:, (rr >= 10) & (rr <= 16)].mean()), n=len(h)))
    res["gate"] = rows
    X = np.array([[x["open_s"], x["wall_s"], 1.0] for x in rows]); y = np.array([x["hits_px"] for x in rows])
    w = 1 / np.array([max(x["hits_px_err"], 1e-7) for x in rows])
    Xw = X * w[:, None]
    coef, *_ = np.linalg.lstsq(Xw, y * w, rcond=None)
    cov = np.linalg.inv(Xw.T @ Xw)
    chi2 = float((((X @ coef - y) * w) ** 2).sum())
    sc = np.sqrt(max(chi2 / max(len(y) - 3, 1), 1))
    res["gate_fit"] = dict(model="hits/px/readout = k_open*open_s + k_wall*wall_s + c",
                           k_open=float(coef[0]), k_wall=float(coef[1]), c=float(coef[2]),
                           errs=[float(e) for e in np.sqrt(np.diag(cov)) * sc], chi2=chi2, ndof=len(y) - 3)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    s10 = [x for x in rows if x["hz"] == 10000]
    ax[0].errorbar([x["open_s"] for x in s10], [x["hits_px"] for x in s10], [x["hits_px_err"] for x in s10],
                   fmt="o", label="internal 10 kHz x 10k acc (1 s wall)")
    xx = np.linspace(0, .55, 50)
    ax[0].plot(xx, coef[0] * xx + coef[1] * 1.0 + coef[2], "k--", lw=.8, label="joint fit")
    ax[0].set_xlabel("gate-open time per readout (s)  [width 3 ns .. 50 us]")
    ax[0].set_ylabel("hit pixels / pixel / readout (ROI)")
    ax[0].set_title("Gate width at fixed wall time (gain 9000)", fontsize=10)
    ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
    cols = {1000: "C1", 10000: "C0", 100000: "C2"}
    for x in rows:
        closed = x["width_ns"] < 5
        ax[1].errorbar(x["wall_s"], x["hits_px"], x["hits_px_err"], fmt="s" if closed else "o",
                       color=cols[int(x["hz"])], mfc="none" if closed else cols[int(x["hz"])])
    for hz, col in cols.items():
        ax[1].plot([], [], "o", color=col, label=f"internal {hz // 1000} kHz")
    ax[1].plot([], [], "ks", mfc="none", label="3 ns gate (closed)")
    ax[1].plot([], [], "ko", label="gate open (0.1-0.5 s per readout)")
    ax[1].set_xlabel("wall time per readout (s)"); ax[1].set_ylabel("hit pixels / pixel / readout (ROI)")
    ax[1].set_title(f"fit: {coef[0]:.3f} hits/px per s gate-open + {coef[1]:.4f} per s wall time", fontsize=10)
    ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
    fig.tight_layout(); fig.savefig(out / "gate.png", dpi=130); plt.close(fig)

    def closed_partner(r):
        return [q for q in C if A(q)["gate_width_ns"] < 5 and A(q)["internal_trigger_hz"] == A(r)["internal_trigger_hz"]
                and A(q)["accumulations"] == A(r)["accumulations"]]

    beam = []
    for r in C:
        if A(r)["gate_width_ns"] < 5:
            continue
        cl = closed_partner(r)
        net = r["roi_stack"].astype(float).mean(0) - np.mean([q["roi_stack"].astype(float).mean(0) for q in cl], 0)
        rr = radii(r)
        per = np.array([core_minus_annulus(f, rr) for f in r["roi_stack"].astype(float)])
        beam.append(dict(label=lab(r), open_s=open_s(r),
                         core_minus_ann_adu_per_s=core_minus_annulus(net, rr) / open_s(r),
                         err=float(per.std(ddof=1) / np.sqrt(len(per)) / open_s(r)),
                         core_net_adu_px_per_s=float(net[rr < 6].mean() / open_s(r)),
                         ann10_16_net_adu_px_per_s=float(net[(rr >= 10) & (rr <= 16)].mean() / open_s(r)),
                         ring20_40_net_adu_px_per_s=float(net[(rr >= 20) & (rr <= 40)].mean() / open_s(r)),
                         far_net_adu_px_per_s=float(net[rr > 60].mean() / open_s(r))))
    res["beam"] = beam
    try:
        refs = {}
        for k, pth in (("1221_high_pump", REF_1221), ("1225_blocked", REF_1225)):
            q = ht.load_run(ht.DATA_ROOT / pth)
            m = q["roi_stack"].astype(float).mean(0)
            refs[k] = dict(core_minus_ann_adu_per_s=core_minus_annulus(m, radii(q)) / open_s(q), open_s=open_s(q),
                           gain=A(q)["gain"])
        res["beam_ref"] = refs
    except Exception as e:
        res["beam_ref"] = str(e)

    # full frame: longest-open 10 kHz run minus the 10 kHz closed runs
    c10 = [r for r in C if A(r)["internal_trigger_hz"] == 10000]
    lo = max(c10, key=open_s)
    cl = closed_partner(lo)
    net = lo["sum"] / ht.n_rd(lo) - np.mean([q["sum"] / ht.n_rd(q) for q in cl], 0)
    sm = ndi.gaussian_filter(net, 2)
    dog = sm - ndi.gaussian_filter(net, 12)
    bg = 1.4826 * np.median(np.abs(dog - np.median(dog)))
    cand, dd = [], dog.copy()
    dd[:20], dd[-20:], dd[:, :20], dd[:, -20:] = -np.inf, -np.inf, -np.inf, -np.inf
    for _ in range(6):
        y, x = np.unravel_index(np.argmax(dd), dd.shape)
        cand.append(dict(x=int(x), y=int(y), z=float(dog[y, x] / bg), adu_per_readout=float(sm[y, x])))
        dd[max(0, y - 25):y + 25, max(0, x - 25):x + 25] = -np.inf
    res["compact_spot_candidates_fullframe"] = cand
    res["fullframe_net_adu_px_per_s_median"] = float(np.median(net) / open_s(lo))
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.6))
    v = np.percentile(sm, [1, 99.7])
    im = ax[0].imshow(sm, cmap="magma", vmin=v[0], vmax=v[1])
    fig.colorbar(im, ax=ax[0], label="ADU / readout (2 px smoothed)")
    ax[0].plot(*SPOT, "c+", ms=14)
    ax[0].set_title(f"Full frame: {A(lo)['gate_width_ns'] / 1e3:g} us gate x {int(A(lo)['accumulations']):,} "
                    f"({open_s(lo):g} s open) minus 3 ns gate, same wall time\ninternal 10 kHz, gain 9000; "
                    "cyan + = old beam spot (501, 454); green = most compact spots", fontsize=9)
    for c in cand[:3]:
        ax[0].add_patch(plt.Circle((c["x"], c["y"]), 15, ec="lime", fc="none"))
        ax[0].text(c["x"] + 18, c["y"], f"z {c['z']:.0f}", color="lime", fontsize=8)
    x0, y0, n = lo["meta"]["roi"]
    zs = ndi.gaussian_filter(net[y0:y0 + n, x0:x0 + n], 1.5)
    v = np.percentile(zs, [1, 99.7])
    im = ax[1].imshow(zs, cmap="magma", vmin=v[0], vmax=v[1], extent=(x0, x0 + n, y0 + n, y0))
    fig.colorbar(im, ax=ax[1], label="ADU / readout")
    ax[1].add_patch(plt.Circle(SPOT, 6, ec="c", fc="none"))
    ax[1].add_patch(plt.Circle(SPOT, 16, ec="c", fc="none", ls="--"))
    ax[1].set_title("ROI 150x150, same data; circles r = 6 and 16 px at (501, 454)", fontsize=9)
    fig.tight_layout(); fig.savefig(out / "full_frame.png", dpi=130); plt.close(fig)

    # single raw readouts with events marked
    so = [r for r in c10 if abs(A(r)["gate_width_ns"] - 10000) < 1][0]
    scl = [r for r in c10 if A(r)["gate_width_ns"] < 5][0]
    fig, ax = plt.subplots(2, 3, figsize=(13.5, 9))
    for row, r, name in ((0, so, "gate OPEN: 10 us x 10,000 = 0.1 s per readout"),
                         (1, scl, "gate CLOSED: 3 ns x 10,000 (same 1 s wall time)")):
        ab = above(r); x0, y0, n = r["meta"]["roi"]
        for k in range(2):
            f = ab[k]; ev = clusters(f)
            ax[row, k].imshow(f, cmap="gray", vmin=-50, vmax=400, extent=(x0, x0 + n, y0 + n, y0))
            for e in ev:
                ax[row, k].add_patch(plt.Circle((x0 + e[1], y0 + e[0]), 3, ec="r", fc="none", lw=.7))
            ax[row, k].set_title(f"single raw readout #{k + 1}: {len(ev)} events\n{name}", fontsize=8)
        ev = clusters(ab[0])
        ax[row, 2].imshow(ab[0][45:85, 45:85], cmap="gray", vmin=-50, vmax=400,
                          extent=(x0 + 45, x0 + 85, y0 + 85, y0 + 45))
        for e in ev:
            if 45 <= e[0] < 85 and 45 <= e[1] < 85:
                ax[row, 2].add_patch(plt.Circle((x0 + e[1], y0 + e[0]), 2, ec="r", fc="none", lw=1))
        ax[row, 2].set_title("zoom 40x40 px of readout #1 (ADU above pedestal, -50..400)", fontsize=8)
    fig.suptitle("PI-MAX4 single-photon events with the camera's INTERNAL trigger (UNHERALDED: stray/room light + "
                 "intensifier background), gain 9000\nred circles = events (connected pixels >= 60 ADU containing "
                 "a pixel >= 200 ADU above pedestal)", fontsize=10)
    fig.tight_layout(); fig.savefig(out / "single_photons.png", dpi=140); plt.close(fig)

    evo = [e for f in above(so) for e in clusters(f)]
    evc = [e for f in above(scl) for e in clusters(f)]
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    for E, nm in ((evo, "gate open"), (evc, "gate closed")):
        ax[0].hist([e[2] for e in E], bins=np.linspace(150, 2000, 75), histtype="step", label=f"{nm} ({len(E)})")
        ax[1].hist([e[3] for e in E], bins=np.linspace(0, 6000, 61), histtype="step", label=nm)
        ax[2].hist([e[4] for e in E], bins=np.arange(0.5, 40.5), histtype="step", label=nm)
    ax[0].set_xlabel("event peak pixel (ADU above pedestal)"); ax[0].set_yscale("log")
    ax[1].set_xlabel("event summed charge (ADU, pixels >= 60)"); ax[1].set_yscale("log")
    ax[2].set_xlabel("event size (pixels >= 60 ADU)")
    for x_ in ax:
        x_.legend(fontsize=8); x_.grid(alpha=.3)
    fig.suptitle(f"Event pulse heights and sizes, gain 9000, internal 10 kHz, {len(above(so))} readouts each", fontsize=10)
    fig.tight_layout(); fig.savefig(out / "events.png", dpi=130); plt.close(fig)
    res["events"] = {nm: dict(per_readout=len(E) / len(above(r)), peak_median=float(np.median([e[2] for e in E])),
                              sum_median=float(np.median([e[3] for e in E])),
                              size_median=float(np.median([e[4] for e in E])),
                              frac_size_gt20=float(np.mean([e[4] > 20 for e in E])),
                              frac_peak_gt1000=float(np.mean([e[2] > 1000 for e in E])))
                     for nm, E, r in (("open", evo, so), ("closed", evc, scl)) if E}


def closed_bg(runs, out, res):
    """Phase D: closed gate (3 ns) vs wall time per readout; intensifier off; gain 1."""
    D = [r for r in runs if lab(r).startswith("D_")]
    if not D:
        return
    rows = []
    for r in D:
        ab = above(r)
        ev = [e for f in ab[:4] for e in clusters(f)]
        rows.append(dict(label=lab(r), gain=A(r)["gain"], intensifier=r["meta"]["readback"].get("Enable Intensifier"),
                         wall_s=ht.readout_interval_s(r), hits_px=float((ab >= 200).mean()),
                         roi_mean_adu=float(r["roi_stack"].astype(float).mean()),
                         event_size_median=float(np.median([e[4] for e in ev])) if ev else None,
                         event_peak_median=float(np.median([e[2] for e in ev])) if ev else None))
    res["closed_bg"] = rows
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for sel, mk, nm in ((lambda x: x["gain"] == 9000 and x["intensifier"], "o", "gain 9000, intensifier on"),
                        (lambda x: x["gain"] == 9000 and not x["intensifier"], "s", "gain 9000, intensifier OFF"),
                        (lambda x: x["gain"] == 1, "^", "gain 1")):
        R = [x for x in rows if sel(x)]
        if R:
            ax[0].plot([x["wall_s"] for x in R], [x["hits_px"] for x in R], mk, label=nm)
            ax[1].plot([x["wall_s"] for x in R], [x["roi_mean_adu"] for x in R], mk, label=nm)
    on = [x for x in rows if x["gain"] == 9000 and x["intensifier"] and x["hits_px"] > 0]
    if len(on) > 2:
        k, c = np.polyfit(np.log([x["wall_s"] for x in on]), np.log([x["hits_px"] for x in on]), 1)
        res["closed_bg_powerlaw"] = dict(exponent=float(k), hits_px_at_1s=float(np.exp(c)))
        tt = np.linspace(0.2, 4.2, 50); ax[0].plot(tt, np.exp(c) * tt ** k, "k--", lw=.8, label=f"~ t^{k:.2f}")
    ax[0].set_xscale("log"); ax[0].set_yscale("log")
    ax[0].set_xlabel("wall time per readout (s)"); ax[0].set_ylabel("hit pixels / pixel / readout (3 ns gate)")
    ax[0].set_title("Closed-gate background vs time (internal 1 kHz)", fontsize=10)
    ax[1].set_xlabel("wall time per readout (s)"); ax[1].set_ylabel("ROI mean (ADU)")
    ax[1].set_title("Mean level (pedestal + diffuse background)", fontsize=10)
    for x_ in ax:
        x_.legend(fontsize=8); x_.grid(alpha=.3, which="both")
    fig.tight_layout(); fig.savefig(out / "closed_background.png", dpi=130); plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("folders", nargs="+"); p.add_argument("--out", default=None)
    a = p.parse_args()
    runs = [r for f in a.folders for r in ht.load_runs(f)]
    out = pathlib.Path(a.out or pathlib.Path(a.folders[-1]) / "figures")
    out.mkdir(parents=True, exist_ok=True)
    res = {}
    timing(runs, out, res)
    gain_series(runs, out, res)
    gate_series(runs, out, res)
    closed_bg(runs, out, res)
    (out / "health.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(res, indent=1, default=float))
    print("OUT", out)


if __name__ == "__main__":
    main()
