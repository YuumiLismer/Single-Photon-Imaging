"""Pump-dependence of the light on the sensor, from two or more beam_check.py experiments (no camera).

    python agent_scripts/pump_dependence.py LOW_EXP HIGH_EXP [MORE_EXP ...] --out FIG_DIR

SPDC signal light (heralded or not) scales with pump power, stray room/setup light does not. For each
beam_check experiment the arm frames (W = 27 ns + 10 us, Woff = 10027 ns + 10 us, C = 3 ns closed
gate) are combined as the per-pixel MEDIAN of the per-run mean frames (robust to single-readout
transients), converted to ADU/px per s gate-open, and differenced HIGH - LOW for W-C and W-Woff.
max|z| is the largest |value| of the 2-px smoothed difference in units of its robust (MAD) spread. Prints and saves the top compact candidates (difference of Gaussians),
core(r<6)-minus-annulus(10-16) at --spot, and the global median ratio of the stray scene.
Writes pump_dependence.png / .json into --out.
"""
import argparse, json, sys, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht

p = argparse.ArgumentParser()
p.add_argument("exps", nargs="+"); p.add_argument("--out", required=True)
p.add_argument("--spot", type=int, nargs=2, default=(501, 454))
a = p.parse_args()
out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)


def arms(exp):
    exp = pathlib.Path(exp)
    js = exp / "figures" / "beam.json"
    info = json.loads(js.read_text()) if js.exists() else {}
    rs = {"W": [], "Woff": [], "C": []}
    for r in ht.load_runs(exp / "runs"):
        lab = pathlib.Path(r["meta"]["file"]).stem
        for k in rs:
            if f"_{k}_c" in lab:
                rs[k].append(r)
    acc = int(rs["W"][0]["meta"]["actual"]["accumulations"]); wid = float(rs["W"][0]["meta"]["actual"]["gate_width_ns"])
    open_s = acc * wid * 1e-9
    F = {}
    for k, v in rs.items():
        per = np.array([r["sum"].astype(float) / ht.n_rd(r) for r in v]) / open_s
        F[k] = np.median(per, axis=0) if len(per) >= 3 else per.mean(0)
    return dict(name=exp.name, rate=info.get("rate_hz", np.nan), F=F, open_s=open_s)


def cma(img, cx, cy):
    yy, xx = np.indices(img.shape); rr = np.hypot(xx - cx, yy - cy)
    return float((img[rr < 6] - img[(rr >= 10) & (rr <= 16)].mean()).sum())


def cands(net, n=5):
    sm = ndi.gaussian_filter(net, 2); dog = sm - ndi.gaussian_filter(net, 12)
    bg = 1.4826 * np.median(np.abs(dog - np.median(dog)))
    dd = dog.copy(); dd[:30], dd[-30:], dd[:, :30], dd[:, -30:] = -np.inf, -np.inf, -np.inf, -np.inf
    c = []
    for _ in range(n):
        y, x = np.unravel_index(np.argmax(dd), dd.shape)
        c.append(dict(x=int(x), y=int(y), dog_z=float(dog[y, x] / bg), cma_adu_s=cma(net, x, y)))
        dd[max(0, y - 25):y + 25, max(0, x - 25):x + 25] = -np.inf
    return c


E_ = [arms(e) for e in a.exps]
E_.sort(key=lambda d: d["rate"])
lo = E_[0]
res = {"experiments": [dict(name=d["name"], rate_hz=d["rate"]) for d in E_], "pairs": []}
fig, ax = plt.subplots(len(E_) - 1, 3, figsize=(18, 5.6 * (len(E_) - 1)), squeeze=False)
for i, hi in enumerate(E_[1:]):
    pr = dict(low=lo["name"], high=hi["name"], low_rate=lo["rate"], high_rate=hi["rate"])
    for j, (A, B) in enumerate((("W", "C"), ("W", "Woff"))):
        d = (hi["F"][A] - hi["F"][B]) - (lo["F"][A] - lo["F"][B])
        sm = ndi.gaussian_filter(d, 2)
        z = (sm - np.median(sm)) / (1.4826 * np.median(np.abs(sm - np.median(sm))))  # robust global z of the 2-px smoothed map
        c = cands(d)
        key = f"{A}-{B}"
        pr[key] = dict(candidates=c, spot_cma_adu_s=cma(d, *a.spot), max_abs_z_smoothed=float(np.abs(z[30:-30, 30:-30]).max()))
        v = np.percentile(ndi.gaussian_filter(d, 2), [1, 99.8])
        h = ax[i, j].imshow(ndi.gaussian_filter(d, 2), cmap="RdBu_r", vmin=-max(abs(v)), vmax=max(abs(v)))
        fig.colorbar(h, ax=ax[i, j], label="ADU/px per s gate-open")
        ax[i, j].plot(*a.spot, "k+", ms=14)
        for cc in c[:3]:
            ax[i, j].add_patch(plt.Circle((cc["x"], cc["y"]), 15, ec="lime", fc="none"))
        ax[i, j].set_title(f"({key}) at {hi['rate'] / 1e3:.0f} kHz minus at {lo['rate'] / 1e3:.0f} kHz\n"
                           f"top ({c[0]['x']},{c[0]['y']}) DoG z {c[0]['dog_z']:.1f}; spot {a.spot}: "
                           f"{pr[key]['spot_cma_adu_s']:,.0f} ADU/s", fontsize=9)
    # stray-scene ratio: W-C high / low over the brightest 2 % pixels (smoothed)
    wl = ndi.gaussian_filter(lo["F"]["W"] - lo["F"]["C"], 4); wh = ndi.gaussian_filter(hi["F"]["W"] - hi["F"]["C"], 4)
    m = wl > np.percentile(wl, 98)
    pr["stray_ratio_high_over_low"] = float(np.median(wh[m] / wl[m]))
    sc = ax[i, 2]; sc.scatter(wl[::8, ::8].ravel(), wh[::8, ::8].ravel(), s=1, alpha=0.3)
    mx = np.percentile(wl, 99.9); sc.plot([0, mx], [0, mx], "k--", lw=1)
    sc.set_xlabel(f"W-C at {lo['rate'] / 1e3:.0f} kHz (ADU/px/s, 4-px smoothed)")
    sc.set_ylabel(f"W-C at {hi['rate'] / 1e3:.0f} kHz"); sc.set_title(
        f"per-pixel scene brightness, high vs low pump; bright-pixel ratio {pr['stray_ratio_high_over_low']:.3f}\n"
        f"(SPDC light would lie above the dashed 1:1 line by ~{hi['rate'] / lo['rate']:.1f}x)", fontsize=9)
    res["pairs"].append(pr)
fig.suptitle("Pump dependence of the light on the sensor (beam_check experiments; median of runs per arm)", fontsize=11)
fig.tight_layout(); fig.savefig(out / "pump_dependence.png", dpi=100)
(out / "pump_dependence.json").write_text(json.dumps(ht._jsonable(res), indent=1))
for pr in res["pairs"]:
    print(f"{pr['high']} ({pr['high_rate'] / 1e3:.0f} kHz) - {pr['low']} ({pr['low_rate'] / 1e3:.0f} kHz): "
          f"stray ratio {pr['stray_ratio_high_over_low']:.3f}")
    for key in ("W-C", "W-Woff"):
        q = pr[key]
        print(f"  {key}: spot {q['spot_cma_adu_s']:,.0f} ADU/s, max|z| {q['max_abs_z_smoothed']:.1f}, top",
              [(c["x"], c["y"], round(c["dog_z"], 1), round(c["cma_adu_s"])) for c in q["candidates"][:4]])
print("OUT", out)
