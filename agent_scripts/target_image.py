"""Heralded single-photon target image (agent session 5): long interleaved run at the found gate.

    python agent_scripts/target_image.py --delay 542 --width 40 --acc 80000 --cycles 4 --block 50

ON = take_run at the gate, OFF = ht.take_offdelay_background (same gate +5 us), alternated in blocks
so drifts (herald rate, background) hit both alike; blocks are merged per arm. Saves ht.plot_run,
ht.plot_raw_vs_accumulated, ht.plot_wiener, ht.threshold_sweep and a clean 'heralded single-photon
image' (ON-OFF photon hits per readout, ROI 150x150) with the significance of the spot (beam-profile
matched filter minus 4 control apertures, Welch errors, as window_figures.py). Stops if the readout
interval exceeds --max-readout-s (Tristan: readouts <= ~1 s wall clock).
Re-analyse a finished run with --analyse <experiment folder> (no camera).
"""
import argparse, json, sys, pathlib, datetime
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import hoi_troubleshoot as ht
from analyze_ab import beam_template

p = argparse.ArgumentParser()
p.add_argument("--delay", type=float, default=None); p.add_argument("--width", type=float, default=None)
p.add_argument("--acc", type=int, default=80000); p.add_argument("--cycles", type=int, default=4)
p.add_argument("--block", type=int, default=50); p.add_argument("--thr", type=float, default=200)
p.add_argument("--max-readout-s", type=float, default=1.05); p.add_argument("--tag", default="target")
p.add_argument("--analyse", default=None, help="only re-analyse this experiment folder")
a = p.parse_args()
SPOT = (501, 454)


def merge(rs):
    m = dict(rs[0]); m["meta"] = json.loads(json.dumps(rs[0]["meta"]))
    m["sum"] = sum(r["sum"].astype(float) for r in rs)
    for k in ("roi_stack", "roi_totals", "full_totals"):
        m[k] = np.concatenate([r[k] for r in rs])
    m["t"] = np.concatenate([np.asarray(r["t"], float) - np.asarray(r["t"], float)[0] + 1e4 * i for i, r in enumerate(rs)])
    m["meta"]["n_readouts_done"] = sum(ht.n_rd(r) for r in rs)
    return m


if a.analyse:
    exp = pathlib.Path(a.analyse); ht.set_experiment(exp)
    rs = ht.load_runs(str(exp))
    ons = [r for r in rs if r["meta"]["settings"]["kind"] == "signal"]
    offs = [r for r in rs if r["meta"]["settings"]["kind"] == "dark_offdelay"]
else:
    exp = ht.new_experiment(f"agent_{a.tag}")
    base = dict(gain=9000, trigger_source="External", trigger_threshold_v=0.5, roi_center=(500, 455), kind="signal",
                gate_delay_ns=a.delay, gate_width_ns=a.width, accumulations=a.acc, n_readouts=a.block,
                notes="agent session 5 target: heralded single-photon image")
    ons, offs = [], []
    cam, limits = ht.open_camera()
    try:
        ht.cool_down(cam, setpoint_c=-20.0, poll_s=30)
        for c in range(a.cycles):
            for arm in (("on", "off") if c % 2 == 0 else ("off", "on")):
                if arm == "on":
                    r = ht.take_run(cam, **{**base, "label": f"{a.tag}_on_c{c + 1}"})
                else:
                    r = ht.take_offdelay_background(cam, {**base, "label": f"{a.tag}_c{c + 1}"}, shift_ns=5000.0,
                                                    label=f"{a.tag}_off_c{c + 1}")
                if r is None:
                    raise RuntimeError("no readouts")
                if float((r["first"] >= ht.SATURATION_ADU).mean()) > 1e-3:
                    raise RuntimeError("saturation > 0.1%: stopping")
                (ons if arm == "on" else offs).append(r)
                dt = ht.readout_interval_s(r)
                print(f"  {arm} c{c + 1}: {dt:.2f} s/readout", flush=True)
                if dt > a.max_readout_s:
                    raise RuntimeError(f"readout interval {dt:.2f} s > {a.max_readout_s} s: stopping")
    finally:
        ht.shutdown(cam)

S, D = merge(ons), merge(offs)
w, rr = beam_template(S["meta"]["roi"], SPOT)
eff = (w * w / w.sum()).sum()
ctrl = [np.roll(np.roll(w, dy, 0), dx, 1) for dx, dy in ((45, 45), (-45, 45), (45, -45), (-45, -45))]


def hits_all(rs):
    H, M = 0, []
    for r in rs:
        st = ht.above_pedestal(r); hit = (st - np.median(st, 0)) >= a.thr
        H = H + hit.sum(0)
        M.append(((hit * w).sum((1, 2)) - np.mean([(hit * c).sum((1, 2)) for c in ctrl], 0)) / eff)
    M = np.concatenate(M)
    return H / len(M), M


Hon, Mon = hits_all(ons)
Hoff, Moff = hits_all(offs)
net = Mon.mean() - Moff.mean()
err = np.hypot(Mon.std(ddof=1) / np.sqrt(len(Mon)), Moff.std(ddof=1) / np.sqrt(len(Moff)))
act = S["meta"]["actual"]
acc = int(act["accumulations"]); delay = float(act["gate_delay_ns"]); width = float(act["gate_width_ns"])
dts = [ht.readout_interval_s(r) for r in ons]
res = dict(delay=delay, width=width, acc=acc, readouts_on=len(Mon), readouts_off=len(Moff),
           net_spot_hits_per_readout=net, err=err, z=net / err, per_trigger=net / acc,
           on_spot=Mon.mean(), off_spot=Moff.mean(), gates_per_s=float(acc / np.median(dts)),
           s_per_readout_on=float(np.median(dts)), heralded_spot_hits_per_s_wall=float(net / np.median(dts)))
print(json.dumps(ht._jsonable(res), indent=1))
fig_dir = exp / "figures"; fig_dir.mkdir(exist_ok=True, parents=True)
(fig_dir / "target.json").write_text(json.dumps(ht._jsonable(res), indent=1))

for fn, kw in ((ht.plot_run, dict(dark=D, spot=SPOT, save_as="target_plot_run")),
               (ht.plot_raw_vs_accumulated, dict(dark=D, save_as="target_raw_vs_sum")),
               (ht.plot_wiener, dict(dark=D, save_as="target_wiener"))):
    try:
        out = fn(S, **kw); plt.close("all")
        if fn is ht.plot_wiener:
            res["wiener_snr_before_after"] = out
    except Exception as e:
        print(f"  {fn.__name__} failed: {type(e).__name__}: {e}")
try:
    df, best = ht.threshold_sweep(S, D, spot=SPOT, save_as="target_threshold_sweep"); plt.close("all")
    res["best_threshold_adu"] = best
    print(df.to_string()); print("best threshold", best)
except Exception as e:
    print(f"  threshold_sweep failed: {type(e).__name__}: {e}")
(fig_dir / "target.json").write_text(json.dumps(ht._jsonable(res), indent=1))

x0, y0, n = S["meta"]["roi"]
ext = (x0 - .5, x0 + n - .5, y0 + n - .5, y0 - .5)
fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.4))
vmax = np.percentile(gaussian_filter(Hon, 1), 99.7)
panels = ((ax[0], Hon, f"ON gate {delay:g}-{delay + width:g} ns after the herald"),
          (ax[1], Hoff, f"OFF gate {delay + 5000:g}-{delay + 5000 + width:g} ns (off-delay background)"))
for a_, im, t in panels:
    h = a_.imshow(gaussian_filter(im, 1), cmap="inferno", vmin=0, vmax=vmax, extent=ext)
    a_.set_title(t + "\nphoton hits per pixel per readout (1 px smooth)", fontsize=10)
    fig.colorbar(h, ax=a_, fraction=0.046)
dimg = gaussian_filter(Hon - Hoff, 2); lim = np.abs(dimg).max()
h = ax[2].imshow(dimg, cmap="RdBu_r", vmin=-lim, vmax=lim, extent=ext)
ax[2].set_title(f"HERALDED single photons = ON - OFF (2 px smooth)\nnet {net:+.2f} +- {err:.2f} spot hits/readout, "
                f"z = {net / err:+.1f}", fontsize=10)
fig.colorbar(h, ax=ax[2], fraction=0.046)
for a_ in ax:
    a_.set_xlabel("x (px)"); a_.set_ylabel("y (px)")
fig.suptitle(f"Heralded single-photon image, PI-MAX4 emICCD, {S['meta']['started'][:16]}: External herald trigger "
             f"(~{res['gates_per_s'] / 1e3:.0f} k gates/s), gain 9000, gate {delay:g} ns + {width:g} ns, {acc:,} acc/readout, "
             f"{len(Mon)} ON + {len(Moff)} OFF readouts, hit threshold {a.thr:g} ADU, ROI 150x150", fontsize=10)
fig.tight_layout(); fig.savefig(fig_dir / "heralded_single_photon_image.png", dpi=130)
print("EXP", exp)
