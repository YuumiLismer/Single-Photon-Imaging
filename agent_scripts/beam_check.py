"""Signal-beam check with the EXTERNAL trigger (agent session 5, after Tristan re-coupled the signal arm).

    python agent_scripts/beam_check.py --tag beam --cycles 3 --block 10

Measures the trigger rate (short readouts), picks accumulations so a readout takes <= ~0.8 s, then
interleaves three arms (gain 9000):
  W    = 27 ns + 10 us gate      (beam light inside the gate: unheralded + heralded)
  Woff = 10027 ns + 10 us gate   (same gate, shifted by its own width: in-gate stray/unheralded light)
  C    = 27 ns + 3 ns gate       (closed gate, same wall time: ungated CCD/EM background)
Analysis: full-frame W - C per s gate-open, compact-spot search (difference of Gaussians), core(r<6)
minus annulus(10-16 px) at the old spot (501, 454) and at the brightest compact candidate, in ADU/s
and photons/s (532 ADU per photon event at gain 9000, session 3 H2). Saves figures/beam.json + beam.png.
"""
import argparse, json, sys, pathlib, datetime
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht

ADU_PER_PHOTON = 532.0
p = argparse.ArgumentParser()
p.add_argument("--tag", default="beam"); p.add_argument("--cycles", type=int, default=3)
p.add_argument("--block", type=int, default=10); p.add_argument("--width", type=float, default=10000.0)
p.add_argument("--target-s", type=float, default=0.8, help="target wall-clock time per readout")
p.add_argument("--spot", type=int, nargs=2, default=(501, 454))
a = p.parse_args()
exp = ht.new_experiment(f"agent_{a.tag}"); fig_dir = exp / "figures"
base = dict(gain=9000, trigger_source="External", trigger_threshold_v=0.5, roi_center=tuple(a.spot),
            notes=f"agent {a.tag}", frame_timeout_s=30)
W = a.width
arms = [("W", 27.0, W), ("Woff", 27.0 + W, W), ("C", 27.0, 3.0)]
runs = {k: [] for k, _, _ in arms}
info = {}
cam, limits = ht.open_camera()
try:
    ht.cool_down(cam, setpoint_c=-20.0, poll_s=30)
    rate = ht.measure_trigger_rate(cam, base, accumulations=20_000, n_readouts=4)
    # a gate blocks triggers while open (dead time): gates/s ~ R / (1 + R * width); seen 14:00, 68 kHz -> 40 k gates/s at 10 us
    eff = rate / (1 + rate * (W + 50) * 1e-9)
    acc = int(min(a.target_s * eff, ht.MAX_GATE_OPEN_S / (W * 1e-9)) // 1000 * 1000)
    info.update(rate_hz=rate, acc=acc, start=str(datetime.datetime.now()))
    print(f"RATE {rate / 1e3:.1f} kHz -> {acc} acc", flush=True)
    for c in range(a.cycles):
        order = arms[c % 3:] + arms[:c % 3]
        for name, d, w in order:
            r = ht.take_run(cam, **{**base, "label": f"{a.tag}_{name}_c{c + 1}", "gate_delay_ns": d,
                                    "gate_width_ns": w, "accumulations": acc, "n_readouts": a.block})
            if r is None:
                raise RuntimeError("no readouts (triggers stopped?)")
            if float((r["first"] >= ht.SATURATION_ADU).mean()) > 1e-3:
                raise RuntimeError("saturation > 0.1%: stopping")
            runs[name].append(r)
            dt = ht.readout_interval_s(r)
            print(f"  {name} c{c + 1}: {dt:.2f} s/readout", flush=True)
            if dt > 1.2:
                raise RuntimeError(f"readout interval {dt:.2f} s > 1.2 s: trigger rate dropped, stopping")
    info["rate_end_hz"] = ht.measure_trigger_rate(cam, base, accumulations=20_000, n_readouts=4)
finally:
    ht.shutdown(cam)


def mean_frame(rs):
    return sum(r["sum"].astype(float) for r in rs) / sum(ht.n_rd(r) for r in rs)


def cma(img, cx, cy):
    yy, xx = np.indices(img.shape); rr = np.hypot(xx - cx, yy - cy)
    return float((img[rr < 6] - img[(rr >= 10) & (rr <= 16)].mean()).sum())


open_s = W * 1e-9 * info["acc"]
F = {k: mean_frame(v) for k, v in runs.items() if v}
res = dict(info, open_s=open_s, n_readouts={k: sum(ht.n_rd(r) for r in v) for k, v in runs.items()},
           dt={k: float(np.mean([ht.readout_interval_s(r) for r in v])) for k, v in runs.items() if v})
out = {}
for diff, (A_, B_) in {"W-C": ("W", "C"), "Woff-C": ("Woff", "C"), "W-Woff": ("W", "Woff")}.items():
    if A_ not in F or B_ not in F:
        continue
    net = (F[A_] - F[B_]) / open_s
    sm = ndi.gaussian_filter(net, 2); dog = sm - ndi.gaussian_filter(net, 12)
    bg = 1.4826 * np.median(np.abs(dog - np.median(dog)))
    dd = dog.copy(); dd[:30], dd[-30:], dd[:, :30], dd[:, -30:] = -np.inf, -np.inf, -np.inf, -np.inf
    cand = []
    for _ in range(5):
        y, x = np.unravel_index(np.argmax(dd), dd.shape)
        cand.append(dict(x=int(x), y=int(y), dog_z=float(dog[y, x] / bg), cma_adu_s=cma(net, x, y)))
        dd[max(0, y - 25):y + 25, max(0, x - 25):x + 25] = -np.inf
    # block-to-block scatter of core-minus-annulus at the old spot and the top candidate
    def blocks(cx, cy):
        v = [cma((ra["sum"] / ht.n_rd(ra) - rb["sum"] / ht.n_rd(rb)).astype(float) / open_s, cx, cy)
             for ra, rb in zip(runs[A_], runs[B_])]
        return float(np.mean(v)), float(np.std(v, ddof=1) / np.sqrt(len(v))) if len(v) > 1 else float("nan")
    old = blocks(*a.spot); top = blocks(cand[0]["x"], cand[0]["y"])
    out[diff] = dict(candidates=cand, old_spot_cma_adu_s=old, top_cma_adu_s=top,
                     top_photons_s=top[0] / ADU_PER_PHOTON, old_photons_s=old[0] / ADU_PER_PHOTON,
                     median_adu_px_s=float(np.median(net)))
res["diffs"] = out
(fig_dir / "beam.json").write_text(json.dumps(ht._jsonable(res), indent=1))
print(json.dumps(ht._jsonable(res), indent=1))

fig, ax = plt.subplots(2, 3, figsize=(17, 10.5))
for j, diff in enumerate(["W-C", "Woff-C", "W-Woff"]):
    if diff not in out:
        continue
    A_, B_ = diff.split("-")
    net = ndi.gaussian_filter((F[A_] - F[B_]) / open_s, 2)
    v = np.percentile(net, [1, 99.8])
    h = ax[0, j].imshow(net, cmap="magma", vmin=v[0], vmax=v[1]); fig.colorbar(h, ax=ax[0, j], label="ADU/px per s gate-open")
    ax[0, j].plot(*a.spot, "c+", ms=14)
    for c_ in out[diff]["candidates"][:3]:
        ax[0, j].add_patch(plt.Circle((c_["x"], c_["y"]), 15, ec="lime", fc="none"))
    t = out[diff]["top_cma_adu_s"]
    ax[0, j].set_title(f"{diff} full frame; top ({out[diff]['candidates'][0]['x']},{out[diff]['candidates'][0]['y']}) "
                       f"core-ann {t[0]:,.0f}+-{t[1]:,.0f} ADU/s\n(~{t[0] / ADU_PER_PHOTON:,.0f} photons/s); old spot "
                       f"{out[diff]['old_spot_cma_adu_s'][0]:,.0f} ADU/s", fontsize=9)
    cx, cy = out["W-C"]["candidates"][0]["x"], out["W-C"]["candidates"][0]["y"]
    z = net[max(0, cy - 75):cy + 75, max(0, cx - 75):cx + 75]
    v = np.percentile(z, [1, 99.9])
    h = ax[1, j].imshow(z, cmap="magma", vmin=v[0], vmax=v[1], extent=(max(0, cx - 75), cx + 75, cy + 75, max(0, cy - 75)))
    fig.colorbar(h, ax=ax[1, j]); ax[1, j].set_title(f"{diff}: 150x150 around ({cx},{cy})", fontsize=9)
fig.suptitle(f"Beam check {res['start'][:16]}: External trigger {info['rate_hz'] / 1e3:.1f} kHz "
             f"(end {info.get('rate_end_hz', float('nan')) / 1e3:.1f}), gain 9000, {info['acc']:,} acc x 10 us = "
             f"{open_s:.2f} s open/readout; W = 27 ns+10 us, Woff = {27 + W:g} ns+10 us, C = 3 ns gate", fontsize=10)
fig.tight_layout(); fig.savefig(fig_dir / "beam.png", dpi=110)
print("EXP", exp)
