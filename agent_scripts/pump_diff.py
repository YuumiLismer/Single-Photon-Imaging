"""Where does the signal light land? Pump-dependent light = (W - C) at high herald rate minus (W - C) at low rate,
both per second of gate-open time, from two beam_check.py experiments (no camera).
    python agent_scripts/pump_diff.py HIGH_EXP LOW_EXP
    python agent_scripts/pump_diff.py --pairs HI1,LO1 HI2,LO2 ... [--bin 64]   # alternating runs: per-pair differences,
        binned into bin x bin px blocks, mean +- error over pairs (z map), so room-light drift averages out"""
import sys, pathlib, json, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter


def w_minus_c(exp):
    per = {}
    for r in ht.load_runs(exp):
        a = r["meta"]["actual"]; lab = r["meta"]["settings"]["label"]
        kind = "C" if a["gate_width_ns"] < 10 else ("Woff" if a["gate_delay_ns"] > 5000 else "W")
        n = ht.n_rd(r)
        per.setdefault(kind, []).append((r["sum"] / n, a["accumulations"] * a["gate_width_ns"] * 1e-9))
    img = lambda k: np.median([s for s, _ in per[k]], 0)
    open_s = np.median([o for _, o in per["W"]])
    return (img("W") - img("C")) / open_s, (img("Woff") - img("C")) / open_s


if sys.argv[1] == "--pairs":
    args = sys.argv[2:]; b = 64
    if "--bin" in args:
        i = args.index("--bin"); b = int(args[i + 1]); args = args[:i] + args[i + 2:]
    diffs = []
    for pr in args:
        h, l = pr.split(",")
        diffs.append(w_minus_c(h)[0] - w_minus_c(l)[0])
    D = np.array(diffs)
    nb = 1024 // b
    B = D[:, :nb * b, :nb * b].reshape(len(D), nb, b, nb, b).mean((2, 4))
    B = B - np.median(B, axis=(1, 2), keepdims=True)          # drop the uniform room-light drift per pair
    m, e = B.mean(0), B.std(0, ddof=1) / np.sqrt(len(B))
    z = m / e
    tot = D.mean((1, 2))
    print(f"{len(D)} pairs; full-frame pump-dependent mean per pair: " + ", ".join(f"{t:+.1f}" for t in tot) + " ADU/px/s")
    iy, ix = np.unravel_index(np.argmax(z), z.shape)
    print(f"bin {b} px: max z {z.max():+.1f} at block ({ix * b + b // 2}, {iy * b + b // 2}), {m[iy, ix]:+.1f} +- {e[iy, ix]:.1f} ADU/px/s;"
          f" blocks with z > 3: {int((z > 3).sum())} of {z.size}")
    fig, ax = plt.subplots(1, 2, figsize=(12, 5.2))
    h0 = ax[0].imshow(m, cmap="RdBu_r", vmin=-np.abs(m).max(), vmax=np.abs(m).max(), extent=[0, nb * b, nb * b, 0])
    plt.colorbar(h0, ax=ax[0], label="pump-dependent ADU/px per s gate-open (median removed)")
    ax[0].set_title(f"mean over {len(D)} high-low pairs, {b} px blocks")
    h1 = ax[1].imshow(z, cmap="RdBu_r", vmin=-5, vmax=5, extent=[0, nb * b, nb * b, 0])
    plt.colorbar(h1, ax=ax[1], label="z (mean / error over pairs)"); ax[1].set_title("significance")
    out = pathlib.Path(args[0].split(",")[0]) / "figures" / f"pump_pairs_bin{b}.png"
    fig.suptitle("Pump-dependent light (W - C, high minus low herald rate), alternating runs"); fig.tight_layout(); fig.savefig(out, dpi=100)
    print(out); sys.exit()
hi, lo = sys.argv[1], sys.argv[2]
(wh, woh), (wl, wol) = w_minus_c(hi), w_minus_c(lo)
d = gaussian_filter(wh - wl, 8); doff = gaussian_filter(woh - wol, 8)
fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
for a_, im, t in ((ax[0], gaussian_filter(wh, 8), "W - C at high pump"), (ax[1], d, "pump-dependent: (W-C)high - (W-C)low"),
                  (ax[2], doff, "same for Woff (late gate, 10 us after herald)")):
    h = a_.imshow(im, cmap="magma" if a_ is ax[0] else "RdBu_r",
                  **({} if a_ is ax[0] else dict(vmin=-np.abs(d).max(), vmax=np.abs(d).max())))
    a_.set_title(t); plt.colorbar(h, ax=a_, label="ADU/px per s gate-open (8 px smooth)")
y, x = np.unravel_index(np.argmax(d), d.shape)
ax[1].plot(x, y, "c+", ms=16)
fig.suptitle(f"{pathlib.Path(hi).name} minus {pathlib.Path(lo).name}; max pump-dependent at ({x}, {y})")
fig.tight_layout(); out = pathlib.Path(hi) / "figures" / "pump_diff.png"; fig.savefig(out, dpi=100)
yy, xx = np.mgrid[:1024, :1024]
for name, im in (("W", wh - wl), ("Woff", woh - wol)):
    core = im[np.hypot(xx - x, yy - y) < 60].mean(); ring = im[(np.hypot(xx - x, yy - y) > 150) & (np.hypot(xx - x, yy - y) < 250)].mean()
    print(f"{name}: pump-dependent full-frame mean {im.mean():+.1f} ADU/px/s | core(r<60 at max) {core:+.1f} vs ring {ring:+.1f}")
print("max at", x, y, "->", out)
