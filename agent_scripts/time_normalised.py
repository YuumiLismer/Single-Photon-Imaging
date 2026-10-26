"""Re-analyse a target run with the background modelled per second of readout wall time.
Per readout: spot score M (beam template minus 4 control apertures, as target_image.py) and the
readout duration dt. Model: M = S*[ON] + b*dt (b fitted on OFF readouts, through the origin), and
for the image: hits(pixel) = s(pixel)*[ON] + rate(pixel)*dt."""
import sys, pathlib, json, numpy as np
sys.path.insert(0, r"C:/Users/ReschLab/Single-Photon-Imaging"); sys.path.insert(0, r"C:/Users/ReschLab/Single-Photon-Imaging/agent_scripts")
import hoi_troubleshoot as ht
from analyze_ab import beam_template
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
exp, SPOT, TEMPLATE = sys.argv[1], (462, 460), sys.argv[2]
rs = ht.load_runs(exp)
ons = [r for r in rs if r["meta"]["settings"]["kind"] == "signal"]
offs = [r for r in rs if r["meta"]["settings"]["kind"] == "dark_offdelay"]
w, _ = beam_template(ons[0]["meta"]["roi"], SPOT, template_exp=TEMPLATE)
eff = (w * w / w.sum()).sum()
ctrl = [np.roll(np.roll(w, dy, 0), dx, 1) for dx, dy in ((45, 45), (-45, 45), (45, -45), (-45, -45))]
def per_readout(rr):
    M, DT, H = [], [], 0
    for r in rr:
        hit = ht.above_pedestal(r) >= 200
        t = np.asarray(r["t"], float); dt = np.diff(t); dt = np.concatenate([[np.median(dt)], dt])[:len(hit)]
        M.append(((hit * w).sum((1, 2)) - np.mean([(hit * c).sum((1, 2)) for c in ctrl], 0)) / eff)
        spot_raw = (hit * w).sum((1, 2)) / eff
        DT.append(dt); H = H + hit.sum(0)
    return np.concatenate(M), np.concatenate(DT), H
Mon, Ton, Hon = per_readout(ons); Moff, Toff, Hoff = per_readout(offs)
print(f"ON  {len(Mon)} readouts, mean dt {Ton.mean():.3f} s | OFF {len(Moff)} readouts, mean dt {Toff.mean():.3f} s")
naive = Mon.mean() - Moff.mean()
b = (Moff * Toff).sum() / (Toff * Toff).sum()          # OFF: M = b*dt
res = Mon - b * Ton
S = res.mean(); err = np.hypot(res.std(ddof=1) / np.sqrt(len(res)), Ton.mean() * np.std(Moff - b * Toff, ddof=1) / np.sqrt(len(Moff)) / Toff.mean())
print(f"naive per-readout ON-OFF: {naive:+.2f} | OFF spot score per s b = {b:+.3f}/s")
print(f"time-normalised net heralded spot hits/readout: {S:+.2f} +- {err:.2f}  z {S/err:+.1f}")
# image: per-pixel background rate from OFF, scaled to ON wall time
img = Hon - Hoff * (Ton.sum() / Toff.sum())
img_n = img / len(Mon)
naive_img = Hon / len(Mon) - Hoff / len(Moff)
roi = ons[0]["meta"]["roi"]; x0, y0 = roi[0], roi[1]
ext = [x0, x0 + 150, y0 + 150, y0]
fig, ax = plt.subplots(1, 2, figsize=(12, 5.2))
for a_, im, title in ((ax[0], naive_img, f"per-readout ON - OFF (old)\nnet {naive:+.2f}/readout"),
                      (ax[1], img_n, f"ON - OFF x (ON wall time / OFF wall time)\nnet {S:+.2f} +- {err:.2f}/readout, z {S/err:+.1f}")):
    sm = gaussian_filter(im, 2); v = np.abs(sm).max()
    h = a_.imshow(sm, cmap="RdBu_r", vmin=-v, vmax=v, extent=ext); a_.set_title(title); plt.colorbar(h, ax=a_, label="hits/px/readout (2 px smooth)")
    a_.set_xlabel("x (px)"); a_.set_ylabel("y (px)")
m0 = ons[0]["meta"]["actual"]
fig.suptitle(f"{pathlib.Path(exp).name}: gate {m0['gate_delay_ns']:.0f} ns + {m0['gate_width_ns']:.0f} ns, {len(Mon)} ON + {len(Moff)} OFF readouts; background scales with readout wall time")
fig.tight_layout(); out = pathlib.Path(exp) / "figures" / "heralded_time_normalised.png"; fig.savefig(out, dpi=120); print(out)
