"""Locate an alignment diode / laser spot on the sensor, light-safely (same settings as live_view.py:
Internal trigger at 10 kHz, emICCD gain 1, short gates). Starts at the shortest gate and 1 accumulation and
doubles the exposure (accumulations first, then gate) until the spot peak is well above the background;
backs off as soon as more than 25 pixels saturate. Saves the frame and reports the spot centroid and size.

    python agent_scripts/diode_locate.py [--tag diode]
"""
import argparse
import json
import pathlib
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht
import camera_control as cc

SAT_ADU, TRIGGER_HZ, MAX_ACC = 60_000, 10_000, 5_000
p = argparse.ArgumentParser()
p.add_argument("--tag", default="diode")
p.add_argument("--peak-adu", type=float, default=8_000, help="stop once the smoothed spot peak is this far above the median")
a = p.parse_args()

folder = ht.new_experiment(f"agent_{a.tag}_locate") / "figures"
folder.mkdir(parents=True, exist_ok=True)
cam, limits = ht.open_camera()
w_min = max(3.0, float(limits["width_min_ns"]))
base = {**ht.DEFAULT_SETTINGS, "label": "diode_locate", "gain": 1, "trigger_source": "Internal",
        "internal_trigger_hz": TRIGGER_HZ, "gate_delay_ns": 100.0, "photon_mode": "Disabled"}
gate, acc, img, log = w_min, 1, None, []


def grab(gate, acc):
    try:
        cam.stop_acquisition()
    except Exception:
        pass
    ht.configure(cam, {**base, "gate_width_ns": gate, "accumulations": acc}, verbose=False)
    cam.start_acquisition()
    frames, t0 = [], time.time()
    while len(frames) < 3 and time.time() - t0 < 10:          # drop the first frame (settings change)
        if cc.poll_frame(cam, 0.5):
            frames += cam.read_multiple_images()
    cam.stop_acquisition()
    return np.asarray(frames[-1], float)


try:
    for _ in range(30):
        img = grab(gate, acc)
        n_sat = int((img >= SAT_ADU).sum())
        sm = gaussian_filter(img - np.median(img), 3)
        peak = float(sm.max())
        log.append(dict(gate_ns=gate, acc=acc, n_sat=n_sat, peak_above_median=peak, max_adu=float(img.max())))
        print(f"gate {gate:.1f} ns x {acc} acc: max {img.max():.0f} ADU, smoothed peak +{peak:.0f}, saturated px {n_sat}", flush=True)
        if n_sat > 25:
            if acc > 1:
                acc //= 2
            elif gate > w_min * 1.01:
                gate = max(w_min, gate / 2)
            else:
                print("SATURATING at the shortest exposure: stopping (dim the diode)"); break
            img = grab(gate, acc); break
        if peak >= a.peak_adu:
            break
        if acc < MAX_ACC:
            acc = min(MAX_ACC, acc * 2)
        elif gate < 10_000:
            gate = min(10_000, gate * 2)
        else:
            break
finally:
    try:
        cam.stop_acquisition()
    except Exception:
        pass
    ht.shutdown(cam)

bg = np.median(img)
sm = gaussian_filter(img - bg, 3)
y0, x0 = np.unravel_index(np.argmax(sm), sm.shape)
mask = sm > 0.5 * sm[y0, x0]                                   # half-maximum region around the peak
from scipy.ndimage import label
lab, _ = label(mask)
reg = lab == lab[y0, x0]
w = np.clip(img - bg, 0, None) * reg
yy, xx = np.mgrid[:img.shape[0], :img.shape[1]]
cx, cy = float((w * xx).sum() / w.sum()), float((w * yy).sum() / w.sum())
sx, sy = float(np.sqrt((w * (xx - cx) ** 2).sum() / w.sum())), float(np.sqrt((w * (yy - cy) ** 2).sum() / w.sum()))
res = dict(centroid_xy=[round(cx, 1), round(cy, 1)], peak_xy=[int(x0), int(y0)], rms_size_xy_px=[round(sx, 1), round(sy, 1)],
           fwhm_area_px=int(reg.sum()), peak_above_median_adu=float(sm[y0, x0]), median_adu=float(bg),
           gate_ns=gate, acc=acc, steps=log)
(folder / "diode_locate.json").write_text(json.dumps(res, indent=1))
np.save(folder / "diode_frame.npy", img)
fig, ax = plt.subplots(1, 2, figsize=(14, 6.5))
lo, hi = np.percentile(img, [1, 99.95])
ax[0].imshow(img, cmap="gray", vmin=lo, vmax=hi); ax[0].plot(cx, cy, "c+", ms=20, mew=1.5)
ax[0].axhline(512, color="y", lw=0.5); ax[0].axvline(512, color="y", lw=0.5)
ax[0].set_title(f"full frame, gain 1, gate {ht.fmt_ns(gate)} x {acc} acc; spot centroid ({cx:.0f}, {cy:.0f})")
r = int(max(40, 4 * max(sx, sy)))
sl = (slice(max(0, int(cy) - r), int(cy) + r), slice(max(0, int(cx) - r), int(cx) + r))
ax[1].imshow(img[sl], cmap="magma", extent=[sl[1].start, sl[1].stop, sl[0].stop, sl[0].start])
ax[1].contour(xx[sl], yy[sl], reg[sl], [0.5], colors="c", linewidths=0.8)
ax[1].set_title(f"zoom: rms size {sx:.1f} x {sy:.1f} px, half-max area {reg.sum()} px")
fig.tight_layout(); fig.savefig(folder / "diode_locate.png", dpi=100)
print(json.dumps({k: v for k, v in res.items() if k != "steps"}))
print(folder / "diode_locate.png")
