"""One uncropped, unsubtracted full-sensor picture with a wide gate (Tristan: "show me one picture with
large gate width, don't crop anything"). External herald trigger, gain 9000 by default; accumulations
are chosen from the measured herald rate so each readout takes ~--target-s and gate x accumulations
stays <= 1 s.

    python agent_scripts/full_frame_snapshot.py [--delay 27] [--width 10000] [--readouts 10]
"""
import argparse
import pathlib
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht

p = argparse.ArgumentParser()
p.add_argument("--delay", type=float, default=27.0, help="gate delay (ns)")
p.add_argument("--width", type=float, default=10_000.0, help="gate width (ns), <= 10,000")
p.add_argument("--gain", type=int, default=9000)
p.add_argument("--readouts", type=int, default=10)
p.add_argument("--target-s", type=float, default=0.8, help="wall-clock seconds per readout")
p.add_argument("--tag", default="")
a = p.parse_args()
assert a.width <= 10_000 and a.gain <= 9000

folder = ht.new_experiment("agent_full_frame" + (f"_{a.tag}" if a.tag else "")) / "figures"
base = {**ht.DEFAULT_SETTINGS, "gain": a.gain, "trigger_source": "External", "photon_mode": "Disabled"}
cam, _ = ht.open_camera()
try:
    rate = ht.measure_trigger_rate(cam, {**base, "gain": 1}, accumulations=100_000, n_readouts=2)
    acc = int(min(rate * a.target_s * 0.9, 1.0 / (a.width * 1e-9)))
    r = ht.take_run(cam, **{**base, "label": "full_frame_snapshot", "gate_delay_ns": a.delay,
                            "gate_width_ns": a.width, "accumulations": acc, "n_readouts": a.readouts})
finally:
    ht.shutdown(cam)

img = r["sum"] / a.readouts                                  # mean ADU per readout, every pixel
act = r["meta"]["actual"]
lo, hi = np.percentile(img, [0.5, 99.9])
fig, ax = plt.subplots(figsize=(10, 9.2))
im = ax.imshow(img, cmap="magma", vmin=lo, vmax=hi, origin="upper", interpolation="nearest")
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label="ADU per readout (raw, no background subtracted)")
ax.set_xlabel("x (px)")
ax.set_ylabel("y (px)")
ax.set_title(f"Full sensor 1024x1024, raw, uncropped: gate {act['gate_delay_ns']:g} ns + {act['gate_width_ns'] / 1e3:g} us after "
             f"each herald\n{act['accumulations']:,} gates/readout x {a.readouts} readouts, gain {act['gain']}, "
             f"herald ~{rate / 1e3:.0f} kHz, {act['temp_c']} C", fontsize=10)
fig.tight_layout()
out = folder / "full_frame_raw.png"
fig.savefig(out, dpi=130)
np.save(folder / "full_frame_mean.npy", img)
print(f"herald {rate / 1e3:.1f} kHz, {acc:,} acc, max pixel {img.max():.0f} ADU/readout")
print(out)
