"""Check that the XPS pump waveplate (S.Pos2) controls the herald rate seen at TRIGGER IN.
Steps the waveplate through offsets from its current angle, measures the camera-timed trigger
rate at each (gain 1, 3 ns gate), then returns it to the start angle and leaves the group
enabled/disabled as found.

    python agent_scripts/waveplate_rate_check.py --offsets 0 5 10 15 22.5 0 --readings 6
"""
import argparse, pathlib, sys, time
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht
import trigger_rate_monitor as trm
import xps_waveplate as xw

p = argparse.ArgumentParser()
p.add_argument("--offsets", type=float, nargs="+", default=[0, 5, 10, 15, 22.5, 0])
p.add_argument("--readings", type=int, default=6, help="rate readings per angle")
a = p.parse_args()

wp = xw.Waveplate()
start, was_disabled = wp.position(), 20 <= wp.initial_status <= 39
print(f"waveplate {wp.positioner} at {start:.3f} deg, group status {wp.initial_status}")
cam, _ = ht.open_camera()
rows = []
try:
    for off in a.offsets:
        ang = wp.move_to(start + off)
        time.sleep(1.0)
        rates = trm.monitor(cam, trm.parse_args(["--average", str(a.readings)]), max_readings=a.readings + 2)
        m = np.mean(rates) if rates else float("nan")
        s = np.std(rates) if rates else float("nan")
        rows.append((off, ang, m, s))
        print(f"\n>>> offset {off:+6.1f} deg (at {ang:.3f}): {m / 1e3:7.1f} +/- {s / 1e3:5.1f} kHz\n", flush=True)
finally:
    ht.shutdown(cam)
    end = wp.move_to(start)
    wp.close(disable=was_disabled)
    print(f"waveplate returned to {end:.3f} deg" + (", group disabled again" if was_disabled else ""))
print("\noffset_deg  angle_deg  rate_kHz  sd_kHz")
for off, ang, m, s in rows:
    print(f"{off:+9.1f}  {ang:9.3f}  {m / 1e3:8.1f}  {s / 1e3:6.1f}")
