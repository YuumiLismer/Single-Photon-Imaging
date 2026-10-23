"""Check that the XPS pump waveplate (S.Pos2) controls the herald rate seen at TRIGGER IN.
Steps the waveplate through offsets from its current angle, measures the camera-timed trigger
rate at each (gain 1, 3 ns gate), then returns it to the start angle and leaves the group
enabled/disabled as found. Every reading is logged (with time and angle) to
<experiment>/herald_rate_readings.csv and the day's rate-vs-time plot is redrawn.

    python agent_scripts/waveplate_rate_check.py --offsets 0 5 10 15 22.5 0 --readings 6
    python agent_scripts/waveplate_rate_check.py --offsets 0 2 0      # first small test move
"""
import argparse, datetime, pathlib, sys, time
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht
import hold_herald_rate as hh
import xps_waveplate as xw

p = argparse.ArgumentParser()
p.add_argument("--offsets", type=float, nargs="+", default=[0, 5, 10, 15, 22.5, 0])
p.add_argument("--readings", type=int, default=6, help="good rate readings per angle")
p.add_argument("--target-s", type=float, default=1.0)
p.add_argument("--max-rate", type=float, default=250e3, help="stop stepping (and return) above this rate (Hz)")
a = p.parse_args()
if max(abs(o) for o in a.offsets) > hh.MAX_MOVE_DEG:
    raise SystemExit("offsets must stay within 50 deg of the start angle")

wp = xw.Waveplate()
start, st = wp.position(), wp.initial_status
was_disabled = 20 <= st <= 39
print(f"waveplate {wp.positioner} at {start:.3f} deg, group status {st}")
if not (was_disabled or 11 <= st <= 19):
    wp.close()
    raise SystemExit(f"group status {st} is neither disabled nor ready: not touching it")
experiment = ht.new_experiment("waveplate_rate_check")
log = hh.ReadingLog(experiment)
log.phase, log.angle = "response check", start
cam, _ = ht.open_camera()
rows = []
try:
    for off in a.offsets:
        ang = hh.move(wp, start + off, log)
        time.sleep(1.0)
        m, s = hh.measure(cam, a.readings, log, a.target_s)
        rows.append((off, ang, m, s))
        print(f"\n>>> offset {off:+6.1f} deg (at {ang:.3f}): {m / 1e3:7.2f} +/- {s / 1e3:5.2f} kHz\n", flush=True)
        if m > a.max_rate:
            print(f"rate above the {a.max_rate / 1e3:.0f} kHz ceiling: stopping, returning to start", flush=True)
            break
finally:
    ht.shutdown(cam)
    log.phase = "return"
    end = hh.move(wp, start, log)
    wp.close(disable=was_disabled)
    log.close()
    print(f"waveplate returned to {end:.3f} deg" + (", group disabled again" if was_disabled else ""))
    png = hh.plot_rate_history(hh.session_csvs(datetime.datetime.combine(datetime.date.today(), datetime.time())),
                               experiment / "figures" / "herald_rate_vs_time.png")
    print(f"readings: {log.path}\nrate vs time: {png}")
print("\noffset_deg  angle_deg  rate_kHz  sd_kHz")
for off, ang, m, s in rows:
    print(f"{off:+9.1f}  {ang:9.3f}  {m / 1e3:8.2f}  {s / 1e3:6.2f}")
