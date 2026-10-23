"""
trigger_rate_monitor.py -- live read-out of how many trigger pulses per second reach the
PI-MAX4's TRIGGER IN, using only the camera (no counter software).

    python trigger_rate_monitor.py                      # 0.5 V threshold, rising edge
    python trigger_rate_monitor.py --threshold 0.3      # try another threshold
    python trigger_rate_monitor.py --csv rates.csv      # also log every reading
Stop with Ctrl-C. The camera is closed (lowest gain, intensifier off) on exit.

How it works: in "Readout Per Trigger" mode with N on-CCD accumulations, the camera fires
one gate per trigger and reads out after N triggers. So the time between readouts is
N / (trigger rate), and rate = N / interval. N is adjusted automatically so a readout takes
about --target-s seconds; if N were too small, the CCD readout time (~0.23 s) would set
the interval instead and the number would only be a lower bound (flagged with ">=").

Light safety: emICCD gain is set to its minimum and the gate to 3 ns, so the intensifier
is open for a negligible fraction of the time. Close the notebook's camera first (only one
program can hold the camera).
"""

import argparse
import collections
import csv
import datetime
import sys
import time
import pathlib

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import hoi_troubleshoot as ht


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Live trigger-rate monitor for the PI-MAX4 TRIGGER IN.")
    p.add_argument("--threshold", type=float, default=0.5, help="Trigger Threshold (V)")
    p.add_argument("--edge", default="Rising Edge", choices=["Rising Edge", "Falling Edge"])
    p.add_argument("--target-s", type=float, default=1.0, help="aim for this many seconds per readout")
    p.add_argument("--start-accumulations", type=int, default=20_000)
    p.add_argument("--average", type=int, default=10, help="readouts in the rolling average")
    p.add_argument("--timeout-s", type=float, default=10.0, help="report 'no triggers' after this long")
    p.add_argument("--csv", default=None, help="optional CSV file to log every reading")
    return p.parse_args(argv)


def _configure(cam, args, accumulations, quiet=False):
    s = {**ht.DEFAULT_SETTINGS, "label": "trigger_rate_monitor", "gain": 1, "gate_delay_ns": 100.0,
         "gate_width_ns": 3.0, "accumulations": int(accumulations), "trigger_source": "External",
         "trigger_threshold_v": args.threshold, "trigger_edge": args.edge, "photon_mode": "Disabled"}
    ht.configure(cam, s, verbose=not quiet)
    return int(cam.get_attribute_value("Accumulations"))


def monitor(cam=None, args=None, max_readings=None, clock=time.time, on_reading=None, stop_after_good=None):
    """Run the monitor. Returns the last --average good (not readout-limited) rates in Hz.
    on_reading(dict(time, accumulations, interval_s, rate_hz, readout_limited)) is called for every
    reading; stop_after_good=N returns as soon as N good readings in a row (same accumulations) are in."""
    a = args or parse_args([])
    opened_here = cam is None
    if cam is None:
        cam, _ = ht.open_camera()
    import camera_control as cc
    attrs = cam.get_all_attributes()
    acc_max = getattr(attrs.get("Accumulations"), "max", None) or 10_000_000
    ro = cam.get_attribute_value("Readout Time Calculation") if "Readout Time Calculation" in attrs else None
    ro_s = float(ro) / 1000 if isinstance(ro, (int, float)) and ro > 0 else 0.23

    n_acc = _configure(cam, a, min(a.start_accumulations, acc_max))
    writer = None
    if a.csv:
        f = open(a.csv, "a", newline="")
        writer = csv.writer(f)
        writer.writerow(["time", "accumulations", "interval_s", "rate_hz", "readout_limited"])
    recent = collections.deque(maxlen=a.average)
    readings = 0
    print(f"\nTrigger threshold {a.threshold} V, {a.edge}. Readout time {ro_s:.3f} s. Ctrl-C to stop.\n")
    print(f"{'time':>8}  {'kHz':>9}  {'avg kHz':>16}  accumulations")
    def done():
        return ((max_readings is not None and readings >= max_readings)
                or (stop_after_good is not None and len(recent) >= stop_after_good))

    try:
        while not done():
            cam.start_acquisition()
            last = None
            changed = False
            try:
                while not changed and not done():
                    t_wait = clock()
                    got = False
                    while clock() - t_wait < a.timeout_s:
                        if cc.poll_frame(cam, 0.5):
                            got = True
                            break
                    if not got:
                        print(f"{datetime.datetime.now():%H:%M:%S}  no triggers in {a.timeout_s:.0f} s "
                              f"(cable, APD, threshold {a.threshold} V, edge?)")
                        recent.clear()
                        last = None
                        if n_acc > 1000:                  # maybe just very slow: ask for fewer
                            n_acc = _configure_after_stop(cam, a, max(n_acc // 10, 1000))
                            changed = True
                        continue
                    imgs = cam.read_multiple_images()
                    if len(imgs) != 1:            # fell behind: arrival times unknown, re-sync
                        last = clock()
                        continue
                    for _ in imgs:
                        now = clock()
                        if last is not None:
                            dt = now - last
                            limited = dt < 1.5 * ro_s
                            rate = n_acc / dt
                            readings += 1
                            if not limited:
                                recent.append(rate)
                            avg = (f"{np.mean(recent) / 1e3:7.2f} +/- {np.std(recent) / 1e3:5.2f}"
                                   if len(recent) > 1 else "")
                            print(f"{datetime.datetime.now():%H:%M:%S}  {'>=' if limited else '  '}"
                                  f"{rate / 1e3:7.2f}  {avg:>16}  {n_acc:,}", flush=True)
                            if on_reading:
                                on_reading(dict(time=datetime.datetime.now().isoformat(timespec="milliseconds"),
                                                accumulations=n_acc, interval_s=dt, rate_hz=rate,
                                                readout_limited=bool(limited)))
                            if writer:
                                writer.writerow([datetime.datetime.now().isoformat(timespec="seconds"),
                                                 n_acc, f"{dt:.4f}", f"{rate:.1f}", int(limited)])
                            # keep each readout near target_s so readout time never limits it
                            want = int(np.clip(rate * a.target_s, 1000, acc_max))
                            if limited or not (0.5 * n_acc <= want <= 2 * n_acc):
                                if limited:
                                    want = int(min(n_acc * 4, acc_max))
                                if want != n_acc:
                                    n_acc = _configure_after_stop(cam, a, want)
                                    recent.clear()
                                    changed = True
                                    break
                        last = now
            finally:
                try:
                    cam.stop_acquisition()
                except Exception:
                    pass
    except KeyboardInterrupt:
        print("\nstopped")
    finally:
        if writer:
            f.close()
        if opened_here:
            ht.shutdown(cam)
    return list(recent)


def _configure_after_stop(cam, a, n):
    try:
        cam.stop_acquisition()
    except Exception:
        pass
    time.sleep(0.2)
    n_set = _configure(cam, a, n, quiet=True)
    print(f"           (accumulations -> {n_set:,})")
    return n_set


if __name__ == "__main__":
    monitor(args=parse_args())
