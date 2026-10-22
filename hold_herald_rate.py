"""
hold_herald_rate.py -- keep the herald (idler) rate at TRIGGER IN inside a target band by
rotating the pump half-wave plate on the XPS (positioner S.Pos2), with the rate measured by
the PI-MAX4 itself (gain 1, 3 ns gate: light-safe).

    python hold_herald_rate.py --calibrate                      # sweep 0..45 deg from the current angle, fit, return
    python hold_herald_rate.py --hold --low 15e3 --high 20e3 --minutes 30
    python hold_herald_rate.py --calibrate --hold --minutes 30  # both

Calibration fits  rate(theta) = B + A cos^2(2 (theta - theta0))  (HWP in front of a polarizer/PBS)
and saves it to <experiment>/figures/waveplate_calibration.json plus a plot. The hold loop
measures the rate, rescales the fitted curve to the measured value (the rate drifts with
coupling and pump power), and moves to the nearest angle that gives the band centre.
A move is capped at --max-step degrees. Every reading goes to a CSV in the experiment folder.
The camera is closed on exit. The waveplate stays where the lock left it (use --restore to
return it to the start angle).
"""

import argparse
import csv
import datetime
import json
import pathlib
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import hoi_troubleshoot as ht
import trigger_rate_monitor as trm
import xps_waveplate as xw


def model(theta, A, B, theta0):
    return B + A * np.cos(np.deg2rad(2 * (np.asarray(theta) - theta0))) ** 2


def fit(theta, rate):
    """Least-squares fit of model(): grid over theta0 (period 90 deg), linear in A and B."""
    theta, rate = np.asarray(theta, float), np.asarray(rate, float)
    best = None
    for t0 in np.arange(0, 90, 0.05):
        c = np.cos(np.deg2rad(2 * (theta - t0))) ** 2
        M = np.column_stack([c, np.ones_like(c)])
        (A, B), *_ = np.linalg.lstsq(M, rate, rcond=None)
        if A <= 0:
            continue
        r = np.sum((M @ [A, B] - rate) ** 2)
        if best is None or r < best[0]:
            best = (r, A, B, t0)
    if best is None:
        raise ValueError("rate does not vary with angle: is this the pump waveplate?")
    return dict(A=float(best[1]), B=float(best[2]), theta0=float(best[3]))


def angle_for(target, cal, now, measured):
    """Angle nearest `now` where the fitted curve, rescaled to match `measured` at `now`, gives `target`."""
    scale = measured / model(now, **cal)
    c = (target / scale - cal["B"]) / cal["A"]
    if c <= 0 or c >= 1:                       # out of reach: go to the max / min
        x = 0.0 if c >= 1 else 90.0            # 2 (theta - theta0) at the max / min
        cands = [cal["theta0"] + (x + 180 * k) / 2 for k in range(-3, 4)]
    else:
        x = np.rad2deg(np.arccos(np.sqrt(c)))
        cands = [cal["theta0"] + (s * x + 180 * k) / 2 for s in (1, -1) for k in range(-3, 4)]
    return float(min(cands, key=lambda t: abs(t - now)))


def measure(cam, n):
    rates = trm.monitor(cam, trm.parse_args(["--average", str(n)]), max_readings=n + 2)
    return (float(np.mean(rates)), float(np.std(rates))) if rates else (float("nan"), float("nan"))


def calibrate(cam, wp, start, step, span, n, folder):
    th, r, s = [], [], []
    for off in np.arange(0, span + 1e-9, step):
        ang = wp.move_to(start + off)
        time.sleep(0.5)
        m, sd = measure(cam, n)
        print(f">>> {ang:8.3f} deg: {m / 1e3:7.1f} +/- {sd / 1e3:5.1f} kHz", flush=True)
        th.append(ang); r.append(m); s.append(sd)
    ok = np.isfinite(r)
    cal = fit(np.array(th)[ok], np.array(r)[ok])
    out = {**cal, "positioner": wp.positioner, "theta": th, "rate": r, "sd": s,
           "time": datetime.datetime.now().isoformat(timespec="seconds")}
    print(f"fit: max {(cal['A'] + cal['B']) / 1e3:.1f} kHz at {cal['theta0'] % 90:.2f} deg (mod 90), "
          f"min {cal['B'] / 1e3:.1f} kHz at {(cal['theta0'] + 45) % 90:.2f} deg (mod 90)")
    (folder / "waveplate_calibration.json").write_text(json.dumps(out, indent=1))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.errorbar(th, np.array(r) / 1e3, np.array(s) / 1e3, fmt="o", label="measured")
    x = np.linspace(min(th), max(th), 300)
    ax.plot(x, model(x, **cal) / 1e3, label=r"B + A cos$^2$(2($\theta-\theta_0$))")
    ax.set_xlabel(f"{wp.positioner} angle (deg)")
    ax.set_ylabel("trigger rate at TRIGGER IN (kHz)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(folder / "waveplate_calibration.png", dpi=150)
    return cal


def hold(cam, wp, cal, low, high, minutes, n, max_step, folder):
    target = 0.5 * (low + high)
    path = folder / "herald_rate_lock.csv"
    t_end = time.time() + 60 * minutes
    with open(path, "a", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time", "angle_deg", "rate_hz", "sd_hz", "in_band", "next_angle_deg"])
        while time.time() < t_end:
            now = wp.position()
            m, sd = measure(cam, n)
            inside = bool(low <= m <= high)
            nxt = now
            if np.isfinite(m) and not inside:
                nxt = now + float(np.clip(angle_for(target, cal, now, m) - now, -max_step, max_step))
            w.writerow([datetime.datetime.now().isoformat(timespec="seconds"), f"{now:.3f}",
                        f"{m:.0f}", f"{sd:.0f}", int(inside), f"{nxt:.3f}"])
            f.flush()
            print(f">>> {datetime.datetime.now():%H:%M:%S}  {now:8.3f} deg  {m / 1e3:7.1f} kHz  "
                  + ("in band" if inside else f"-> {nxt:.3f} deg"), flush=True)
            if nxt != now:
                wp.move_to(nxt)
    print(f"lock log: {path}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--positioner", default="S.Pos2")
    p.add_argument("--calibrate", action="store_true")
    p.add_argument("--hold", action="store_true")
    p.add_argument("--step", type=float, default=5.0, help="calibration step (deg)")
    p.add_argument("--span", type=float, default=45.0, help="calibration span from the start angle (deg)")
    p.add_argument("--readings", type=int, default=5, help="rate readings per point")
    p.add_argument("--low", type=float, default=15e3)
    p.add_argument("--high", type=float, default=20e3)
    p.add_argument("--minutes", type=float, default=10.0)
    p.add_argument("--max-step", type=float, default=5.0, help="largest single lock move (deg)")
    p.add_argument("--calibration", default=None, help="JSON from an earlier --calibrate (default: newest)")
    p.add_argument("--restore", action="store_true", help="return to the start angle at the end")
    a = p.parse_args(argv)
    if not (a.calibrate or a.hold):
        p.error("give --calibrate and/or --hold")

    folder = ht.new_experiment("herald_rate_lock") / "figures"
    folder.mkdir(parents=True, exist_ok=True)
    wp = xw.Waveplate(a.positioner)
    start, was_disabled = wp.position(), 20 <= wp.initial_status <= 39
    print(f"{a.positioner} at {start:.3f} deg (group status {wp.initial_status})")
    cam, _ = ht.open_camera()
    try:
        cal = calibrate(cam, wp, start, a.step, a.span, a.readings, folder) if a.calibrate else None
        if a.hold:
            if cal is None:
                src = a.calibration or max(ht.DATA_ROOT.glob("*/*herald_rate_lock*/figures/waveplate_calibration.json"),
                                           key=lambda q: q.stat().st_mtime)
                cal = json.loads(pathlib.Path(src).read_text())
                print(f"using calibration {src}")
            hold(cam, wp, {k: cal[k] for k in ("A", "B", "theta0")}, a.low, a.high, a.minutes,
                 a.readings, a.max_step, folder)
        elif not a.restore:
            wp.move_to(start)                       # calibration only: put it back
    finally:
        ht.shutdown(cam)
        if a.restore:
            wp.move_to(start)
        print(f"{a.positioner} now at {wp.position():.3f} deg")
        wp.close(disable=was_disabled and not a.hold)


if __name__ == "__main__":
    main()
