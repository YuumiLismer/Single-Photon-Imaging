"""
hold_herald_rate.py -- keep the herald (idler) rate at TRIGGER IN inside a target band by
rotating the pump half-wave plate on the XPS (positioner S.Pos2), with the rate measured by
the PI-MAX4 itself (gain 1, 3 ns gate: light-safe).

    python hold_herald_rate.py --measure                        # rate at the current angle, no move
    python hold_herald_rate.py --calibrate                      # sweep 0..45 deg from the current angle, fit, return
    python hold_herald_rate.py --hold --low 15e3 --high 20e3 --minutes 30
    python hold_herald_rate.py --calibrate --hold --minutes 30  # both

Calibration fits  rate(theta) = B + A cos^2(2 (theta - theta0))  (HWP in front of a polarizer/PBS)
and saves it to <experiment>/figures/waveplate_calibration.json plus a plot. The hold loop
measures the rate, rescales the fitted curve to the measured value (the rate drifts with
coupling and pump power), and moves to the nearest angle that gives the band centre.
A move is capped at --max-step degrees (and any move at 50 deg).

Every single rate reading and every move goes to <experiment>/herald_rate_readings.csv (time, phase,
angle). At the end all such files of the day (or since --history-since) are plotted as rate and angle
vs wall-clock time in <experiment>/figures/herald_rate_vs_time.png (re-plot any time with
agent_scripts/herald_rate_history.py --lock-logs). The camera is closed on exit. The waveplate stays
where the lock left it (use --restore to return it to the start angle).
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

READINGS_CSV = "herald_rate_readings.csv"
READ_FIELDS = ["time", "phase", "angle_deg", "accumulations", "interval_s", "rate_hz", "readout_limited"]
MAX_MOVE_DEG = 50.0                            # safety: largest single waveplate move


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


class ReadingLog:
    """Append-only CSV of every rate reading and every waveplate move (rate empty), with time stamps."""

    def __init__(self, experiment):
        self.path = pathlib.Path(experiment) / READINGS_CSV
        new = not self.path.exists()
        self.f = open(self.path, "a", newline="")
        self.w = csv.DictWriter(self.f, READ_FIELDS)
        if new:
            self.w.writeheader()
        self.phase, self.angle = "", float("nan")

    def reading(self, r):
        self.w.writerow({"time": r["time"], "phase": self.phase, "angle_deg": f"{self.angle:.3f}",
                         "accumulations": r["accumulations"], "interval_s": f"{r['interval_s']:.4f}",
                         "rate_hz": f"{r['rate_hz']:.1f}", "readout_limited": int(r["readout_limited"])})
        self.f.flush()

    def move(self, angle):
        self.angle = angle
        self.w.writerow({**{k: "" for k in READ_FIELDS}, "phase": "move", "angle_deg": f"{angle:.3f}",
                         "time": datetime.datetime.now().isoformat(timespec="milliseconds")})
        self.f.flush()

    def close(self):
        self.f.close()


def move(wp, deg, log=None):
    """Absolute move with the 50 deg single-move limit; logged. No move if already there."""
    now = wp.position()
    if abs(deg - now) > MAX_MOVE_DEG:
        raise ValueError(f"refusing a {deg - now:+.1f} deg move (limit {MAX_MOVE_DEG} deg)")
    ang = wp.move_to(deg) if abs(deg - now) > 1e-3 else now
    if log:
        log.move(ang)
    return ang


_last_rate = [None]                            # carry the accumulation count from one point to the next


def measure(cam, n, log=None, target_s=1.0):
    """Mean and sd of n good (not readout-limited) readings in a row at the current angle."""
    start = int(np.clip(_last_rate[0] * target_s, 1000, 10_000_000)) if _last_rate[0] else 20_000
    args = trm.parse_args(["--average", str(n), "--start-accumulations", str(start), "--target-s", str(target_s)])
    rates = trm.monitor(cam, args, max_readings=4 * n + 6, stop_after_good=n,
                        on_reading=log.reading if log else None)
    if not rates:
        return float("nan"), float("nan")
    _last_rate[0] = float(np.mean(rates))
    return float(np.mean(rates)), float(np.std(rates))


def plot_rate_history(csvs, out_png, low=15e3, high=20e3, title=None):
    """Rate (kHz, log scale) and waveplate angle vs wall-clock time from herald_rate_readings.csv files."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    rows = []
    for c in csvs:
        with open(c, newline="") as f:
            rows += list(csv.DictReader(f))
    rows.sort(key=lambda r: r["time"])
    if not rows:
        return None
    T = lambda rr: [datetime.datetime.fromisoformat(r["time"]) for r in rr]
    rd = [r for r in rows if r["rate_hz"]]
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True, gridspec_kw=dict(height_ratios=[3, 1.4]))
    ax.axhspan(low / 1e3, high / 1e3, color="g", alpha=0.15, label=f"band {low / 1e3:g}-{high / 1e3:g} kHz")
    for i, ph in enumerate(dict.fromkeys(r["phase"] for r in rd)):
        for lim in (False, True):
            S = [r for r in rd if r["phase"] == ph and bool(int(r["readout_limited"])) == lim]
            if S:
                ax.scatter(T(S), [float(r["rate_hz"]) / 1e3 for r in S], s=18 if lim else 12,
                           fc="none" if lim else f"C{i}", ec=f"C{i}",
                           label=(ph or "reading") + (" (readout-limited: lower bound)" if lim else ""))
    ax.set_yscale("log")
    ax.set_ylabel("herald rate at TRIGGER IN (kHz)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7, ncol=2, loc="best")
    ang = [r for r in rows if r["angle_deg"] not in ("", "nan")]
    ax2.step(T(ang), [float(r["angle_deg"]) for r in ang], where="post", color="k", lw=1)
    ax2.plot(T(ang), [float(r["angle_deg"]) for r in ang], ".", color="k", ms=3)
    ax2.set_ylabel("S.Pos2 angle (deg)")
    ax2.grid(alpha=0.3)
    loc = mdates.AutoDateLocator()
    ax2.xaxis.set_major_locator(loc)
    ax2.xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc, show_offset=False))
    ax2.set_xlabel(f"time on {rows[0]['time'][:10]}")
    inb = [low <= float(r["rate_hz"]) <= high for r in rd if r["phase"] == "hold"]
    ax.set_title(title or ("Herald rate vs time (camera-timed, gain 1, 3 ns gate); bottom: pump waveplate angle"
                           + (f"\nlock: {np.mean(inb) * 100:.0f} % of {len(inb)} single readings in band" if inb else "")),
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    return out_png


def session_csvs(since):
    """Every herald_rate_readings.csv under DATA_ROOT modified since `since` (datetime)."""
    return sorted(q for q in ht.DATA_ROOT.glob(f"*/*/{READINGS_CSV}")
                  if datetime.datetime.fromtimestamp(q.stat().st_mtime) >= since)


def calibrate(cam, wp, start, step, span, n, folder, log=None, target_s=1.0, max_rate=250e3):
    """Sweep start .. start+span (span may be negative) in |step| steps. Stops the sweep as soon as a
    point reads above max_rate (the camera is never left at a rate above the ceiling for long)."""
    th, r, s = [], [], []
    if log:
        log.phase = "calibrate"
    step = abs(step) * (1 if span >= 0 else -1)
    for off in np.arange(0, span + 1e-9 * np.sign(step), step):
        ang = move(wp, start + off, log)
        time.sleep(0.5)
        m, sd = measure(cam, n, log, target_s)
        print(f">>> {ang:8.3f} deg: {m / 1e3:7.1f} +/- {sd / 1e3:5.1f} kHz", flush=True)
        th.append(ang); r.append(m); s.append(sd)
        if m > max_rate:
            print(f"rate {m / 1e3:.0f} kHz above the {max_rate / 1e3:.0f} kHz ceiling: sweep stopped", flush=True)
            break
    ok = np.isfinite(r)
    if ok.sum() < 4:
        raise RuntimeError("fewer than 4 calibration points with a rate: laser off?")
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
    ax.axhspan(15, 20, color="g", alpha=0.15)
    ax.set_xlabel(f"{wp.positioner} angle (deg)")
    ax.set_ylabel("trigger rate at TRIGGER IN (kHz)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(folder / "waveplate_calibration.png", dpi=150)
    plt.close(fig)
    return cal


def hold(cam, wp, cal, low, high, minutes, n, max_step, folder, log=None, target_s=1.0):
    target = 0.5 * (low + high)
    path = folder / "herald_rate_lock.csv"
    t_end = time.time() + 60 * minutes
    if log:
        log.phase = "hold"
    with open(path, "a", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time", "angle_deg", "rate_hz", "sd_hz", "in_band", "next_angle_deg"])
        while time.time() < t_end:
            now = wp.position()
            if log:
                log.angle = now
            m, sd = measure(cam, n, log, target_s)
            inside = bool(low <= m <= high)
            nxt = now
            if np.isfinite(m) and not inside:
                nxt = now + float(np.clip(angle_for(target, cal, now, m) - now, -max_step, max_step))
                if abs(nxt - now) < 0.01:          # already at the best angle (band out of reach): no move
                    nxt = now
                    print(f"    band out of reach: {m / 1e3:.1f} kHz is the fitted extreme", flush=True)
            w.writerow([datetime.datetime.now().isoformat(timespec="seconds"), f"{now:.3f}",
                        f"{m:.0f}", f"{sd:.0f}", int(inside), f"{nxt:.3f}"])
            f.flush()
            print(f">>> {datetime.datetime.now():%H:%M:%S}  {now:8.3f} deg  {m / 1e3:7.1f} kHz  "
                  + ("in band" if inside else f"-> {nxt:.3f} deg"), flush=True)
            if nxt != now:
                move(wp, nxt, log)
    print(f"lock log: {path}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--positioner", default="S.Pos2")
    p.add_argument("--measure", action="store_true", help="only measure the rate at the current angle (no move)")
    p.add_argument("--calibrate", action="store_true")
    p.add_argument("--hold", action="store_true")
    p.add_argument("--step", type=float, default=5.0, help="calibration step (deg)")
    p.add_argument("--span", type=float, default=45.0,
                   help="calibration span from the start angle (deg); negative sweeps the other way")
    p.add_argument("--max-rate", type=float, default=250e3,
                   help="calibration stops at a point above this rate (Hz); safety ceiling")
    p.add_argument("--readings", type=int, default=5, help="good rate readings per point")
    p.add_argument("--target-s", type=float, default=1.0, help="seconds per readout (raise if readouts show '>=')")
    p.add_argument("--low", type=float, default=15e3)
    p.add_argument("--high", type=float, default=20e3)
    p.add_argument("--minutes", type=float, default=10.0)
    p.add_argument("--max-step", type=float, default=5.0, help="largest single lock move (deg)")
    p.add_argument("--calibration", default=None, help="JSON from an earlier --calibrate (default: newest)")
    p.add_argument("--restore", action="store_true", help="return to the start angle at the end")
    p.add_argument("--start-angle", type=float, default=None,
                   help="move here first and sweep from here (e.g. to stay below --max-rate); within 50 deg")
    p.add_argument("--park", action="store_true",
                   help="after a calibration-only run, go to the fitted band-centre angle (or the minimum) "
                        "instead of back to the start angle")
    p.add_argument("--history-since", default=None,
                   help="rate-vs-time plot covers readings files modified since this ISO time (default: today 00:00)")
    a = p.parse_args(argv)
    if not (a.measure or a.calibrate or a.hold):
        p.error("give --measure, --calibrate and/or --hold")

    wp = xw.Waveplate(a.positioner)
    start, st = wp.position(), wp.initial_status
    was_disabled = 20 <= st <= 39
    print(f"{a.positioner} at {start:.3f} deg (group status {st})")
    if not (was_disabled or 11 <= st <= 19):
        wp.close()
        raise SystemExit(f"group {wp.group} status {st} is neither disabled nor ready: not touching it")
    experiment = ht.new_experiment("herald_rate_lock")
    folder = experiment / "figures"
    log = ReadingLog(experiment)
    log.angle = start
    cam, _ = ht.open_camera()
    try:
        if a.measure:
            log.phase = "measure"
            m, sd = measure(cam, a.readings, log, a.target_s)
            print(f">>> {start:8.3f} deg: {m / 1e3:7.2f} +/- {sd / 1e3:5.2f} kHz", flush=True)
        sweep_start = start
        if a.start_angle is not None:
            log.phase = "calibrate"
            sweep_start = move(wp, a.start_angle, log)
        cal = (calibrate(cam, wp, sweep_start, a.step, a.span, a.readings, folder, log, a.target_s, a.max_rate)
               if a.calibrate else None)
        if a.hold:
            if cal is None:
                src = a.calibration or max(ht.DATA_ROOT.glob("*/*herald_rate_lock*/figures/waveplate_calibration.json"),
                                           key=lambda q: q.stat().st_mtime)
                cal = json.loads(pathlib.Path(src).read_text())
                print(f"using calibration {src}")
            hold(cam, wp, {k: cal[k] for k in ("A", "B", "theta0")}, a.low, a.high, a.minutes,
                 a.readings, a.max_step, folder, log, a.target_s)
        elif a.calibrate and a.park:
            log.phase = "park"
            now = wp.position()
            goal = angle_for(0.5 * (a.low + a.high), cal, now, model(now, **cal))
            print(f"parking at {goal:.3f} deg (fit: {model(goal, **cal) / 1e3:.1f} kHz)")
            move(wp, goal, log)
        elif a.calibrate and not a.restore:
            log.phase = "return"
            move(wp, start, log)                    # calibration only: put it back
    finally:
        ht.shutdown(cam)
        if a.restore:
            log.phase = "return"
            move(wp, start, log)
        print(f"{a.positioner} now at {wp.position():.3f} deg")
        wp.close(disable=was_disabled and not a.hold)
        log.close()
        since = (datetime.datetime.fromisoformat(a.history_since) if a.history_since
                 else datetime.datetime.combine(datetime.date.today(), datetime.time()))
        png = plot_rate_history(session_csvs(since), folder / "herald_rate_vs_time.png", a.low, a.high)
        print(f"readings: {log.path}\nrate vs time (all readings since {since:%Y-%m-%d %H:%M}): {png}")


if __name__ == "__main__":
    main()
