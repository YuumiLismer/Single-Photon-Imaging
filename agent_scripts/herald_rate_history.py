"""Herald (external trigger) rate vs wall-clock time from the readout time stamps of every saved
External-trigger run (agent session 5; Tristan: "heralds as a function of time to see the decay").

    python agent_scripts/herald_rate_history.py [--since 2026-10-21] [--out <folder>]

Per run: gates/s G = accumulations / median readout interval. A gate blocks the trigger input while
it is pending/open (dead time tau = gate delay + width), so the herald rate is R = G / (1 - G*tau)
(checked 21 Oct 14:00: 68 kHz with a 3 ns gate gave 40 k gates/s with a 10 us gate). The correction
over-corrects for late/long gates (5-10 us: the input re-arms before the gate ends), so only runs with
tau <= 1.2 us (correction < ~15 %) are plotted as the herald rate; the rest are shown as gates/s. Runs whose
interval is within 1.3x the 0.226 s readout floor are readout-limited (rate is only a lower bound)
and are drawn hollow. Costs no camera time. Writes herald_rate_vs_time.png + .json.

    python agent_scripts/herald_rate_history.py --lock-logs [--since 2026-10-23T11:00] [--out <folder>]
        every single reading logged by hold_herald_rate.py / waveplate_rate_check.py
        (herald_rate_readings.csv files modified since --since): rate vs time with the 15-20 kHz band
        (--low/--high) and the pump waveplate angle in a second panel -> herald_rate_vs_time_lock.png
"""
import argparse, sys, pathlib, json, datetime
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht

READOUT_FLOOR_S = 0.226
p = argparse.ArgumentParser()
p.add_argument("--since", default="2026-10-20"); p.add_argument("--out", default=None)
p.add_argument("--lock-logs", action="store_true", help="plot the waveplate/lock reading logs instead of saved runs")
p.add_argument("--low", type=float, default=15e3); p.add_argument("--high", type=float, default=20e3)
a = p.parse_args()
if a.lock_logs:
    import hold_herald_rate as hh
    csvs = hh.session_csvs(datetime.datetime.fromisoformat(a.since))
    out = pathlib.Path(a.out) if a.out else (ht.new_experiment("agent_herald_rate", quiet=True) / "figures")
    out.mkdir(parents=True, exist_ok=True)
    png = hh.plot_rate_history(csvs, out / "herald_rate_vs_time_lock.png", a.low, a.high)
    print(len(csvs), "reading logs ->", png)
    sys.exit()
rows = []
for f in sorted(ht.DATA_ROOT.rglob("*.npz")):
    try:
        with np.load(f, allow_pickle=False) as z:
            meta = json.loads(str(z["meta"])); t = np.asarray(z["t"], float)
    except Exception:
        continue
    act = meta.get("actual", {})
    if str(act.get("trigger_source", "")).lower() != "external" or len(t) < 3:
        continue
    start = datetime.datetime.fromisoformat(meta["started"])
    if start < datetime.datetime.fromisoformat(a.since):
        continue
    dt = float(np.median(np.diff(t))); acc = int(act["accumulations"])
    tau = (float(act["gate_delay_ns"]) + float(act["gate_width_ns"])) * 1e-9
    G = acc / dt; R = G / (1 - G * tau) if G * tau < 0.95 else np.nan
    rows.append(dict(time=(start + datetime.timedelta(seconds=float(t[len(t) // 2] - t[0]))).isoformat(timespec="seconds"),
                     gates_per_s=G, herald_hz=R, dt=dt, tau_ns=tau * 1e9, acc=acc, limited=bool(dt < 1.3 * READOUT_FLOOR_S),
                     label=meta["settings"]["label"]))
out = pathlib.Path(a.out) if a.out else (ht.new_experiment("agent_herald_rate", quiet=True) / "figures")
out.mkdir(parents=True, exist_ok=True)
(out / "herald_rate_vs_time.json").write_text(json.dumps(rows, indent=1))
days = sorted({r["time"][:10] for r in rows})
fig, axes = plt.subplots(1, len(days), figsize=(7.5 * len(days), 4.8), squeeze=False)
for ax, day in zip(axes[0], days):
    R = [r for r in rows if r["time"][:10] == day]
    for lim, kw in ((False, dict(fc="C0", ec="C0", label="herald rate (dead-time corrected)")),
                    (True, dict(fc="none", ec="C0", label="readout-limited (lower bound)"))):
        S = [r for r in R if r["limited"] == lim and np.isfinite(r["herald_hz"]) and r["tau_ns"] <= 1200]
        ax.scatter([datetime.datetime.fromisoformat(r["time"]) for r in S], [r["herald_hz"] / 1e3 for r in S], s=14, **kw)
    S = [r for r in R if not r["limited"] and r["tau_ns"] > 1200]
    ax.scatter([datetime.datetime.fromisoformat(r["time"]) for r in S], [r["gates_per_s"] / 1e3 for r in S], s=6,
               c="0.6", label="gates/s, late/long gates (dead-time limited)")
    ax.axhline(15, color="g", ls="--", lw=1, label="idler target 15 kHz"); ax.axhline(20, color="r", ls=":", lw=1, label="20 kHz (trigger errors above)")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M")); ax.set_xlabel(f"time on {day}")
    ax.set_ylabel("kHz"); ax.grid(alpha=.3); ax.set_ylim(bottom=0)
ax.legend(fontsize=8)
fig.suptitle("Herald (External trigger, 0.5 V rising, 50 ohm) rate vs time, from camera readout time stamps", fontsize=11)
fig.tight_layout(); fig.savefig(out / "herald_rate_vs_time.png", dpi=120)
print(len(rows), "external runs ->", out / "herald_rate_vs_time.png")
