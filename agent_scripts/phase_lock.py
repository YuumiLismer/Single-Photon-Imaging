"""Two-path interference phase: calibrate the phase-shifter slides, record the drift, and hold the phase.

Phase readout: a reference picture at t = 0 sets the fringe direction and spacing. Every frame is summed along the
stripes (one value per row of the fringe pattern) into a profile across the fringes; the fringe part of that
profile (profile minus its smooth envelope) is compared with the t = 0 profile: phase = angle between their
fringe components, r = correlation (+1 same picture, -1 a pi shift; the scatter of profile(t) against
profile(0) is a straight line of slope +1 when the phase is back).

    # rad per degree of each slide (moves it, returns it to the start)
    python agent_scripts/phase_lock.py --calibrate P.Pos1 P.Pos2
    # record for 3 h, correcting every 5 min with the first slide (second slide takes over if the first runs out)
    python agent_scripts/phase_lock.py --minutes 180 --check-every 300 --lock P.Pos1 P.Pos2
    # free-running drift, no corrections
    python agent_scripts/phase_lock.py --minutes 60 --free

Stop early by creating <experiment>/STOP. Live window: zoomed picture, phase vs time, profile(t) vs profile(0), r vs time.
"""
import argparse
import json
import pathlib
import sys
import time

import numpy as np
import matplotlib
from scipy.ndimage import gaussian_filter, gaussian_filter1d

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht
from xps_waveplate import XPS, Waveplate

SAT_ADU = 60_000
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--gate", type=float, default=10_000); p.add_argument("--acc", type=int, default=4096)
p.add_argument("--gain", type=int, default=1)
p.add_argument("--nref", type=int, default=10, help="frames averaged for the t = 0 reference")
p.add_argument("--calibrate", nargs="+", default=None, metavar="POSITIONER")
p.add_argument("--minutes", type=float, default=0)
p.add_argument("--check-every", type=float, default=300, help="s between corrections")
p.add_argument("--lock", nargs="+", default=[], metavar="POSITIONER", help="slides to correct with (first is primary)")
p.add_argument("--cal", default=None, help="calibration json (default: the newest phase_cal.json)")
p.add_argument("--free", action="store_true", help="record only, no corrections")
p.add_argument("--tol", type=float, default=0.15, help="rad: phase error accepted at a check")
p.add_argument("--max-step", type=float, default=1.0, help="deg: largest single slide move")
p.add_argument("--range", type=float, default=3.0, help="deg: a slide may move this far from its start")
p.add_argument("--backlash", type=float, default=0.2, help="deg: every correction approaches its target from below by this much")
p.add_argument("--navg", type=int, default=4, help="frames averaged per phase measurement at a check")
p.add_argument("--no-window", action="store_true")
a = p.parse_args()
if not a.no_window:
    matplotlib.use("TkAgg")
import matplotlib.pyplot as plt   # noqa: E402
import camera_control as cc       # noqa: E402

wrap = lambda x: (x + np.pi) % (2 * np.pi) - np.pi


class Fringes:
    """Reference picture -> fringe direction/period, beam mask, and the per-frame profile/phase/correlation."""

    def __init__(self, ref):
        env = gaussian_filter(ref, 12)
        self.mask = env > 0.3 * env.max()
        res = (gaussian_filter(ref, 1) - env) * self.mask
        F = np.abs(np.fft.fft2(res)) ** 2
        ky, kx = np.meshgrid(np.fft.fftfreq(ref.shape[0]), np.fft.fftfreq(ref.shape[1]), indexing="ij")
        k = np.hypot(kx, ky)
        F[(k < 1 / 120) | (k > 1 / 4)] = 0
        iy, ix = np.unravel_index(F.argmax(), F.shape)
        self.period = 1 / k[iy, ix]
        ux, uy = kx[iy, ix] / k[iy, ix], ky[iy, ix] / k[iy, ix]
        yy, xx = np.nonzero(self.mask)
        u = xx * ux + yy * uy                                  # coordinate across the stripes (px)
        self.u0 = np.floor(u.min())
        self.bins = (u - self.u0).astype(int)
        self.idx = (yy, xx)
        self.nb = self.bins.max() + 1
        self.ucoord = self.u0 + np.arange(self.nb)
        self.cnt = np.bincount(self.bins, minlength=self.nb)
        ys, xs = yy, xx
        self.box = (max(0, ys.min() - 30), ys.max() + 30, max(0, xs.min() - 30), xs.max() + 30)
        self.angle = float(np.degrees(np.arctan2(uy, ux)))
        self.p0, self.c0 = None, None
        self.p0, self.c0, _, _, _ = self.measure(ref)

    def profile(self, img):
        prof = np.bincount(self.bins, img[self.idx], minlength=self.nb)      # sum along each stripe row
        ok = self.cnt >= 3
        smooth = gaussian_filter1d(prof, self.period * 0.8)
        return np.where(ok, prof - smooth, 0.0), np.where(ok, smooth, 0.0)

    def measure(self, img):
        fr, env = self.profile(img)
        c = np.sum(fr * np.exp(-2j * np.pi * self.ucoord / self.period))
        vis = float(2 * np.abs(c) / max(env.sum(), 1e-9))
        if self.c0 is None:
            return fr, c, 0.0, 1.0, vis
        dphi = float(np.angle(c * np.conj(self.c0)))
        r = float(np.corrcoef(fr, self.p0)[0, 1])
        return fr, c, dphi, r, vis


cam, limits = ht.open_camera()
base = {**ht.DEFAULT_SETTINGS, "label": "phase", "gain": a.gain, "trigger_source": "Internal",
        "internal_trigger_hz": 10_000, "gate_delay_ns": 100.0, "photon_mode": "Disabled",
        "gate_width_ns": a.gate, "accumulations": a.acc}


def start():
    try:
        cam.stop_acquisition()
    except Exception:
        pass
    ht.configure(cam, base, verbose=False)
    cam.start_acquisition()


def frame(fresh_after=None):
    """Next frame (median-subtracted). fresh_after: discard frames that may have started before this time."""
    while True:
        try:
            if not cc.poll_frame(cam, 0.5):
                if not a.no_window:
                    plt.pause(0.01)
                continue
            fr = cam.read_multiple_images()
        except Exception as e:                                 # acquisition dropped: restart it
            print(f"{time.strftime('%H:%M:%S')} camera error, restarting: {e}", flush=True)
            start(); continue
        if not fr:
            continue
        t = time.time()
        if fresh_after is not None and t - fresh_after < 1.2 * a.acc / 10_000 + 0.1:
            continue
        img = np.asarray(fr[-1], dtype=float)
        if (img >= SAT_ADU).sum() > 25:
            print("WARNING: saturating pixels (reduce the light or --acc)", flush=True)
        return t, img - np.median(img)


def mean_phase(n, fresh_after=None):
    cs = []
    for _ in range(n):
        t, img = frame(fresh_after); fresh_after = None
        cs.append(np.exp(1j * fr_.measure(img)[2]))
    return float(np.angle(np.mean(cs))), img


xps = XPS() if (a.calibrate or not a.free) else None
exp = ht.new_experiment("agent_phase_cal" if a.calibrate else "agent_phase")
out = exp / "figures"; out.mkdir(parents=True, exist_ok=True)
start()
print("taking the t = 0 reference ...", flush=True)
ref = np.mean([frame()[1] for _ in range(a.nref)], axis=0)
fr_ = Fringes(ref)
np.save(out / "reference.npy", ref)
print(f"fringes: period {fr_.period:.1f} px, direction {fr_.angle:.0f} deg, {fr_.mask.sum()} px in the beam, "
      f"visibility of the summed profile {fr_.measure(ref)[4]:.2f}", flush=True)

try:
    if a.calibrate:
        cal = {}
        for pos in a.calibrate:
            wp = Waveplate(pos, xps)
            th0 = wp.position()
            wp.move_to(th0 - 0.2); wp.move_to(th0)            # always approach from below (backlash)
            ph0, _ = mean_phase(3, time.time())
            angs, phs, step, unw = [th0], [0.0], 0.02, 0.0
            while abs(unw) < 2.5 * np.pi and abs(angs[-1] - th0) < a.range and len(angs) < 80:
                th = wp.move_to(angs[-1] + step)
                ph, _ = mean_phase(3, time.time())
                phs.append(phs[-1] + wrap((ph - ph0) - phs[-1])); angs.append(th)   # unwrapped along the scan
                print(f"  {pos} {th:9.4f} deg  phase {phs[-1]:+7.3f} rad", flush=True)
                unw = phs[-1]
                if len(angs) >= 4 and abs(phs[-1] - phs[-4]) < 0.15:   # too little per step: bigger steps
                    step = min(step * 2, 0.5)
            wp.move_to(th0 - 0.2); wp.move_to(th0)
            angs, phs = np.array(angs), np.array(phs)
            slope = float(np.polyfit(angs - th0, phs, 1)[0])
            cal[pos] = dict(rad_per_deg=slope, deg_per_pi=float(np.pi / abs(slope)) if slope else None,
                            start_deg=th0, angles=angs.tolist(), phases=phs.tolist())
            print(f"{pos}: {slope:+.3f} rad/deg -> pi needs {np.pi / abs(slope):.4f} deg" if slope else f"{pos}: no phase change",
                  flush=True)
        (out / "phase_cal.json").write_text(json.dumps(cal, indent=1))
        fig, ax = plt.subplots(1, len(cal), figsize=(6 * len(cal), 4.5), squeeze=False)
        for axi, (pos, c) in zip(ax[0], cal.items()):
            axi.plot(np.array(c["angles"]) - c["start_deg"], c["phases"], "o-")
            axi.axhline(np.pi, color="gray", ls=":"); axi.axhline(-np.pi, color="gray", ls=":")
            axi.set_xlabel(f"{pos} angle - {c['start_deg']:.3f} deg"); axi.set_ylabel("fringe phase (rad)")
            axi.set_title(f"{pos}: {c['rad_per_deg']:+.2f} rad/deg, pi = {c['deg_per_pi']:.3f} deg" if c["deg_per_pi"] else pos)
        fig.tight_layout(); fig.savefig(out / "phase_cal.png", dpi=100)
        print(out / "phase_cal.json", flush=True)
    else:
        slopes, wps, th_start = {}, {}, {}
        if not a.free:
            calf = pathlib.Path(a.cal) if a.cal else max(ht.DATA_ROOT.glob("*/*_agent_phase_cal/figures/phase_cal.json"),
                                                         key=lambda q: q.stat().st_mtime)
            cal = json.loads(pathlib.Path(calf).read_text()); print("calibration:", calf, flush=True)
            for pos in a.lock:
                slopes[pos] = cal[pos]["rad_per_deg"]; wps[pos] = Waveplate(pos, xps); th_start[pos] = wps[pos].position()
        log = open(out / "phase_log.csv", "w")
        # drift_rad: the phase drift itself, unwrapped frame to frame, with the jumps made by our own corrections left out
        log.write("t_s,phase_rad,drift_rad,r,visibility,correction_rad," + ",".join(f"{q}_deg" for q in a.lock) + "\n")
        t0 = time.time(); last_check = t0; prev = 0.0; unw = 0.0; applied = 0.0
        T, PH, R, V, UNW, CHK = [], [], [], [], [], []
        if not a.no_window:
            plt.ion()
            fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))
            y0, y1, x0, x1 = fr_.box
            imh = ax[0, 0].imshow(ref[y0:y1, x0:x1], cmap="inferno", extent=[x0, x1, y1, y0])
            ax[0, 0].set_title(f"live (period {fr_.period:.1f} px)")
            lph, = ax[0, 1].plot([], [], ".", ms=3, label="measured phase")
            lfree, = ax[0, 1].plot([], [], "-", color="gray", lw=1, label="drift (unwrapped, corrections removed)")
            ax[0, 1].set_xlabel("time (min)"); ax[0, 1].set_ylabel("phase vs t = 0 (rad)"); ax[0, 1].legend(loc="upper left")
            lsc, = ax[1, 0].plot([], [], ".", ms=4); ax[1, 0].set_xlabel("profile at t = 0 (summed along the stripes)")
            ax[1, 0].set_ylabel("profile now"); m = np.abs(fr_.p0).max() * 1.3
            ax[1, 0].plot([-m, m], [-m, m], "k:", lw=0.8); ax[1, 0].set_xlim(-m, m); ax[1, 0].set_ylim(-m, m)
            lr, = ax[1, 1].plot([], [], ".", ms=3); ax[1, 1].set_ylim(-1.05, 1.05)
            ax[1, 1].set_xlabel("time (min)"); ax[1, 1].set_ylabel("correlation r with t = 0")
            fig.tight_layout()
        stop = exp / "STOP"
        while time.time() - t0 < a.minutes * 60 and not stop.exists():
            t, img = frame()
            prof, _, dphi, r, vis = fr_.measure(img)
            unw += wrap(dphi - prev); prev = dphi
            T.append((t - t0) / 60); PH.append(dphi); UNW.append(unw); R.append(r); V.append(vis)
            angs = [wps[q].position() for q in a.lock] if wps else []
            log.write(f"{t - t0:.2f},{dphi:.4f},{unw:.4f},{r:.4f},{vis:.4f},{applied:.4f}," + ",".join(f"{v:.4f}" for v in angs) + "\n")
            log.flush()
            if not a.free and t - last_check >= a.check_every:
                last_check = t
                for it in range(6):                                 # correct until back at the t = 0 picture
                    err, img = mean_phase(a.navg)
                    if it == 0:                                     # drift up to the moment we start moving
                        unw += wrap(err - prev)
                    if abs(err) < a.tol:
                        break
                    for pos in a.lock:                              # first slide with travel left
                        d = float(np.clip(-err / slopes[pos], -a.max_step, a.max_step))
                        tgt = wps[pos].position() + d
                        if abs(tgt - th_start[pos]) <= a.range:
                            if d < 0:                               # same approach direction as the calibration
                                wps[pos].move_to(tgt - a.backlash)
                            wps[pos].move_to(tgt); applied += d * slopes[pos]
                            break
                    else:
                        print("all slides at the end of their range", flush=True); break
                    frame(time.time())
                else:
                    err, img = mean_phase(a.navg)
                prev = err                                          # the correction jump is not drift
                print(f"{time.strftime('%H:%M:%S')} check: error after {it} moves {err:+.3f} rad, r {fr_.measure(img)[3]:+.3f}, "
                      f"slides {[round(wps[q].position(), 4) for q in a.lock]}", flush=True)
                CHK.append(T[-1])
                f2, a2 = plt.subplots(figsize=(5, 5))
                a2.plot(fr_.p0, fr_.measure(img)[0], ".", ms=4); a2.set_xlabel("profile at t = 0"); a2.set_ylabel("profile now")
                a2.set_title(f"t = {T[-1]:.1f} min after correction, r = {fr_.measure(img)[3]:+.3f}")
                f2.tight_layout(); f2.savefig(out / f"scatter_{T[-1]:06.1f}min.png", dpi=90); plt.close(f2)
            if not a.no_window:
                y0, y1, x0, x1 = fr_.box
                imh.set_data(img[y0:y1, x0:x1]); imh.set_clim(*np.percentile(img[y0:y1, x0:x1], [1, 99.9]))
                lph.set_data(T, PH); lfree.set_data(T, UNW)
                ax[0, 1].relim(); ax[0, 1].autoscale_view()
                lsc.set_data(fr_.p0, prof); lr.set_data(T, R); ax[1, 1].set_xlim(0, max(T[-1], 1))
                ax[0, 0].set_title(f"live: r = {r:+.2f}, phase {dphi:+.2f} rad, visibility {vis:.2f}")
                plt.pause(0.01)
        log.close()
        d = np.genfromtxt(out / "phase_log.csv", delimiter=",", names=True)
        f3, ax3 = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
        tm = d["t_s"] / 60
        ax3[0].plot(tm, d["drift_rad"], "-", color="gray", lw=1, label="drift (unwrapped; corrections removed)")
        ax3[0].plot(tm, d["phase_rad"], ".", ms=2, label="measured phase vs t = 0 (wrapped)")
        for c in CHK:
            ax3[0].axvline(c, color="C2", lw=0.5, alpha=0.5)
        ax3[0].set_ylabel("phase (rad)"); ax3[0].legend()
        ax3[0].set_title(("FREE RUNNING (no corrections)" if a.free else f"LOCKED: corrections every {a.check_every:.0f} s (green)")
                         + f"; drift rms {np.std(d['drift_rad']):.2f} rad, range {np.ptp(d['drift_rad']):.2f} rad")
        ax3[1].plot(tm, d["r"], ".", ms=2); ax3[1].set_ylabel("r with t = 0 picture"); ax3[1].set_ylim(-1.05, 1.05)
        ax3[2].plot(tm, d["visibility"], ".", ms=2); ax3[2].set_ylabel("fringe visibility"); ax3[2].set_xlabel("time (min)")
        f3.tight_layout(); f3.savefig(out / "phase_vs_time.png", dpi=100)
        print(out / "phase_vs_time.png", flush=True)
finally:
    try:
        cam.stop_acquisition()
    except Exception:
        pass
    ht.shutdown(cam)
    if xps:
        xps.close()
