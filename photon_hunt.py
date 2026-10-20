"""
photon_hunt.py -- unattended search for heralded single photons on the PI-MAX4, then a
long 3 ns image. Needs nobody in the lab: no step blocks the signal arm.

    python photon_hunt.py                      # full hunt with defaults
    python photon_hunt.py --target-readouts 500
    python photon_hunt.py --help

Steps (each logged to report.md, each figure saved as PNG in the same folder):
  0. open the camera, wait for -20 C Locked, check that settings can be written
  1. trigger-rate test (camera timing only)
  2. wide-gate reference: is the beam spot on the sensor at all? (abort if not)
  3. coarse delay scan with 1 us windows: which window holds the heralded photons?
     (unheralded light is the same in every window, so only the herald stands out)
  4. halving ladder from that window down to ~10 ns (early vs late half each step;
     a small delay scan resolves a step whose photons straddle the boundary)
  5. fine delay scan with the final gate (3 ns) to find the peak
  6. long target image at the peak + off-delay background, with all figures
The camera is left in a safe state at the end (lowest gain, intensifier off).

Output: Desktop/HOI_freshstaret/data/processed/troubleshoot_figures/hunt_<timestamp>/
Runs (.npz) go to the usual troubleshoot data folder, so compare_runs.ipynb sees them.
"""

import argparse
import datetime
import json
import sys
import time
import traceback
import pathlib

import numpy as np
import matplotlib
if __name__ == "__main__":
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import hoi_troubleshoot as ht


class HuntStop(Exception):
    """Raised to end the hunt early with a clear reason (not a crash)."""


class Report:
    def __init__(self, folder):
        self.folder = folder
        self.path = folder / "report.md"
        self.data = {"steps": []}
        self.lines = []

    def log(self, text=""):
        print(text, flush=True)
        self.lines.append(text)
        self.path.write_text("\n".join(self.lines) + "\n", encoding="utf-8")

    def step(self, name, **info):
        self.data["steps"].append({"step": name, **{k: ht._jsonable(v) for k, v in info.items()}})
        (self.folder / "report.json").write_text(json.dumps(self.data, indent=1), encoding="utf-8")


# --------------------------------------------------------------------------------------
# Scoring helpers: compare time windows using the spot aperture (no dark frame needed)
# --------------------------------------------------------------------------------------

def window_series(run, spot, thr):
    """(hits per readout, ring-subtracted spot ADU per readout) for one run."""
    ap = ht.spot_masks(run, spot)[0]
    return ht.hit_series(run, thr, ap), ht.spot_series(run, spot)


def pick_best(runs, spot, thr, min_hits=20):
    """Index of the brightest run, and its significance against all the others pooled.
    Uses photon hits when there are enough of them, otherwise spot ADU."""
    H, A = zip(*[window_series(r, spot, thr) for r in runs])
    use_hits = sum(h.sum() for h in H) >= min_hits
    S = H if use_hits else A
    means = np.array([s.mean() for s in S])
    i = int(np.argmax(means))
    others = np.concatenate([S[j] for j in range(len(S)) if j != i]) if len(S) > 1 else np.array([0.0, 0.0])
    z = ht._diff_z(S[i], others)
    return i, float(z), ("hits" if use_hits else "ADU"), means


def merge_runs(a, b):
    """Pool the per-readout data of two runs with identical settings (for a second pass)."""
    out = dict(a)
    out["roi_stack"] = np.concatenate([a["roi_stack"], b["roi_stack"]])
    out["roi_totals"] = np.concatenate([a["roi_totals"], b["roi_totals"]])
    out["sum"] = a["sum"] + b["sum"]
    out["meta"] = json.loads(json.dumps(a["meta"]))
    out["meta"]["n_readouts_done"] = a["meta"]["n_readouts_done"] + b["meta"]["n_readouts_done"]
    return out


def check_saturation(run, rep, limit=1e-3):
    sat = float((run["first"] >= ht.SATURATION_ADU).mean())
    if sat > limit:
        raise HuntStop(f"{100 * sat:.2f}% of pixels saturated in '{run['meta']['settings']['label']}'. "
                       "Stopping to protect the intensifier; lower the gain or accumulations.")


# --------------------------------------------------------------------------------------
# The hunt
# --------------------------------------------------------------------------------------

def hunt(cam=None, args=None):
    a = args or parse_args([])
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    folder = ht.FIG_DIR / f"hunt_{stamp}"
    folder.mkdir(parents=True, exist_ok=True)
    fig_dir_before = ht.FIG_DIR
    ht.FIG_DIR = folder
    rep = Report(folder)
    rep.log(f"# Heralded-photon hunt {stamp}")
    rep.log(f"Settings: {vars(a)}")
    opened_here = cam is None
    base = dict(gain=a.gain, gate_delay_ns=309, gate_width_ns=100_000, accumulations=a.accumulations,
                n_readouts=20, trigger_source="External", trigger_threshold_v=a.trigger_threshold,
                trigger_edge="Rising Edge", trigger_termination="Fifty Ohms", trigger_coupling="DC",
                photon_mode="Disabled", roi_center=tuple(a.roi_center), roi_size=150,
                idler_rate_khz=a.idler_khz, notes="photon_hunt")
    result = {"found": False}
    try:
        # 0. camera ------------------------------------------------------------------
        rep.log("\n## 0. Camera")
        if cam is None:
            cam, limits = ht.open_camera()
        else:
            import camera_control as cc
            limits = cc.get_repetitive_gate_constraint(cam)
        if not ht.cool_down(cam, setpoint_c=-20.0, poll_s=30):
            raise HuntStop("Sensor did not lock at -20 C.")
        ht.check_settable(cam)
        state = ht.read_state(cam)
        rep.log(f"Camera state: {state}")
        min_w = float(limits["width_min_ns"])

        # 1. trigger rate --------------------------------------------------------------
        rep.log("\n## 1. Trigger rate")
        rate = ht.measure_trigger_rate(cam, base, accumulations=a.rate_test_accumulations, width_ns=3,
                                       n_readouts=4, gate_delay_ns=309)
        rep.log(f"Camera-timed trigger rate: {rate / 1e3:.1f} kHz (lower bound if readout-limited). "
                + (f"Idler counter: {a.idler_khz:g} kHz." if a.idler_khz else "Idler counter: not given."))
        rep.step("trigger_rate", rate_hz=rate)
        if rate and rate > 3 * 20e3:
            rep.log("WARNING: the camera triggers far faster than ~15 kHz of idler photons. Most gates may "
                    "have no partner photon; the hunt continues, but check TRIGGER IN on a scope.")

        # 2. wide-gate reference -------------------------------------------------------
        rep.log("\n## 2. Wide-gate reference (is the beam on the sensor?)")
        ref_acc = int(min(a.accumulations, ht.MAX_GATE_OPEN_S / (base["gate_width_ns"] * 1e-9)))
        ref = ht.take_run(cam, **{**base, "label": "hunt_ref_wide", "n_readouts": a.ref_readouts,
                                  "accumulations": ref_acc})
        check_saturation(ref, rep)
        spot = tuple(a.spot) if a.spot else ht.find_spot(ref, None)
        s = ht.spot_series(ref, spot)
        z_ref = s.mean() / (s.std(ddof=1) / np.sqrt(len(s)))
        contrast = ht.spot_contrast(ref, spot)
        rep.log(f"Spot at {spot}: {s.mean():.0f} ring-subtracted ADU per readout, {z_ref:.1f} sigma; "
                f"compactness contrast {contrast:.1f} (beam ~40 on 20 Oct, stray light only ~3-4).")
        ht.plot_run(ref, None, threshold=a.hit_thr, spot=spot, save_as="2_wide_gate_reference")
        rep.step("reference", spot=spot, spot_adu=s.mean(), z=z_ref, contrast=contrast)
        if contrast < a.min_contrast:
            msg = (f"WARNING: no compact beam spot at the wide gate (contrast {contrast:.1f} < {a.min_contrast:g}). "
                   "Either the beam is blocked, or the pump is low so there is little unheralded light to see "
                   "(the heralded signal per readout does not depend on pump power).")
            if a.spot is None:
                raise HuntStop(msg + " Re-run with the known beam position, e.g. --spot 501 454, "
                               "or with more --ref-readouts.")
            rep.log(msg + f" Continuing at the given spot {spot}.")

        # 3. coarse delay scan ---------------------------------------------------------
        rep.log("\n## 3. Coarse delay scan")
        W = float(a.coarse_width)
        delays = np.arange(max(a.coarse_start, limits["delay_min_ns"]), a.coarse_stop, W)
        coarse = ht.delay_scan(cam, base, delays, W, n_readouts=a.coarse_readouts, label="hunt_coarse")
        for r in coarse:
            check_saturation(r, rep)
        noise = np.median([ht.dark_level(r)[1] for r in coarse])
        thr = max(a.hit_thr, 6 * noise)
        rep.log(f"Narrow-gate pixel noise {noise:.1f} ADU -> photon threshold {thr:g} ADU above pedestal.")
        i, z, metric, means = pick_best(coarse, spot, thr)
        if z < a.min_z and a.coarse_passes > 1:
            rep.log(f"Pass 1 best window {delays[i]:g} ns at only {z:.1f} sigma ({metric}); "
                    f"taking a second pass with {a.coarse_readouts * 3} readouts per window.")
            more = ht.delay_scan(cam, base, delays, W, n_readouts=a.coarse_readouts * 3, label="hunt_coarse2")
            coarse = [merge_runs(x, y) for x, y in zip(coarse, more)]
            i, z, metric, means = pick_best(coarse, spot, thr)
        ht.plot_delay_scan(coarse, threshold=thr, spot=spot, save_as="3_coarse_delay_scan")
        rep.log(f"Best window: {delays[i]:g} ns + {ht.fmt_ns(W)}, {z:.1f} sigma above the other windows "
                f"({metric}). Per-window means: {np.round(means, 2).tolist()}")
        rep.step("coarse", delays=delays, best=delays[i], z=z, metric=metric, means=means, threshold=thr)
        if z < a.min_z:
            raise HuntStop(f"No window stands out (best {z:.1f} sigma < {a.min_z:g}). No heralded photons "
                           f"found between {delays[0]:g} and {delays[-1] + W:g} ns. See the recommendations.")

        # 4. halving ladder ------------------------------------------------------------
        rep.log("\n## 4. Halving ladder")
        d, w = float(delays[i]), W
        steps = 0
        while w / 2 >= max(a.ladder_stop_width, min_w) and steps < 20:
            steps += 1
            ea, la = ht.take_halves(cam, base, d, w, n_readouts=a.ladder_readouts, label="hunt_ladder")
            res = ht.compare_halves(ea, la, z_min=3.0, threshold=thr, spot=spot,
                                    save_as=f"4_ladder_{steps:02d}_{d:g}_{w:g}")
            rep.log(f"Window {d:g} ns + {ht.fmt_ns(w)}: {res['pick']} (z {res['z']:+.1f}, hits z "
                    f"{res['z_hits'] if res['z_hits'] is None else round(res['z_hits'], 1)})")
            rep.step("ladder", delay=d, width=w, **res)
            if res["pick"] != "unclear":
                d, w = res["next_delay_ns"], res["next_width_ns"]
                continue
            # photons may straddle the boundary: scan half-width gates across the window
            sub = np.arange(d - w / 4, d + w + 1e-9, w / 4)
            scan = ht.delay_scan(cam, base, sub, w / 2, n_readouts=a.ladder_readouts, label=f"hunt_straddle{steps:02d}")
            j, zj, metric_j, _ = pick_best(scan, spot, thr)
            ht.plot_delay_scan(scan, threshold=thr, spot=spot, save_as=f"4_ladder_{steps:02d}_straddle_scan")
            rep.log(f"  straddle scan: best {sub[j]:g} ns + {ht.fmt_ns(w / 2)} at {zj:.1f} sigma ({metric_j})")
            if zj < 3:
                raise HuntStop(f"Lost the signal at {ht.fmt_ns(w)} windows (straddle scan {zj:.1f} sigma). "
                               "Re-run with more --ladder-readouts.")
            d, w = float(sub[j]), w / 2

        # 5. fine delay scan -----------------------------------------------------------
        rep.log("\n## 5. Fine delay scan")
        fw = max(a.final_width, min_w)
        centre = d + w / 2
        fine_delays = np.arange(centre - w / 2 - 2 * fw, centre + w / 2 + 2 * fw + 1e-9, fw / 2)
        fine = ht.delay_scan(cam, base, fine_delays, fw, n_readouts=a.fine_readouts, label="hunt_fine")
        j, zf, metric_f, means_f = pick_best(fine, spot, thr)
        ht.plot_delay_scan(fine, threshold=thr, spot=spot, save_as="5_fine_delay_scan")
        best_delay = float(fine_delays[j])
        rep.log(f"Peak at {best_delay:g} ns with a {ht.fmt_ns(fw)} gate, {zf:.1f} sigma ({metric_f}).")
        rep.step("fine", delays=fine_delays, best=best_delay, z=zf, metric=metric_f, means=means_f)
        result.update(found=True, delay_ns=best_delay, width_ns=fw, spot=spot)

        # 6. target image --------------------------------------------------------------
        if a.target_readouts > 0:
            rep.log("\n## 6. Target image")
            target = {**base, "gate_delay_ns": best_delay, "gate_width_ns": fw,
                      "accumulations": a.target_accumulations, "n_readouts": a.target_readouts}
            runT = ht.take_run(cam, **{**target, "label": "hunt_target"})
            check_saturation(runT, rep)
            bg = ht.take_offdelay_background(cam, runT, shift_ns=a.offdelay_shift,
                                             n_readouts=max(a.target_readouts // 2, 10), label="hunt_target_offdelay")
            ht.plot_run(runT, bg, threshold=thr, spot=spot, save_as="6_target_image")
            ht.plot_raw_vs_accumulated(runT, bg, save_as="6_target_raw_vs_sum")
            ht.plot_wiener(runT, bg, save_as="6_target_wiener")
            ht.threshold_sweep(runT, bg, spot=spot, save_as="6_target_threshold_sweep")
            m = ht.run_metrics(runT, bg, spot=spot, threshold=thr)
            rep.log(f"Target: {m.get('spot_net_adu_per_readout', float('nan')):.0f} net spot ADU per readout "
                    f"({m.get('spot_z', float('nan')):.1f} sigma above the off-delay background); hits "
                    f"{m.get('spot_hits_per_readout', float('nan')):.2f} vs "
                    f"{m.get('spot_dark_hits_per_readout', float('nan')):.2f} per readout "
                    f"({m.get('spot_hits_z', float('nan')):.1f} sigma).")
            rep.step("target", **{k: v for k, v in m.items() if isinstance(v, (int, float, str, type(None)))})
    except HuntStop as e:
        rep.log(f"\n## Stopped\n{e}")
        result["stopped"] = str(e)
    except KeyboardInterrupt:
        rep.log("\n## Stopped by user")
        result["stopped"] = "user"
    except Exception as e:
        rep.log(f"\n## Error\n```\n{traceback.format_exc()}\n```")
        result["stopped"] = f"error: {e}"
    finally:
        plt.close("all")
        rep.log("\n## Result")
        rep.log(json.dumps(ht._jsonable(result)))
        if not result.get("found"):
            rep.log("\nRecommendations if nothing was found:\n"
                    "- Compare the camera-timed trigger rate with the idler counter; scope TRIGGER IN.\n"
                    "- Measure the idler-signal coincidence delay on the time tagger with the far end of "
                    "the delay fibre on an APD; then set the gate from that delay minus ~35 ns.\n"
                    "- Re-run with more readouts (--coarse-readouts) or a wider range (--coarse-stop).")
        rep.step("result", **result)
        if cam is not None and opened_here and not a.no_shutdown:
            try:
                ht.shutdown(cam)
            except Exception as e:
                rep.log(f"(shutdown: {e})")
        ht.FIG_DIR = fig_dir_before
        rep.log(f"\nReport and figures: {folder}")
    return result, folder


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Unattended heralded-photon hunt on the PI-MAX4.")
    p.add_argument("--gain", type=int, default=9000)
    p.add_argument("--accumulations", type=int, default=5000, help="per readout during the search")
    p.add_argument("--trigger-threshold", type=float, default=0.5)
    p.add_argument("--idler-khz", type=float, default=None, help="idler counter reading, for the report")
    p.add_argument("--roi-center", type=int, nargs=2, default=(500, 455))
    p.add_argument("--spot", type=int, nargs=2, default=None, help="beam spot (x y); default: find it")
    p.add_argument("--hit-thr", type=float, default=200, help="photon threshold, ADU above pedestal")
    p.add_argument("--ref-readouts", type=int, default=20, help="readouts for the wide-gate beam check")
    p.add_argument("--min-contrast", type=float, default=10,
                   help="beam compactness needed at the wide gate (beam ~40, stray light ~3-4)")
    p.add_argument("--rate-test-accumulations", type=int, default=100_000)
    p.add_argument("--coarse-start", type=float, default=27)
    p.add_argument("--coarse-stop", type=float, default=10_000)
    p.add_argument("--coarse-width", type=float, default=1000)
    p.add_argument("--coarse-readouts", type=int, default=40)
    p.add_argument("--coarse-passes", type=int, default=2)
    p.add_argument("--min-z", type=float, default=5.0, help="significance needed to accept a window")
    p.add_argument("--ladder-readouts", type=int, default=30)
    p.add_argument("--ladder-stop-width", type=float, default=10)
    p.add_argument("--final-width", type=float, default=3)
    p.add_argument("--fine-readouts", type=int, default=40)
    p.add_argument("--target-accumulations", type=int, default=50_000)
    p.add_argument("--target-readouts", type=int, default=200, help="0 skips the target image")
    p.add_argument("--offdelay-shift", type=float, default=5000)
    p.add_argument("--no-shutdown", action="store_true", help="leave the camera configured at the end")
    return p.parse_args(argv)


if __name__ == "__main__":
    result, folder = hunt(args=parse_args())
    sys.exit(0 if result.get("found") else 1)
