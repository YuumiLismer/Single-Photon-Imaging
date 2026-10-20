"""
hoi_troubleshoot.py -- acquisition, bookkeeping and plotting for troubleshoot.ipynb
and compare_runs.ipynb (heralded single-photon imaging on the PI-MAX4 emICCD).

Design
------
* Camera control reuses the hardware-tested helpers in
  Desktop/HOI_freshstaret/src/camera_control.py (raw-DLL Repetitive Gate setter,
  PICam race-condition retries, stop_acquisition patch). Nothing camera-level is
  re-implemented here.
* Every acquisition is a "run": a dict with the summed full frame, the first raw
  readout, a stack of per-readout ROI crops, per-readout totals, and a metadata dict
  holding BOTH the requested settings and what the camera reports back.
  Runs are saved as .npz (one file per run) so every plot can be rebuilt later
  without the camera (compare_runs.ipynb).
* The ROI is a software crop of a full-sensor readout, so each run has both the
  full frame and the zoomed ROI.

Units: camera values are ADU ("counts"). "per readout" means divided by the number
of camera readouts in the run; each readout itself contains `accumulations` gated
exposures (one gate per idler trigger) summed on the CCD.
"""

import sys
import json
import time
import datetime
import re
import pathlib

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
from scipy.ndimage import gaussian_filter

# --------------------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------------------

HOI_ROOT = pathlib.Path(r"C:\Users\ReschLab\Desktop\HOI_freshstaret")
# Every experiment gets its own folder, sorted by date and start time:
#   data/troubleshoot/2026-10-20/1614_agent_exp5_on_off_test/runs/      (.npz runs)
#                                                           /figures/   (.png, reports)
# Call new_experiment("name") before a block of runs. If a run or figure is saved before
# that, an experiment called "session" is started automatically.
DATA_ROOT = HOI_ROOT / "data" / "troubleshoot"
EXPERIMENT_DIR = None           # folder of the current experiment (set by new_experiment)
DATA_DIR = None                 # EXPERIMENT_DIR / "runs"
FIG_DIR = None                  # EXPERIMENT_DIR / "figures"

if str(HOI_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(HOI_ROOT / "src"))

SATURATION_ADU = 65000          # 16-bit ADC; flag pixels at or above this
MAX_GAIN = 9000                 # emICCD gain safety limit (never exceed)
MAX_GATE_OPEN_S = 1.0           # above gain 1000: gate width x accumulations per readout
SESSION_RUNS = []               # every run taken (or loaded) in this kernel, in order



# --------------------------------------------------------------------------------------
# Experiment folders
# --------------------------------------------------------------------------------------

def _slug(name):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name)).strip("_") or "experiment"


def new_experiment(name, when=None, quiet=False):
    """Start a new experiment folder DATA_ROOT/<YYYY-MM-DD>/<HHMM>_<name>/ with runs/ and
    figures/ inside, and send every following run and figure there. Returns the folder."""
    when = when or datetime.datetime.now()
    day = DATA_ROOT / f"{when:%Y-%m-%d}"
    folder = day / f"{when:%H%M}_{_slug(name)}"
    k = 2
    while folder.exists():
        folder = day / f"{when:%H%M}_{_slug(name)}_{k}"
        k += 1
    (folder / "runs").mkdir(parents=True)
    (folder / "figures").mkdir()
    set_experiment(folder, quiet=quiet)
    return folder


def set_experiment(folder, quiet=False):
    """Point runs and figures at an existing experiment folder (None unsets it)."""
    global EXPERIMENT_DIR, DATA_DIR, FIG_DIR
    if folder is None:
        EXPERIMENT_DIR = DATA_DIR = FIG_DIR = None
        return None
    EXPERIMENT_DIR = pathlib.Path(folder)
    DATA_DIR = EXPERIMENT_DIR / "runs"
    FIG_DIR = EXPERIMENT_DIR / "figures"
    if not quiet:
        print(f"experiment folder: {EXPERIMENT_DIR}")
    return EXPERIMENT_DIR


class experiment:
    """`with ht.experiment("name") as folder:` starts a new experiment folder and puts the
    previous one back afterwards (used by photon_hunt and the agent scripts)."""
    def __init__(self, name):
        self.name = name

    def __enter__(self):
        self.before = EXPERIMENT_DIR
        return new_experiment(self.name)

    def __exit__(self, *exc):
        set_experiment(self.before, quiet=True)
        return False


def _ensure_experiment():
    if EXPERIMENT_DIR is None:
        new_experiment("session")


def list_experiments(contains=""):
    """All experiment folders, oldest first."""
    return sorted(p for p in DATA_ROOT.glob("*/*") if p.is_dir() and contains in p.name)

# --------------------------------------------------------------------------------------
# Settings
# --------------------------------------------------------------------------------------

DEFAULT_SETTINGS = dict(
    label="run",
    kind="signal",              # "signal", "dark" (signal arm blocked), "dark_offdelay", "internal"
    gain=9000,                  # emICCD gain
    gate_delay_ns=309.0,        # Repetitive Gate delay after the trigger
    gate_width_ns=100000.0,     # Repetitive Gate width (100 us = the wide gate from LightField)
    accumulations=5000,         # on-CCD accumulations = gated exposures summed per readout
    n_readouts=20,              # readouts summed in software (LightField "exposures per frame")
    trigger_source="External",  # "External" (idler herald) or "Internal"
    trigger_threshold_v=0.5,
    trigger_edge="Rising Edge",
    trigger_termination="Fifty Ohms",
    trigger_coupling="DC",
    internal_trigger_hz=15000,  # only used when trigger_source == "Internal"
    photon_mode="Disabled",     # "Disabled" (raw ADU), "Thresholding" (0/1), "Clipping"
    photon_threshold=None,      # ADU, only used when photon_mode != "Disabled"
    roi_center=(512, 512),      # (x, y) pixel of the ROI centre on the full sensor
    roi_size=150,
    keep_roi_frames=400,        # per-readout ROI crops kept for histograms / threshold sweeps
    max_seconds=None,           # optional wall-clock limit for the run
    frame_timeout_s=60,         # stop if no readout arrives for this long (no triggers?)
    idler_rate_khz=None,        # type the counter reading; the camera cannot measure it
    notes="",
)

READBACK_ATTRS = [
    "Sensor Temperature Reading", "Sensor Temperature Status", "Sensor Temperature Set Point",
    "emICCD Gain", "emICCD Gain Control Mode", "Intensifier Gain", "EM Gain",
    "Accumulations", "Enable Intensifier", "Intensifier Status", "Gating Mode",
    "Trigger Source", "Trigger Response", "Trigger Determination", "Trigger Threshold",
    "Trigger Termination", "Trigger Coupling", "Trigger Frequency",
    "Photon Detection Mode", "Photon Detection Threshold",
    "ADC Speed", "ADC Quality", "ADC Analog Gain", "Readout Control Mode",
    "Readout Time Calculation", "Frame Rate Calculation",
]


def _jsonable(v):
    if isinstance(v, (bool, int, float, str)) or v is None:
        return v
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return float(v)
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, dict):
        return {str(k): _jsonable(x) for k, x in v.items()}
    return str(v)


def fmt_ns(x):
    """3 -> '3 ns', 2500 -> '2.5 us', 100000 -> '100 us'."""
    if x is None:
        return "?"
    x = float(x)
    if x >= 1000:
        return f"{x / 1000:g} us"
    return f"{x:g} ns"


# --------------------------------------------------------------------------------------
# Camera: open, cool, configure, read back, shut down
# --------------------------------------------------------------------------------------

def open_camera():
    """Open the PI-MAX4 (LightField must be closed) and print what it is.
    Returns (cam, gate_limits)."""
    import camera_control as cc
    cam = cc.connect_camera()
    print("Opened:", cam.get_device_info())
    print("Detector size (w, h):", cam.get_detector_size())
    limits = cc.get_repetitive_gate_constraint(cam)
    print(f"Gate limits: delay {limits['delay_min_ns']:g} .. {limits['delay_max_ns']:g} ns, "
          f"width {limits['width_min_ns']:g} .. {limits['width_max_ns']:g} ns")
    return cam, limits


def cool_down(cam, setpoint_c=-20.0, tol_c=0.5, timeout_min=45, poll_s=15):
    """Set the sensor set point and wait until the camera reports 'Locked' within tol_c.
    Plots the cooling curve. Returns True when locked."""
    import camera_control as cc
    cc.safe_set_attribute_value(cam, "Sensor Temperature Set Point", setpoint_c)
    cc.safe_commit_parameters(cam)
    hist, t0, locked = [], time.time(), False
    while True:
        temp = float(cam.get_attribute_value("Sensor Temperature Reading"))
        status = str(cam.get_attribute_value("Sensor Temperature Status"))
        mins = (time.time() - t0) / 60
        hist.append((mins, temp))
        print(f"  {mins:5.1f} min   {temp:7.2f} C   {status}")
        if status.lower() == "locked" and abs(temp - setpoint_c) <= tol_c:
            locked = True
            print(f"LOCKED at {temp:.2f} C.")
            break
        if mins > timeout_min:
            print(f"TIMEOUT: still {temp:.2f} C / {status}. Check the fan and airflow.")
            break
        time.sleep(poll_s)
    if len(hist) > 1:
        xs, ys = zip(*hist)
        plt.figure(figsize=(6, 2.8))
        plt.plot(xs, ys, "o-")
        plt.axhline(setpoint_c, color="r", ls="--", label=f"set point {setpoint_c:g} C")
        plt.xlabel("time (min)"); plt.ylabel("sensor temperature (C)")
        plt.title("Sensor cool-down"); plt.legend(); plt.grid(alpha=.3)
        plt.tight_layout(); plt.show()
    return locked


def read_state(cam):
    """Dict of what the camera currently reports for every attribute in READBACK_ATTRS,
    plus the real Repetitive Gate (read through the raw DLL)."""
    import camera_control as cc
    attrs = cam.get_all_attributes()
    out = {}
    for name in READBACK_ATTRS:
        if name not in attrs:
            continue
        try:
            out[name] = _jsonable(cam.get_attribute_value(name))
        except Exception as e:
            out[name] = f"<unreadable: {e}>"
    try:
        d, w = cc.get_repetitive_gate(cam)
        out["Repetitive Gate delay (ns)"] = float(d)
        out["Repetitive Gate width (ns)"] = float(w)
    except Exception as e:
        out["Repetitive Gate"] = f"<unreadable: {e}>"
    return out


def show_state(cam):
    """Print the camera state as a two-column table (compare against LightField)."""
    for k, v in read_state(cam).items():
        print(f"  {k:32s} {v}")


ENUM_SETTINGS = ["Trigger Source", "Trigger Response", "Gating Mode", "Trigger Determination",
                 "Trigger Termination", "Trigger Coupling", "Photon Detection Mode"]


def _set_enum(cam, attrs, name, label, retries=100, wait_s=0.1, nudge_every=10):
    """Set an enum attribute by label. Retries only on PICam 'acquisition in progress'
    (error 20, seen on this camera on 20 Oct for every trigger/gating setting), re-sending
    stop_acquisition every `nudge_every` tries like camera_control.safe_set_attribute_value;
    falls back to the integer index. Returns (ok, error_text)."""
    import camera_control as cc
    busy = None
    try:
        busy = cc.picam_defs.PicamError.PicamError_AcquisitionInProgress
    except Exception:
        pass
    last = None
    for attempt in range(retries):
        try:
            cam.set_attribute_value(name, label)
            return True, None
        except Exception as e:
            last = e
            if busy is not None and getattr(e, "code", None) == busy:
                if attempt and attempt % nudge_every == 0:
                    try:
                        cam.stop_acquisition()
                    except Exception:
                        pass
                time.sleep(wait_s)
                continue
            break
    idx = getattr(attrs.get(name), "labels", {}) or {}
    if isinstance(idx, dict) and label in idx:
        try:
            cam.set_attribute_value(name, idx[label])
            return True, None
        except Exception as e:
            last = e
    return False, f"{type(last).__name__}: {last}"


def check_settable(cam):
    """Re-send the CURRENT value of each enum setting and print what happens.
    Changes nothing on the camera; shows the real error if a setting cannot be written."""
    attrs = cam.get_all_attributes()
    for name in ENUM_SETTINGS:
        if name not in attrs:
            print(f"  SKIP {name}: not on this camera")
            continue
        cur = cam.get_attribute_value(name)
        ok, err = _set_enum(cam, attrs, name, cur, retries=1)
        print(f"  {'OK  ' if ok else 'FAIL'} {name} = {cur!r}" + ("" if ok else f"   -> {err}"))


def configure(cam, s, verbose=True, strict=True):
    """Push settings dict `s` to the camera, then read every setting back.
    With strict=True, raises RuntimeError (before any data is taken) if the camera did not
    accept a trigger, gating, photon-detection, gain, accumulation or gate setting."""
    import camera_control as cc
    attrs = cam.get_all_attributes()      # cached once (repeated calls once hit RecursionError)
    lines = []
    wanted = {}

    def enum(name, label):
        if name not in attrs:
            lines.append(f"  SKIP {name}: not on this camera")
            return
        ok, err = _set_enum(cam, attrs, name, label)
        wanted[name] = label
        lines.append(f"  {'OK  ' if ok else 'FAIL'} {name} -> {label}" + ("" if ok else f"   ({err})"))

    def num(name, value, integer=False):
        if name not in attrs:
            lines.append(f"  SKIP {name}: not on this camera")
            return
        a, v = attrs[name], value
        lo, hi = getattr(a, "min", None), getattr(a, "max", None)
        if isinstance(lo, (int, float)) and isinstance(hi, (int, float)) and hi > lo:
            v = min(max(v, lo), hi)
        v = int(round(v)) if integer else float(v)
        cc.safe_set_attribute_value(cam, name, v)
        if integer:
            wanted[name] = v
        note = "" if v == value else f"   (asked {value}, clipped to camera range)"
        lines.append(f"  OK   {name} -> {v}{note}")

    enum("Trigger Source", s["trigger_source"])
    enum("Trigger Response", "Readout Per Trigger")
    enum("Gating Mode", "Repetitive")
    if s["trigger_source"] == "External":
        enum("Trigger Determination", s["trigger_edge"])
        enum("Trigger Termination", s["trigger_termination"])
        enum("Trigger Coupling", s["trigger_coupling"])
        num("Trigger Threshold", s["trigger_threshold_v"])
    else:
        num("Trigger Frequency", s["internal_trigger_hz"])
    if "Enable Intensifier" in attrs:
        cc.safe_set_attribute_value(cam, "Enable Intensifier", True)
        lines.append("  OK   Enable Intensifier -> True")
    num("emICCD Gain", s["gain"], integer=True)
    num("Accumulations", s["accumulations"], integer=True)
    enum("Photon Detection Mode", s["photon_mode"])
    if s["photon_mode"] != "Disabled" and s["photon_threshold"] is not None:
        num("Photon Detection Threshold", s["photon_threshold"])
    cc.safe_commit_parameters(cam)

    limits = cc.get_repetitive_gate_constraint(cam)
    d, w = cc.set_repetitive_gate(cam, s["gate_delay_ns"], s["gate_width_ns"], clip_to=limits)
    cc.safe_commit_parameters(cam)
    rd, rw = cc.get_repetitive_gate(cam)
    tol = max(0.5, float(limits.get("width_increment_ns") or 0))
    if abs(rd - d) > tol or abs(rw - w) > tol:
        raise RuntimeError(f"Gate read-back ({rd}, {rw}) ns does not match the request ({d}, {w}) ns.")
    lines.append(f"  OK   Repetitive Gate -> delay {rd:g} ns, width {fmt_ns(rw)} (read back)")

    # Read every setting back: the read-back, not the set call, is what the camera will use.
    wrong = []
    for name, want in wanted.items():
        try:
            got = cam.get_attribute_value(name)
        except Exception as e:
            got = f"<unreadable: {e}>"
        if str(got) != str(want):
            wrong.append(f"{name}: wanted {want!r}, camera has {got!r}")
    if wrong:
        lines.append("  NOT APPLIED: " + "; ".join(wrong))
    if verbose:
        print("\n".join(lines))
    else:   # always surface anything that did not go as asked
        for ln in lines:
            if "FAIL" in ln or "SKIP" in ln or "clipped" in ln or "NOT APPLIED" in ln:
                print(ln)
    if wrong and strict:
        raise RuntimeError(
            "The camera did not take these settings, so no data was taken (the run would be "
            "mislabelled):\n  " + "\n  ".join(wrong) +
            "\nRun ht.check_settable(cam) to see the camera's error for each setting.")
    if any("FAIL" in ln for ln in lines) and not wrong:
        print("  (the FAILed settings already had the requested value, so this run is still valid)")


def shutdown(cam):
    """Safe state: lowest gain, intensifier off, internal trigger, then close."""
    import camera_control as cc
    attrs = cam.get_all_attributes()
    for name, value in (("emICCD Gain", getattr(attrs.get("emICCD Gain"), "min", 1)),
                        ("Enable Intensifier", False)):
        try:
            cc.safe_set_attribute_value(cam, name, value)
            print(f"  {name} -> {value}")
        except Exception as e:
            print(f"  {name}: {e}")
    try:
        cc.set_enum_attr(cam, "Trigger Source", "Internal")
        cc.safe_commit_parameters(cam)
        print("  Trigger Source -> Internal")
    except Exception as e:
        print(f"  Trigger Source: {e}")
    cam.close()
    print("Camera closed.")


# --------------------------------------------------------------------------------------
# Acquisition
# --------------------------------------------------------------------------------------

def roi_box(center, size, det_w, det_h):
    """(x0, y0, size) of a size x size box centred on `center`, clipped to the sensor."""
    size = int(min(size, det_w, det_h))
    x0 = int(np.clip(int(center[0]) - size // 2, 0, det_w - size))
    y0 = int(np.clip(int(center[1]) - size // 2, 0, det_h - size))
    return x0, y0, size


def _acquire(cam, n_readouts, roi, keep_roi_frames, max_seconds, frame_timeout_s):
    """Read `n_readouts` full-sensor readouts. Returns a dict of arrays plus stop reason.
    Ctrl-C / the stop button ends the run early and keeps what was taken."""
    import camera_control as cc
    x0, y0, n = roi
    total = first = None
    roi_frames, roi_tot, full_tot, stamps = [], [], [], []
    stopped_by = "done"
    t0 = last = time.time()
    cam.start_acquisition()
    try:
        while len(roi_tot) < n_readouts:
            if max_seconds and time.time() - t0 > max_seconds:
                stopped_by = "max_seconds"
                break
            try:
                got = cc.poll_frame(cam, 0.5)
                imgs = cam.read_multiple_images() if got else []
            except RuntimeError as e:          # occasional pylablib buffer hiccup
                print(f"\n  buffer hiccup ({e}); restarting acquisition")
                try:
                    cam.stop_acquisition()
                except Exception:
                    pass
                cam.start_acquisition()
                continue
            if not imgs:
                if time.time() - last > frame_timeout_s:
                    print(f"\n  No readout for {frame_timeout_s:.0f} s. Are triggers reaching "
                          "the camera? Check the APD cable, threshold and edge.")
                    stopped_by = "timeout"
                    break
                continue
            for im in imgs:
                im = np.asarray(im, dtype=np.float64)
                if first is None:
                    first = im.copy()
                total = im.copy() if total is None else total + im
                crop = im[y0:y0 + n, x0:x0 + n]
                if len(roi_frames) < keep_roi_frames:
                    roi_frames.append(crop.copy())
                roi_tot.append(float(crop.sum()))
                full_tot.append(float(im.sum()))
                stamps.append(time.time() - t0)
                last = time.time()
                if len(roi_tot) >= n_readouts:
                    break
            print(f"\r  {len(roi_tot)}/{n_readouts} readouts, {time.time() - t0:6.1f} s", end="")
    except KeyboardInterrupt:
        stopped_by = "user"
        print("\n  stopped by user; keeping what was taken")
    finally:
        cam.stop_acquisition()
        time.sleep(0.2)                         # driver settle (see camera_control.acquire_sum)
    print()
    stack = np.array(roi_frames) if roi_frames else np.zeros((0, n, n))
    if stack.size and stack.min() >= 0 and stack.max() <= 65535 and np.all(stack == np.round(stack)):
        stack = stack.astype(np.uint16)
    else:
        stack = stack.astype(np.float32)
    return dict(sum=total, first=first, roi_stack=stack, roi_totals=np.array(roi_tot),
                full_totals=np.array(full_tot), t=np.array(stamps),
                wall_s=time.time() - t0, stopped_by=stopped_by)


def _actual(s, rb):
    """Normalised 'what the camera actually used' dict, from the read-back."""
    def g(k, default):
        v = rb.get(k, default)
        return default if isinstance(v, str) and v.startswith("<unreadable") else v
    return dict(
        gain=g("emICCD Gain", s["gain"]),
        gate_delay_ns=g("Repetitive Gate delay (ns)", s["gate_delay_ns"]),
        gate_width_ns=g("Repetitive Gate width (ns)", s["gate_width_ns"]),
        accumulations=g("Accumulations", s["accumulations"]),
        trigger_source=str(g("Trigger Source", s["trigger_source"])),
        trigger_threshold_v=g("Trigger Threshold", s["trigger_threshold_v"]),
        internal_trigger_hz=g("Trigger Frequency", s["internal_trigger_hz"]),
        photon_mode=str(g("Photon Detection Mode", s["photon_mode"])),
        photon_threshold=g("Photon Detection Threshold", s["photon_threshold"]),
        temp_c=g("Sensor Temperature Reading", None),
        temp_status=str(g("Sensor Temperature Status", "?")),
    )


def take_run(cam, save=True, **settings):
    """Configure the camera, acquire, save, and return a run dict.

    Any key of DEFAULT_SETTINGS can be passed, e.g.
        take_run(cam, label="wide_gate", gain=9000, gate_width_ns=1e5, n_readouts=20)
    """
    unknown = set(settings) - set(DEFAULT_SETTINGS)
    if unknown:
        raise ValueError(f"Unknown setting(s): {sorted(unknown)}. Valid: {sorted(DEFAULT_SETTINGS)}")
    s = {**DEFAULT_SETTINGS, **settings}
    open_s = float(s["gate_width_ns"]) * 1e-9 * int(s["accumulations"])
    if float(s["gain"]) > MAX_GAIN:
        raise ValueError(f"gain {s['gain']} exceeds the safety limit {MAX_GAIN}")
    if float(s["gain"]) > 1000 and open_s > MAX_GATE_OPEN_S + 1e-9:
        raise ValueError(f"{fmt_ns(s['gate_width_ns'])} x {int(s['accumulations']):,} acc = {open_s:.2f} s of "
                         f"gate-open time per readout at gain {s['gain']}; the intensifier limit is "
                         f"{MAX_GATE_OPEN_S:g} s above gain 1000. Lower the accumulations or the width.")
    print(f"=== {s['label']} [{s['kind']}] ===")
    configure(cam, s, verbose=False)
    rb = read_state(cam)
    act = _actual(s, rb)
    if act["temp_status"].lower() != "locked":
        print(f"  WARNING: sensor is {act['temp_c']} C / {act['temp_status']} (not locked)")
    det_w, det_h = cam.get_detector_size()
    roi = roi_box(s["roi_center"], s["roi_size"], det_w, det_h)
    print(f"  gain {act['gain']} | gate {act['gate_delay_ns']:g} ns + {fmt_ns(act['gate_width_ns'])} | "
          f"{act['accumulations']:,} acc x {s['n_readouts']} readouts | {act['trigger_source']} | "
          f"PD {act['photon_mode']} | {act['temp_c']} C")
    started = datetime.datetime.now()
    data = _acquire(cam, s["n_readouts"], roi, s["keep_roi_frames"], s["max_seconds"],
                    s["frame_timeout_s"])
    if data["sum"] is None:
        print("  No readouts were taken; nothing saved.")
        return None
    meta = dict(settings=_jsonable(s), actual=_jsonable(act), readback=rb,
                roi=list(roi), detector=[int(det_w), int(det_h)],
                n_readouts_done=int(len(data["roi_totals"])), wall_s=float(data["wall_s"]),
                stopped_by=data["stopped_by"], started=started.isoformat(timespec="seconds"))
    run = dict(meta=meta, **{k: data[k] for k in ("sum", "first", "roi_stack", "roi_totals",
                                                    "full_totals", "t")})
    if save:
        path = save_run(run)
        print(f"  saved {path.name}")
    SESSION_RUNS.append(run)
    m = run_metrics(run)
    print(f"  done: {meta['n_readouts_done']} readouts in {meta['wall_s']:.0f} s | "
          f"ROI mean {m['roi_mean_per_readout']:.1f} ADU/px/readout | saturated px "
          f"{100 * m['saturated_frac']:.3f}% of first readout")
    print("  " + trigger_rate_note(run))
    try:
        from run_logger import log_run
        log_run(stage="troubleshoot_run", params=meta["actual"] | {"label": s["label"], "kind": s["kind"],
                "n_readouts": meta["n_readouts_done"], "roi": meta["roi"]},
                result={k: m[k] for k in ("roi_mean_per_readout", "roi_total", "saturated_frac")},
                notes=s["notes"])
    except Exception as e:
        print(f"  (run_log not written: {e})")
    return run


# --------------------------------------------------------------------------------------
# Save / load
# --------------------------------------------------------------------------------------

def run_name(run):
    m = run["meta"]
    stamp = m["started"].replace(":", "").replace("-", "")
    return f"{stamp}_{m['settings']['label']}_{m['settings']['kind']}"


def save_run(run, folder=None):
    if folder is None:
        _ensure_experiment()
        folder = DATA_DIR
    folder = pathlib.Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f"{run_name(run)}.npz"
    run["meta"]["file"] = str(path)
    np.savez_compressed(path, meta=json.dumps(run["meta"]),
                        **{k: run[k] for k in ("sum", "first", "roi_stack", "roi_totals",
                                               "full_totals", "t")})
    return path


def load_run(path):
    with np.load(path, allow_pickle=False) as z:
        run = {k: z[k] for k in z.files if k != "meta"}
        run["meta"] = json.loads(str(z["meta"]))
    run["meta"]["file"] = str(path)
    return run


def load_runs(folder=None, contains=""):
    """All runs under `folder` (searched recursively; default: every experiment in
    DATA_ROOT) whose file name contains `contains`, oldest first. `folder` can be an
    experiment folder, a day folder or a runs/ folder. Unreadable files (e.g. one still
    being written) are skipped with a note."""
    out = []
    files = pathlib.Path(folder or DATA_ROOT).rglob(f"*{contains}*.npz")
    for p in sorted(files, key=lambda q: q.name):
        try:
            out.append(load_run(p))
        except Exception as e:
            print(f"  skipped {p.name}: {type(e).__name__}: {e}")
    return out


# --------------------------------------------------------------------------------------
# Timing: what the readout interval says about the trigger rate
# --------------------------------------------------------------------------------------

def readout_interval_s(run):
    """Median time between readouts (s), from the per-readout time stamps."""
    t = np.asarray(run["t"], dtype=float)
    return float(np.median(np.diff(t))) if len(t) > 2 else np.nan


def trigger_rate_note(run):
    """One line on what the readout timing implies about the trigger rate.
    Each readout needs `accumulations` triggers, so gates per second is a LOWER bound on the
    trigger rate. If the interval equals the bare readout time, the camera was waiting on
    readout, not on triggers, and the trigger rate is at least the number shown."""
    a = run["meta"]["actual"]
    dt = readout_interval_s(run)
    if not np.isfinite(dt) or dt <= 0:
        return "trigger rate: not enough readouts to estimate"
    acc = int(a["accumulations"])
    rate = acc / dt
    ro = run["meta"]["readback"].get("Readout Time Calculation")
    ro_s = float(ro) / 1000 if isinstance(ro, (int, float)) and ro > 0 else 0.23
    gate_s = acc * float(a["gate_width_ns"]) * 1e-9
    txt = f"{dt:.3f} s per readout -> at least {rate / 1e3:.1f} k gates/s"
    if abs(dt / max(ro_s, gate_s) - 1) < 0.15:
        limit = "readout time" if ro_s >= gate_s else "gate-open time"
        txt += f" (limited by {limit} {max(ro_s, gate_s):.3f} s, so triggers arrive faster than this)"
    return "trigger rate: " + txt


# --------------------------------------------------------------------------------------
# Analysis
# --------------------------------------------------------------------------------------

def n_rd(run):
    return max(int(run["meta"]["n_readouts_done"]), 1)


def crop(run, img):
    x0, y0, n = run["meta"]["roi"]
    return img[y0:y0 + n, x0:x0 + n]


def gates(run):
    return int(run["meta"]["actual"]["accumulations"]) * n_rd(run)


def match_warnings(run, dark):
    """Settings that differ between a signal run and its dark run (they should not)."""
    a, b = run["meta"]["actual"], dark["meta"]["actual"]
    keys = ["gain", "gate_width_ns", "accumulations", "photon_mode", "photon_threshold"]
    if dark["meta"]["settings"]["kind"] == "dark":
        keys += ["gate_delay_ns", "trigger_source"]
    elif dark["meta"]["settings"]["kind"] == "dark_offdelay":
        keys += ["trigger_source"]
    out = [f"{k}: signal {a[k]} vs dark {b[k]}" for k in keys if a.get(k) != b.get(k)]
    if run["meta"]["roi"] != dark["meta"]["roi"]:
        out.append(f"roi: signal {run['meta']['roi']} vs dark {dark['meta']['roi']}")
    return out


def scaled_dark_sum(run, dark):
    """Dark summed image rescaled to the same number of readouts as `run`."""
    return dark["sum"] * (n_rd(run) / n_rd(dark))


def net_full(run, dark=None):
    """Signal minus matched dark (full frame). Without a dark, minus the frame median."""
    if dark is not None:
        return run["sum"] - scaled_dark_sum(run, dark)
    return run["sum"] - np.median(run["sum"])


def net_significance(run, dark):
    """How many standard deviations the ROI total sits above the dark, using the
    readout-to-readout scatter of each run as the noise estimate."""
    s, d = run["roi_totals"], dark["roi_totals"]
    ns, nd = len(s), len(d)
    if ns < 2 or nd < 2:
        return np.nan
    net = s.sum() - d.sum() * ns / nd
    sig_s = s.std(ddof=1) * np.sqrt(ns)
    sig_d = d.std(ddof=1) * np.sqrt(nd) * ns / nd
    return float(net / np.hypot(sig_s, sig_d))


# --- spot aperture and photon hits -----------------------------------------------------
# The beam is a small spot (~15 px) in a 150 px ROI, so whole-ROI numbers dilute it ~100x.
# Spot numbers use a circle around the spot, minus the median of a ring around it in the
# SAME readout. That ring subtraction also cancels pedestal drift between runs.

SPOT_RADIUS_PX = 10
RING_PX = (20, 35)


def roi_center_xy(run):
    x0, y0, n = run["meta"]["roi"]
    return (x0 + n // 2, y0 + n // 2)


def find_spot(run, dark=None, sigma=3.0):
    """(x, y) on the full sensor of the brightest smoothed spot in the net ROI image."""
    x0, y0, _ = run["meta"]["roi"]
    py, px = find_peak(crop(run, net_full(run, dark)), sigma)
    return (int(x0 + px), int(y0 + py))


def spot_masks(run, spot=None):
    """Boolean (aperture, ring) masks in ROI coordinates. spot = (x, y) on the full sensor;
    default is the ROI centre."""
    x0, y0, n = run["meta"]["roi"]
    sx, sy = spot if spot is not None else roi_center_xy(run)
    yy, xx = np.indices((n, n))
    r = np.hypot(yy - (sy - y0), xx - (sx - x0))
    ap, ring = r <= SPOT_RADIUS_PX, (r >= RING_PX[0]) & (r <= RING_PX[1])
    if ap.sum() == 0 or ring.sum() < 20:
        raise ValueError(f"spot {spot} is too close to / outside the ROI {run['meta']['roi']}")
    return ap, ring


def spot_contrast(run, spot=None, sigma=2.0):
    """How compact and bright the spot is: peak of the 2 px-smoothed mean image inside the
    aperture, above the ring median, in units of the ring's robust spread. A real beam gave ~40
    (20 Oct); stray light / blocked arm gave ~3-4. Needs no background run."""
    ap, ring = spot_masks(run, spot)
    sm = gaussian_filter(crop(run, run["sum"]) / n_rd(run), sigma)
    r = sm[ring]
    spread = 1.4826 * np.median(np.abs(r - np.median(r)))
    return float((sm[ap].max() - np.median(r)) / spread) if spread > 0 else np.inf


def is_thresholded(run):
    return run["meta"]["actual"]["photon_mode"] != "Disabled"


def spot_series(run, spot=None):
    """Per readout: sum over the spot aperture minus (aperture pixels x ring median), ADU.
    For camera-thresholded runs (0/1 pixels) this is net hits in the spot."""
    ap, ring = spot_masks(run, spot)
    st = run["roi_stack"].astype(np.float64)
    if is_thresholded(run):
        return np.array([f[ap].sum() - ap.sum() * f[ring].mean() for f in st])
    return np.array([f[ap].sum() - ap.sum() * np.median(f[ring]) for f in st])


def above_pedestal(run):
    """Per-readout ROI crops minus each readout's own median (ADU above pedestal)."""
    st = run["roi_stack"].astype(np.float64)
    if is_thresholded(run) or not st.size:
        return st
    return st - np.median(st.reshape(len(st), -1), axis=1)[:, None, None]


def hit_series(run, threshold, mask=None):
    """Per readout: pixels at or above `threshold` ADU above that readout's pedestal
    (inside `mask` if given). Camera-thresholded runs: count of 1-pixels."""
    st = above_pedestal(run)
    hit = st >= 0.5 if is_thresholded(run) else st >= threshold
    if mask is not None:
        hit = hit[:, mask]
        return hit.sum(axis=1).astype(float)
    return hit.sum(axis=(1, 2)).astype(float)


def _diff_z(x, y):
    """(mean(x) - mean(y)) / combined standard error; 0 if both are identically zero."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 2 or len(y) < 2:
        return np.nan
    err = np.hypot(x.std(ddof=1) / np.sqrt(len(x)), y.std(ddof=1) / np.sqrt(len(y)))
    d = x.mean() - y.mean()
    if err == 0:
        return 0.0 if d == 0 else np.sign(d) * np.inf
    return float(d / err)


def spot_metrics(run, dark=None, spot=None, threshold=None):
    """Spot numbers per readout, and their significance against the dark if given."""
    out = {}
    try:
        ap, _ = spot_masks(run, spot)
    except ValueError as e:
        print(f"  ({e})")
        return out
    s = spot_series(run, spot)
    out["spot_adu_per_readout"] = float(s.mean())
    if threshold is not None:
        h = hit_series(run, threshold, ap)
        out["spot_hits_per_readout"] = float(h.mean())
    if dark is not None and len(dark["roi_stack"]):
        d = spot_series(dark, spot)
        out["spot_net_adu_per_readout"] = float(s.mean() - d.mean())
        out["spot_z"] = _diff_z(s, d)
        if threshold is not None:
            hd = hit_series(dark, threshold, ap)
            out["spot_dark_hits_per_readout"] = float(hd.mean())
            out["spot_hits_z"] = _diff_z(h, hd)
    return out


def run_metrics(run, dark=None, spot=None, threshold=None):
    m, a, s = run["meta"], run["meta"]["actual"], run["meta"]["settings"]
    roi_sum = crop(run, run["sum"])
    dt = readout_interval_s(run)
    out = dict(
        label=s["label"], kind=s["kind"], started=m["started"],
        gain=a["gain"], gate_delay_ns=a["gate_delay_ns"], gate_width_ns=a["gate_width_ns"],
        accumulations=a["accumulations"], n_readouts=n_rd(run), gates_total=gates(run),
        gate_open_ms=a["gate_width_ns"] * gates(run) * 1e-6,
        wall_s=m["wall_s"],
        trigger=(f"Ext {a['trigger_threshold_v']} V" if a["trigger_source"] == "External"
                 else f"Int {float(a['internal_trigger_hz']) / 1e3:g} kHz"),
        photon_mode=a["photon_mode"],
        photon_threshold=a["photon_threshold"] if a["photon_mode"] != "Disabled" else None,
        temp_c=a["temp_c"], idler_khz=s["idler_rate_khz"],
        roi=f"{m['roi'][2]} px @ ({m['roi'][0] + m['roi'][2] // 2},{m['roi'][1] + m['roi'][2] // 2})",
        full_mean_per_readout=float(run["sum"].mean() / n_rd(run)),
        roi_mean_per_readout=float(roi_sum.mean() / n_rd(run)),
        roi_total=float(roi_sum.sum()),
        saturated_frac=float((run["first"] >= SATURATION_ADU).mean()),
        s_per_readout=dt,
        min_gates_per_s=float(int(a["accumulations"]) / dt) if dt and np.isfinite(dt) else np.nan,
    )
    out.update(spot_metrics(run, dark, spot, threshold))
    if dark is not None:
        dark_roi = crop(run, scaled_dark_sum(run, dark))
        net = roi_sum - dark_roi
        out.update(
            dark_label=dark["meta"]["settings"]["label"],
            net_roi_total=float(net.sum()),
            net_per_1k_gates=float(net.sum() / (gates(run) / 1000)),
            z_net=net_significance(run, dark),
            snr_cameron=float(net.mean() / dark_roi.std()) if dark_roi.std() > 0 else np.nan,
            snr_xu=float((roi_sum.mean() - dark_roi.mean()) / np.sqrt(dark_roi.sum()))
            if dark_roi.sum() > 0 else np.nan,
            mismatched=", ".join(match_warnings(run, dark)),
        )
    return out


def settings_title(run):
    a, s, m = run["meta"]["actual"], run["meta"]["settings"], run["meta"]
    trig = (f"ext trig {a['trigger_threshold_v']:g} V" if a["trigger_source"] == "External"
            else f"int trig {float(a['internal_trigger_hz']) / 1e3:g} kHz")
    pd = "PD off" if a["photon_mode"] == "Disabled" else f"PD {a['photon_mode']} {a['photon_threshold']}"
    temp = f"{a['temp_c']:.0f} C" if isinstance(a["temp_c"], (int, float)) else "? C"
    return (f"{s['label']} [{s['kind']}]  |  gain {a['gain']}  |  gate {a['gate_delay_ns']:g} ns + "
            f"{fmt_ns(a['gate_width_ns'])}  |  {int(a['accumulations']):,} acc x {n_rd(run)} readouts  |  "
            f"{trig}  |  {pd}  |  {temp}")


def find_peak(img, sigma=3.0):
    """(row, col) of the brightest spot after Gaussian smoothing (ignores single hot pixels)."""
    sm = gaussian_filter(np.asarray(img, dtype=np.float64), sigma)
    return np.unravel_index(np.argmax(sm), sm.shape)


def _dark_pixel_noise(run, dark):
    """Per-pixel noise (std) of the net summed image, from the dark's readout-to-readout scatter."""
    st = dark["roi_stack"].astype(np.float64)
    if len(st) < 2:
        return None
    sig_px = float(np.median(st.std(axis=0, ddof=1)))
    ns, nd = n_rd(run), n_rd(dark)
    return sig_px * np.sqrt(ns + ns ** 2 / nd)


def plot_run(run, dark=None, threshold=None, spot=None, band=5, save_as=None):
    """Six-panel figure for one run:
    (a) full frame, (b) ROI zoom, (c) ROI minus dark with the spot aperture,
    (d) histogram of pixel values above each readout's pedestal, signal vs dark,
    (e) row and (f) column profile through the spot.
    `threshold` is in ADU above pedestal (from threshold_sweep); `spot` is (x, y) on the
    full sensor (default: brightest spot of panel c)."""
    n = n_rd(run)
    full = run["sum"] / n
    roi_raw = crop(run, run["sum"]) / n
    net = crop(run, net_full(run, dark))
    x0, y0, size = run["meta"]["roi"]
    if spot is None:
        spot = find_spot(run, dark)
    px, py = int(spot[0]) - x0, int(spot[1]) - y0

    fig, ax = plt.subplots(2, 3, figsize=(17, 10))
    # (a) full frame
    im = ax[0, 0].imshow(full, cmap="inferno", vmin=np.percentile(full, 1), vmax=np.percentile(full, 99.7))
    ax[0, 0].add_patch(Rectangle((x0, y0), size, size, fill=False, ec="cyan", lw=1.5))
    ax[0, 0].set_title("(a) full frame, mean per readout (cyan = ROI)")
    ax[0, 0].set_xlabel("x (px)"); ax[0, 0].set_ylabel("y (px)")
    fig.colorbar(im, ax=ax[0, 0], label="ADU per readout")
    # (b) ROI zoom
    im = ax[0, 1].imshow(roi_raw, cmap="inferno", vmin=np.percentile(roi_raw, 1),
                         vmax=np.percentile(roi_raw, 99.7), extent=(x0, x0 + size, y0 + size, y0))
    ax[0, 1].set_title(f"(b) ROI {size}x{size} px, raw, mean per readout")
    ax[0, 1].set_xlabel("x (px)"); ax[0, 1].set_ylabel("y (px)")
    fig.colorbar(im, ax=ax[0, 1], label="ADU per readout")
    # (c) dark-subtracted ROI with the spot aperture and ring
    lim = np.percentile(np.abs(net), 99.5) or 1
    im = ax[0, 2].imshow(net, cmap="RdBu_r", vmin=-lim, vmax=lim, extent=(x0, x0 + size, y0 + size, y0))
    for rad, ls in ((SPOT_RADIUS_PX, "-"), (RING_PX[0], ":"), (RING_PX[1], ":")):
        ax[0, 2].add_patch(Circle((x0 + px, y0 + py), rad, fill=False, ec="k", lw=.8, ls=ls))
    bgname = ("off-delay background" if dark is not None and dark["meta"]["settings"]["kind"] == "dark_offdelay"
              else "dark")
    sub = (f"minus {bgname} '{dark['meta']['settings']['label']}'" if dark is not None
           else "minus frame median (no background run)")
    ax[0, 2].set_title(f"(c) ROI {sub}\nsum of {n} readouts; circle = spot, dotted = ring")
    ax[0, 2].set_xlabel("x (px)"); ax[0, 2].set_ylabel("y (px)")
    fig.colorbar(im, ax=ax[0, 2], label="net ADU (summed)")
    # (d) histogram of per-readout pixel values above pedestal
    if is_thresholded(run):
        ax[1, 0].text(.5, .5, "camera photon thresholding on:\npixels are 0/1, no histogram",
                      ha="center", va="center", transform=ax[1, 0].transAxes)
    else:
        sv = above_pedestal(run).ravel()
        dv = above_pedestal(dark).ravel() if dark is not None and not is_thresholded(dark) else None
        allv = sv if dv is None else np.concatenate([sv, dv])
        if allv.size:
            lo, hi = np.percentile(allv, 0.001), np.percentile(allv, 99.999)
            bins = np.linspace(lo, max(hi, lo + 1), 200)
            ax[1, 0].hist(sv, bins=bins, histtype="step", lw=1.5, color="C3",
                          weights=np.full(sv.size, 1 / max(len(run["roi_stack"]), 1)),
                          label=f"signal ({len(run['roi_stack'])} readouts)")
            if dv is not None:
                ax[1, 0].hist(dv, bins=bins, histtype="step", lw=1.5, color="k",
                              weights=np.full(dv.size, 1 / max(len(dark["roi_stack"]), 1)),
                              label=f"{bgname} ({len(dark['roi_stack'])} readouts)")
        if threshold is not None:
            ax[1, 0].axvline(threshold, color="C2", ls="--", label=f"threshold {threshold:g}")
        ax[1, 0].set_yscale("log")
        ax[1, 0].legend(fontsize=8)
    ax[1, 0].set_xlabel("ROI pixel value in ONE readout, minus that readout's median (ADU)")
    ax[1, 0].set_ylabel("pixels per readout")
    ax[1, 0].set_title("(d) pixel values above pedestal: photon hits = excess tail")
    # (e, f) profiles through the spot
    noise = _dark_pixel_noise(run, dark) if dark is not None else None
    r0, r1 = max(py - band // 2, 0), min(py + band // 2 + 1, size)
    c0, c1 = max(px - band // 2, 0), min(px + band // 2 + 1, size)
    row = net[r0:r1, :].mean(axis=0)
    col = net[:, c0:c1].mean(axis=1)
    for a_, prof, coord, which, pos in ((ax[1, 1], row, np.arange(x0, x0 + size), "row", y0 + py),
                                        (ax[1, 2], col, np.arange(y0, y0 + size), "column", x0 + px)):
        a_.plot(coord, prof, color="C3", lw=1.2,
                label="signal minus dark" if dark is not None else "signal minus median")
        a_.axhline(0, color="k", lw=.6)
        if noise is not None:
            s2 = 2 * noise / np.sqrt(r1 - r0 if which == "row" else c1 - c0)
            a_.axhspan(-s2, s2, color="grey", alpha=.25, label="+/-2 sigma dark noise")
        a_.set_xlabel(("x" if which == "row" else "y") + " (px)")
        a_.set_ylabel(f"net ADU (mean of {band} {which}s)")
        a_.set_title(f"({'e' if which == 'row' else 'f'}) {which} profile through "
                     f"{'y' if which == 'row' else 'x'} = {pos}")
        a_.legend(fontsize=8)
    title = settings_title(run) + f"  |  spot at ({x0 + px}, {y0 + py})"
    m = run_metrics(run, dark, spot=(x0 + px, y0 + py), threshold=threshold)
    if dark is not None and "spot_z" in m:
        title += (f"\nspot (r = {SPOT_RADIUS_PX} px, ring-subtracted): {m['spot_net_adu_per_readout']:.3g} "
                  f"net ADU per readout, {m['spot_z']:.1f} sigma above {bgname}")
        if "spot_hits_z" in m:
            title += (f"  |  hits >= {threshold:g} ADU: {m['spot_hits_per_readout']:.2f} vs {bgname} "
                      f"{m['spot_dark_hits_per_readout']:.2f} per readout ({m['spot_hits_z']:.1f} sigma)")
        warn = match_warnings(run, dark)
        if warn:
            title += "\nSETTINGS MISMATCH vs dark: " + "; ".join(warn)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    if save_as:
        _savefig(fig, save_as)
    plt.show()
    return fig


def plot_full_frame(run, dark, save_as=None):
    """Full sensor: signal, dark, and signal minus dark on matching scales. Use it to see
    whether features outside the ROI (other spots, glow, streaks) come from the signal arm
    (they vanish in the dark) or not (they are in both)."""
    sig = run["sum"] / n_rd(run)
    drk = dark["sum"] / n_rd(dark)
    net = gaussian_filter(sig - drk, 2)
    vmin, vmax = np.percentile(np.r_[sig.ravel(), drk.ravel()], [1, 99.7])
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
    x0, y0, size = run["meta"]["roi"]
    for a_, img, t in ((ax[0], sig, f"signal '{run['meta']['settings']['label']}'"),
                       (ax[1], drk, f"dark '{dark['meta']['settings']['label']}'")):
        im = a_.imshow(img, cmap="inferno", vmin=vmin, vmax=vmax)
        a_.set_title(t + ", mean per readout")
        fig.colorbar(im, ax=a_, label="ADU per readout")
    lim = np.percentile(np.abs(net), 99.8) or 1
    im = ax[2].imshow(net, cmap="RdBu_r", vmin=-lim, vmax=lim)
    ax[2].set_title("signal minus dark (2 px smoothing): red = from the signal arm")
    fig.colorbar(im, ax=ax[2], label="ADU per readout")
    for a_ in ax:
        a_.add_patch(Rectangle((x0, y0), size, size, fill=False, ec="cyan", lw=1))
        a_.set_xlabel("x (px)"); a_.set_ylabel("y (px)")
    fig.suptitle(settings_title(run), fontsize=10)
    fig.tight_layout()
    if save_as:
        _savefig(fig, save_as)
    plt.show()
    return fig


def plot_raw_vs_accumulated(run, dark=None, save_as=None):
    """One raw readout next to the full accumulation (ROI), like the thesis raw/accumulated
    comparison: single readouts show isolated hits, the sum builds up the beam profile."""
    first = crop(run, run["first"])
    tot = crop(run, net_full(run, dark))
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.8))
    im = ax[0].imshow(first, cmap="inferno", vmin=np.percentile(first, 1), vmax=np.percentile(first, 99.9))
    ax[0].set_title("one raw readout (ROI)")
    fig.colorbar(im, ax=ax[0], label="ADU")
    im = ax[1].imshow(gaussian_filter(tot, 1), cmap="inferno", vmin=np.percentile(tot, 1),
                      vmax=np.percentile(gaussian_filter(tot, 1), 99.8))
    ax[1].set_title(f"sum of {n_rd(run)} readouts " + ("minus dark" if dark is not None else "minus median")
                    + " (1 px smoothing)")
    fig.colorbar(im, ax=ax[1], label="ADU")
    for a_ in ax:
        a_.set_xlabel("x in ROI (px)"); a_.set_ylabel("y in ROI (px)")
    fig.suptitle(settings_title(run), fontsize=10)
    fig.tight_layout()
    if save_as:
        _savefig(fig, save_as)
    plt.show()
    return fig


def suggest_roi_center(run, dark=None, sigma=15, edge=40):
    """Brightest smoothed spot of the full-frame net image (edges ignored). Plots it."""
    img = gaussian_filter(net_full(run, dark), sigma)
    img[:edge, :] = img[-edge:, :] = img[:, :edge] = img[:, -edge:] = -np.inf
    cy, cx = np.unravel_index(np.argmax(img), img.shape)
    show = np.where(np.isfinite(img), img, np.nan)
    plt.figure(figsize=(5.5, 5))
    plt.imshow(show, cmap="inferno")
    plt.plot(cx, cy, "c+", ms=20, mew=2)
    plt.title(f"smoothed net full frame (sigma {sigma} px)\nsuggested ROI centre ({cx}, {cy})")
    plt.xlabel("x (px)"); plt.ylabel("y (px)"); plt.colorbar(label="ADU")
    plt.tight_layout(); plt.show()
    return int(cx), int(cy)


# --------------------------------------------------------------------------------------
# Comparing runs
# --------------------------------------------------------------------------------------

def auto_pair(runs):
    """[(signal_run, dark_or_None), ...]: each non-dark run gets the most recent dark run
    taken before it with the same gain, gate width, accumulations, photon settings, ROI."""
    darks = [r for r in runs if r["meta"]["settings"]["kind"].startswith("dark")]
    out = []
    for r in runs:
        if r["meta"]["settings"]["kind"].startswith("dark"):
            continue
        cands = [d for d in darks if not match_warnings(r, d)]
        before = [d for d in cands if d["meta"]["started"] <= r["meta"]["started"]]
        pick = (before or cands or [None])[-1]
        out.append((r, pick))
    return out


def summary_table(pairs, spot=None, threshold=None):
    """pandas DataFrame, one row per run: settings next to mean / total counts, plus spot
    numbers (aperture around `spot`, default ROI centre) and photon hits at `threshold`
    (ADU above pedestal)."""
    import pandas as pd
    rows = [run_metrics(r, d, spot=spot, threshold=threshold) for r, d in pairs]
    df = pd.DataFrame(rows)
    nice = {
        "label": "label", "kind": "kind", "gain": "gain", "gate_delay_ns": "delay (ns)",
        "gate_width_ns": "width (ns)", "accumulations": "accum", "n_readouts": "readouts",
        "trigger": "trigger", "photon_mode": "photon det", "photon_threshold": "PD thr",
        "temp_c": "T (C)", "idler_khz": "idler (kHz)",
        "s_per_readout": "s / readout", "min_gates_per_s": "gates/s (>=)",
        "roi_mean_per_readout": "ROI mean (ADU/px/readout)", "roi_total": "ROI total (ADU)",
        "spot_adu_per_readout": "spot ADU/readout",
        "spot_net_adu_per_readout": "spot net ADU/readout", "spot_z": "spot sigma",
        "spot_hits_per_readout": "spot hits/readout",
        "spot_dark_hits_per_readout": "dark spot hits/readout", "spot_hits_z": "hits sigma",
        "net_per_1k_gates": "ROI net per 1k gates", "snr_cameron": "SNR Cameron (ROI)",
        "dark_label": "dark used", "saturated_frac": "saturated frac", "wall_s": "time (s)",
        "mismatched": "dark mismatch",
    }
    cols = [c for c in nice if c in df.columns]
    return df[cols].rename(columns=nice)


def plot_compare(pairs, ncols=4, spot=None, save_as=None):
    """Net ROI images side by side, normalised per 1000 gates, on ONE shared colour scale,
    so brightness is directly comparable between runs."""
    imgs = [crop(r, net_full(r, d)) / (gates(r) / 1000) for r, d in pairs]
    lim = max(np.percentile(np.abs(i), 99.5) for i in imgs) or 1
    nrows = int(np.ceil(len(pairs) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 4.4 * nrows), squeeze=False)
    for a_ in axes.ravel():
        a_.axis("off")
    for a_, (r, d), img in zip(axes.ravel(), pairs, imgs):
        a = r["meta"]["actual"]
        a_.imshow(img, cmap="RdBu_r", vmin=-lim, vmax=lim)
        if d is not None:
            z = spot_metrics(r, d, spot).get("spot_z", np.nan)
            ztxt = f", spot {z:.1f} sigma"
        else:
            ztxt = ", NO dark"
        a_.set_title(f"{r['meta']['settings']['label']}\ngain {a['gain']}, gate {a['gate_delay_ns']:g}+"
                     f"{fmt_ns(a['gate_width_ns'])}\n{int(a['accumulations']):,} acc x {n_rd(r)}{ztxt}",
                     fontsize=9)
    sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(-lim, lim))
    fig.colorbar(sm, ax=axes.ravel().tolist(), shrink=0.6, label="net ADU per 1000 gates")
    fig.suptitle("Net ROI (signal minus dark), same colour scale for every panel")
    if save_as:
        _savefig(fig, save_as)
    plt.show()
    return fig


# --------------------------------------------------------------------------------------
# Photon threshold (all thresholds are ADU ABOVE each readout's own pedestal)
# --------------------------------------------------------------------------------------

def dark_level(dark):
    """(median, robust std) of single-readout ROI pixel values in a dark run (absolute ADU)."""
    v = dark["roi_stack"].astype(np.float64).ravel()
    med = float(np.median(v))
    return med, float(1.4826 * np.median(np.abs(v - med)))


def threshold_sweep(run, dark, thresholds=None, spot=None, save_as=None):
    """Software photon counting on the saved per-readout ROI crops (runs must have photon
    detection 'Disabled'). Thresholds are ADU above each readout's own median, so the
    answer carries over to runs whose pedestal drifted (e.g. different gate widths).
    Counts hits in the spot aperture if `spot` is given, else in the whole ROI.
    Returns (DataFrame, best_threshold_above_pedestal)."""
    import pandas as pd
    if is_thresholded(run) or is_thresholded(dark):
        raise ValueError("threshold_sweep needs runs taken with photon_mode='Disabled'")
    S, D = above_pedestal(run), above_pedestal(dark)
    med_abs, sd = dark_level(dark)
    mask = spot_masks(run, spot)[0] if spot is not None else None
    if thresholds is None:
        top = np.percentile(S, 99.995)
        thresholds = np.unique(np.round(np.linspace(3 * max(sd, 1), max(top, 4 * sd), 60)))
    rows = []
    for t in thresholds:
        s = (S[:, mask] >= t).sum(axis=1) if mask is not None else (S >= t).sum(axis=(1, 2))
        d = (D[:, mask] >= t).sum(axis=1) if mask is not None else (D >= t).sum(axis=(1, 2))
        rows.append(dict(threshold=t, sigmas_above_dark=t / sd if sd else np.nan,
                         camera_threshold_adu=med_abs + t,
                         signal_hits=s.mean(), dark_hits=d.mean(), net_hits=s.mean() - d.mean(),
                         z=_diff_z(s, d)))
    df = pd.DataFrame(rows)
    zf = df["z"].replace([np.inf, -np.inf], np.nan)
    best = float(df.loc[zf.idxmax(), "threshold"]) if zf.notna().any() else float(thresholds[0])
    where = "spot aperture" if mask is not None else "whole ROI"

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
    ax[0].semilogy(df.threshold, df.signal_hits, "C3.-", label="signal")
    ax[0].semilogy(df.threshold, df.dark_hits, "k.-", label="dark")
    ax[0].set_xlabel("threshold (ADU above pedestal)")
    ax[0].set_ylabel(f"pixels >= threshold per readout ({where})")
    ax[0].set_title("hits per readout vs threshold"); ax[0].legend()
    ax[1].plot(df.threshold, zf, "C0.-")
    ax[1].axvline(best, color="C2", ls="--",
                  label=f"best {best:g} ADU above pedestal ({best / sd:.1f} dark sigma)")
    ax[1].axhline(0, color="k", lw=.6)
    ax[1].set_xlabel("threshold (ADU above pedestal)"); ax[1].set_ylabel("(signal - dark) hits / error")
    ax[1].set_title("significance of the excess"); ax[1].legend(fontsize=8)
    hit_s = (S >= best).sum(axis=0).astype(float)
    hit_d = (D >= best).sum(axis=0).astype(float) * len(S) / max(len(D), 1)
    netmap = gaussian_filter(hit_s - hit_d, 1.5)
    lim = np.percentile(np.abs(netmap), 99.5) or 1
    im = ax[2].imshow(netmap, cmap="RdBu_r", vmin=-lim, vmax=lim)
    ax[2].set_title(f"net hit map at {best:g} ADU above pedestal ({len(S)} readouts)")
    ax[2].set_xlabel("x in ROI (px)"); ax[2].set_ylabel("y in ROI (px)")
    fig.colorbar(im, ax=ax[2], label="net hits (1.5 px smoothing)")
    fig.suptitle(settings_title(run) + f"\ndark '{dark['meta']['settings']['label']}': pedestal "
                 f"{med_abs:.0f} ADU, robust sigma {sd:.1f} ADU  ->  camera Photon Detection "
                 f"Threshold equivalent ~ {med_abs + best:.0f} ADU", fontsize=10)
    fig.tight_layout()
    if save_as:
        _savefig(fig, save_as)
    plt.show()
    return df, best


# --------------------------------------------------------------------------------------
# Timing: expected delay, expected counts, halving ladder, delay scan
# --------------------------------------------------------------------------------------

def expected_delay_ns(signal_fiber_m=(110, 130), idler_fiber_m=0.0, n_group=1.468,
                      signal_free_space_m=0.0, idler_free_space_m=0.0, apd_latency_ns=(15, 30),
                      trigger_cable_m=(1, 5), cable_ns_per_m=5.0, camera_insertion_ns=35.0):
    """Rough window for the Repetitive Gate delay. Returns (low, high) in ns.
    gate delay = signal arrival - (idler arrival + APD latency + cable + camera insertion delay)."""
    c = 0.299792458                   # m per ns
    def rng(x):
        return (x, x) if np.isscalar(x) else tuple(x)
    sf, ifb, apd, cab = rng(signal_fiber_m), rng(idler_fiber_m), rng(apd_latency_ns), rng(trigger_cable_m)
    sig = lambda L: L * n_group / c + signal_free_space_m / c
    idl = ifb[0] * n_group / c + idler_free_space_m / c
    lo = sig(sf[0]) - idl - apd[1] - cab[1] * cable_ns_per_m - camera_insertion_ns
    hi = sig(sf[1]) - idl - apd[0] - cab[0] * cable_ns_per_m - camera_insertion_ns
    return float(lo), float(hi)


def expected_counts(idler_khz=15.0, signal_khz=18.0, coinc_khz=1.5, signal_apd_eff=0.6,
                    transmission_to_camera=1.0, qe=0.25, accumulations=5000,
                    widths_ns=(100000, 10000, 1000, 100, 30, 10, 3)):
    """Back-of-envelope photoelectrons per readout, heralded vs unheralded, vs gate width.
    Assumes the gate delay is correct. Heralded photons are only those correlated with the
    trigger; unheralded = every other signal-arm photon that happens to fall in the gate."""
    import pandas as pd
    p_herald = (coinc_khz / idler_khz) / signal_apd_eff * transmission_to_camera * qe
    flux_hz = signal_khz * 1e3 / signal_apd_eff * transmission_to_camera
    rows = []
    for w in widths_ns:
        unh = flux_hz * w * 1e-9 * qe
        rows.append({"gate width": fmt_ns(w),
                     "heralded per readout": p_herald * accumulations,
                     "unheralded per readout": unh * accumulations,
                     "heralded fraction": p_herald / (p_herald + unh)})
    return pd.DataFrame(rows)


def take_offdelay_background(cam, settings, shift_ns=5000.0, n_readouts=None, label=None):
    """Background that needs nobody in the lab: identical settings, but the gate is moved
    `shift_ns` later, away from the heralded photon. It still sees the unheralded light,
    dark counts, pedestal and hot pixels, so (signal - this) is the HERALDED part only.
    `settings` is a settings dict or a run (its settings are reused). Only valid when the
    gate width is well below the shift; a 100 us gate cannot have an off-delay background."""
    s = dict(settings["meta"]["settings"]) if isinstance(settings, dict) and "meta" in settings else dict(settings)
    if float(s["gate_width_ns"]) >= shift_ns:
        raise ValueError(f"gate width {fmt_ns(s['gate_width_ns'])} is not much shorter than the "
                         f"shift {fmt_ns(shift_ns)}; the background window would overlap the signal window")
    s.update(kind="dark_offdelay", gate_delay_ns=float(s["gate_delay_ns"]) + shift_ns,
             label=label or f"{s['label']}_offdelay")
    if n_readouts is not None:
        s["n_readouts"] = int(n_readouts)
    return take_run(cam, **{k: v for k, v in s.items() if k in DEFAULT_SETTINGS})


def take_halves(cam, base, delay_ns, width_ns, label="ladder", **overrides):
    """Two runs that split the gate window [delay, delay+width] into an early and a late half.
    Both halves have the same width, so unheralded light, dark counts and pedestal are the
    same in both; any difference is heralded signal."""
    half = width_ns / 2
    common = {**base, **overrides}
    tag = fmt_ns(width_ns).replace(" ", "")
    a = take_run(cam, **{**common, "label": f"{label}_{tag}_early", "kind": "signal",
                         "gate_delay_ns": delay_ns, "gate_width_ns": half})
    b = take_run(cam, **{**common, "label": f"{label}_{tag}_late", "kind": "signal",
                         "gate_delay_ns": delay_ns + half, "gate_width_ns": half})
    return a, b


def compare_halves(a, b, z_min=3.0, threshold=None, spot=None, save_as=None):
    """Early-half vs late-half comparison, scored on the spot aperture (ring-subtracted, so
    pedestal drift cancels) and, if `threshold` is given, on photon hits in the spot.
    Returns a dict with the recommended next window."""
    ap, _ = spot_masks(a, spot)
    sa, sb = spot_series(a, spot), spot_series(b, spot)
    z = _diff_z(sa, sb)
    zt, ha, hb = None, None, None
    if threshold is not None:
        ha, hb = hit_series(a, threshold, ap), hit_series(b, threshold, ap)
        zt = _diff_z(ha, hb)
    decide = zt if (zt is not None and np.isfinite(zt) and (ha.sum() + hb.sum()) >= 10) else z
    aa, ab = a["meta"]["actual"], b["meta"]["actual"]
    if decide >= z_min:
        pick, nxt = "early", (aa["gate_delay_ns"], aa["gate_width_ns"])
    elif decide <= -z_min:
        pick, nxt = "late", (ab["gate_delay_ns"], ab["gate_width_ns"])
    else:
        pick, nxt = "unclear", None

    ia = crop(a, a["sum"]) / n_rd(a)
    ib = crop(b, b["sum"]) / n_rd(b)
    ia_s, ib_s = ia - np.median(ia), ib - np.median(ib)      # remove pedestal drift for display
    d = gaussian_filter(ia_s - ib_s, 1.5)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.8))
    vmin, vmax = np.percentile(np.r_[ia_s.ravel(), ib_s.ravel()], [1, 99.7])
    x0, y0, _ = a["meta"]["roi"]
    cx, cy = (spot if spot is not None else roi_center_xy(a))
    for a_, img, run, name in ((ax[0], ia_s, a, "early"), (ax[1], ib_s, b, "late")):
        im = a_.imshow(img, cmap="inferno", vmin=vmin, vmax=vmax)
        ra = run["meta"]["actual"]
        a_.set_title(f"{name} half: {ra['gate_delay_ns']:g} ns + {fmt_ns(ra['gate_width_ns'])}")
        fig.colorbar(im, ax=a_, label="ADU per readout above median")
    lim = np.percentile(np.abs(d), 99.5) or 1
    im = ax[2].imshow(d, cmap="RdBu_r", vmin=-lim, vmax=lim)
    ax[2].set_title("early minus late (red = early brighter)")
    fig.colorbar(im, ax=ax[2], label="ADU per readout")
    for a_ in ax:
        a_.add_patch(Circle((cx - x0, cy - y0), SPOT_RADIUS_PX, fill=False, ec="c", lw=1))
        a_.set_xlabel("x in ROI (px)"); a_.set_ylabel("y in ROI (px)")
    ztxt = f"spot ADU: {sa.mean():.0f} vs {sb.mean():.0f} per readout, z = {z:+.1f}"
    if zt is not None:
        ztxt += (f"  |  spot hits >= {threshold:g}: {ha.mean():.2f} vs {hb.mean():.2f} per readout, "
                 f"z = {zt:+.1f}")
        if ha.sum() + hb.sum() == 0:
            ztxt += " (NO hits in either half)"
    verdict = (f"keep {pick}: next window {nxt[0]:g} ns + {fmt_ns(nxt[1])}" if nxt else
               f"unclear (|z| < {z_min:g}): no heralded excess in either half, or it straddles the boundary")
    fig.suptitle(f"Halving step {aa['gate_delay_ns']:g} ns + {fmt_ns(2 * aa['gate_width_ns'])}: {verdict}\n"
                 f"{ztxt}\ngain {aa['gain']}, {int(aa['accumulations']):,} acc x {n_rd(a)}/{n_rd(b)} readouts",
                 fontsize=10)
    fig.tight_layout()
    if save_as:
        _savefig(fig, save_as)
    plt.show()
    return dict(z=z, z_hits=zt, pick=pick, next_delay_ns=nxt[0] if nxt else None,
                next_width_ns=nxt[1] if nxt else None)


def delay_scan(cam, base, delays_ns, width_ns, label="dscan", **overrides):
    """One run per delay at a fixed gate width. Returns the list of runs."""
    runs = []
    for d in delays_ns:
        r = take_run(cam, **{**base, **overrides, "label": f"{label}_{d:g}ns", "kind": "signal",
                             "gate_delay_ns": float(d), "gate_width_ns": float(width_ns)})
        if r is not None:
            runs.append(r)
    return runs


def _mean_err(x):
    x = np.asarray(x, float)
    return x.mean(), (x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else 0.0)


def plot_delay_scan(runs, threshold=None, spot=None, save_as=None):
    """Spot ADU per readout (ring-subtracted) and, optionally, photon hits in the spot vs
    gate delay, with standard errors; then the ROI image at each delay on one colour scale.
    Heralded photons show up as a peak at one delay; unheralded light is flat."""
    runs = sorted(runs, key=lambda r: r["meta"]["actual"]["gate_delay_ns"])
    delays = np.array([r["meta"]["actual"]["gate_delay_ns"] for r in runs])
    ap = spot_masks(runs[0], spot)[0]
    adu = [_mean_err(spot_series(r, spot)) for r in runs]
    ncol = 2 if threshold is not None else 1
    fig, ax = plt.subplots(1, ncol, figsize=(7 * ncol, 4.2), squeeze=False)
    ax[0, 0].errorbar(delays, [m for m, _ in adu], [e for _, e in adu], fmt="o-", capsize=3)
    ax[0, 0].set_xlabel("gate delay (ns)"); ax[0, 0].set_ylabel("spot ADU per readout (ring-subtracted)")
    ax[0, 0].set_title("spot signal vs delay")
    if threshold is not None:
        hits = [_mean_err(hit_series(r, threshold, ap)) for r in runs]
        ax[0, 1].errorbar(delays, [m for m, _ in hits], [e for _, e in hits], fmt="o-", capsize=3, color="C3")
        ax[0, 1].set_xlabel("gate delay (ns)")
        ax[0, 1].set_ylabel(f"spot pixels >= {threshold:g} ADU above pedestal, per readout")
        ax[0, 1].set_title("photon hits in the spot vs delay")
    a = runs[0]["meta"]["actual"]
    fig.suptitle(f"Delay scan: gate width {fmt_ns(a['gate_width_ns'])}, gain {a['gain']}, "
                 f"{int(a['accumulations']):,} acc x {n_rd(runs[0])} readouts per point", fontsize=10)
    fig.tight_layout()
    if save_as:
        _savefig(fig, save_as + "_curve")
    plt.show()

    imgs = []
    for r in runs:
        i = crop(r, r["sum"]) / n_rd(r)
        imgs.append(gaussian_filter(i - np.median(i), 1.5))
    lim = max(np.percentile(np.abs(i), 99.5) for i in imgs) or 1
    nc = min(len(runs), 6)
    nr = int(np.ceil(len(runs) / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(2.8 * nc, 3 * nr), squeeze=False)
    for a_ in axes.ravel():
        a_.axis("off")
    for a_, img, d in zip(axes.ravel(), imgs, delays):
        a_.imshow(img, cmap="RdBu_r", vmin=-lim, vmax=lim)
        a_.set_title(f"{d:g} ns", fontsize=9)
    fig.suptitle("ROI per readout minus its own median, same colour scale (1.5 px smoothing)")
    fig.tight_layout()
    if save_as:
        _savefig(fig, save_as + "_images")
    plt.show()


def width_series(cam, base, widths_ns, delay_ns=None, label="wseries", **overrides):
    """One run per gate width at a fixed delay. Returns the list of runs."""
    runs = []
    kind = overrides.pop("kind", "signal")
    for w in widths_ns:
        kw = {**base, **overrides, "label": f"{label}_{fmt_ns(w).replace(' ', '')}", "kind": kind,
              "gate_width_ns": float(w)}
        if delay_ns is not None:
            kw["gate_delay_ns"] = float(delay_ns)
        r = take_run(cam, **kw)
        if r is not None:
            runs.append(r)
    return runs


def measure_trigger_rate(cam, base, accumulations=100_000, width_ns=3.0, n_readouts=4, **overrides):
    """Estimate the trigger rate from the camera alone. With a tiny gate and many
    accumulations, each readout has to wait for `accumulations` triggers, so
    rate ~ accumulations / (time per readout). If the readout time dominates instead,
    only a lower bound is possible. Compare the answer with the idler counter."""
    r = take_run(cam, **{**base, **overrides, "label": "trigger_rate_test", "kind": "signal",
                         "gate_width_ns": float(width_ns), "accumulations": int(accumulations),
                         "n_readouts": int(n_readouts) + 1})
    if r is None:
        return None
    dt = readout_interval_s(r)
    ro = r["meta"]["readback"].get("Readout Time Calculation")
    ro_s = float(ro) / 1000 if isinstance(ro, (int, float)) and ro > 0 else 0.23
    rate = int(r["meta"]["actual"]["accumulations"]) / dt
    if dt > 1.3 * ro_s:
        print(f"\nTrigger rate seen by the camera: about {rate / 1e3:.1f} kHz "
              f"({dt:.2f} s per {int(r['meta']['actual']['accumulations']):,} triggers).")
    else:
        print(f"\nReadout-limited ({dt:.3f} s per readout): trigger rate is AT LEAST {rate / 1e3:.0f} kHz. "
              "Increase `accumulations` (if the camera allows) to get a real number.")
    print("Compare with the idler singles on the counter. A much higher camera rate means the camera "
          "is triggering on something other than idler photons (noise, ringing, double pulses).")
    return rate


def plot_width_series(runs, threshold=None, spot=None, dark_runs=None, save_as=None):
    """Spot signal and photon hits per readout vs total gate-open time per readout
    (accumulations x width), log-log. Light that arrives at random times (unheralded light,
    dark counts) follows the dashed slope-1 line through the widest run. A flat floor at
    small widths means signal that does not depend on the gate (e.g. light leaking past the
    closed intensifier). Heralded photons would add a constant once the gate covers them."""
    runs = sorted(runs, key=lambda r: r["meta"]["actual"]["gate_width_ns"])
    t_open = np.array([int(r["meta"]["actual"]["accumulations"]) * r["meta"]["actual"]["gate_width_ns"] * 1e-9
                       for r in runs])
    ap = spot_masks(runs[0], spot)[0]
    adu = np.array([_mean_err(spot_series(r, spot)) for r in runs])
    ncol = 2 if threshold is not None else 1
    fig, ax = plt.subplots(1, ncol, figsize=(7 * ncol, 4.6), squeeze=False)
    panels = [(ax[0, 0], adu, "spot ADU per readout (ring-subtracted)")]
    if threshold is not None:
        hits = np.array([_mean_err(hit_series(r, threshold, ap)) for r in runs])
        panels.append((ax[0, 1], hits, f"spot hits >= {threshold:g} ADU above pedestal, per readout"))
    for a_, val, ylab in panels:
        y, e = val[:, 0], val[:, 1]
        pos = y > 0
        a_.errorbar(t_open[pos], y[pos], e[pos], fmt="o", capsize=3, color="C0" if a_ is ax[0, 0] else "C3")
        for xx, yy, r in zip(t_open[pos], y[pos], np.array(runs, dtype=object)[pos]):
            a_.annotate(fmt_ns(r["meta"]["actual"]["gate_width_ns"]), (xx, yy), textcoords="offset points",
                        xytext=(4, 4), fontsize=8)
        if (~pos).any():
            a_.plot(t_open[~pos], np.full((~pos).sum(), max(y[pos].min() / 10 if pos.any() else 1e-3, 1e-3)),
                    "v", color="grey", label="zero or below (plotted at floor)")
        if pos.any():
            ref = np.argmax(t_open * pos)
            xs = np.logspace(np.log10(t_open.min()), np.log10(t_open.max()), 50)
            a_.plot(xs, y[ref] * xs / t_open[ref], "k--", lw=.8, label="proportional to gate time")
        if dark_runs:
            dk = sorted(dark_runs, key=lambda r: r["meta"]["actual"]["gate_width_ns"])
            td = np.array([int(r["meta"]["actual"]["accumulations"]) * r["meta"]["actual"]["gate_width_ns"] * 1e-9
                           for r in dk])
            if a_ is ax[0, 0]:
                vd = np.array([_mean_err(spot_series(r, spot)) for r in dk])
            else:
                vd = np.array([_mean_err(hit_series(r, threshold, ap)) for r in dk])
            ok = vd[:, 0] > 0
            a_.errorbar(td[ok], vd[ok, 0], vd[ok, 1], fmt="s", capsize=3, color="grey", label="dark (signal blocked)")
        a_.set_xscale("log"); a_.set_yscale("log")
        a_.set_xlabel("gate-open time per readout (s) = accumulations x width")
        a_.set_ylabel(ylab)
        a_.legend(fontsize=8)
    a = runs[0]["meta"]["actual"]
    fig.suptitle(f"Gate-width series at delay {a['gate_delay_ns']:g} ns: gain {a['gain']}, "
                 f"{int(a['accumulations']):,} acc x {n_rd(runs[0])} readouts", fontsize=10)
    fig.tight_layout()
    if save_as:
        _savefig(fig, save_as)
    plt.show()


# --------------------------------------------------------------------------------------
# Xu-style Wiener post-processing view
# --------------------------------------------------------------------------------------

def plot_wiener(run, dark, window=9, save_as=None):
    """Xu thesis Sec 4.4 two-pass Wiener filter (postprocessing.wiener_postprocess) on the
    ROI, next to plain subtraction. SNR here = mean(net) / std(background) (Cameron Table 3.1)."""
    from postprocessing import wiener_postprocess
    sig = crop(run, run["sum"])
    bkg = crop(run, scaled_dark_sum(run, dark))
    net = sig - bkg
    s_f, n_f = wiener_postprocess(sig, bkg, window=window)
    snr_before = net.mean() / bkg.std() if bkg.std() > 0 else np.nan
    snr_after = (s_f.mean() - n_f.mean()) / n_f.std() if n_f.std() > 0 else np.nan
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    for a_, img, t, cm in ((ax[0], sig / n_rd(run), "raw ROI (ADU per readout)", "inferno"),
                           (ax[1], net, f"minus dark, SNR {snr_before:.3g}", "inferno"),
                           (ax[2], s_f - n_f.mean(), f"two-pass Wiener {window}x{window}, SNR {snr_after:.3g}",
                            "inferno")):
        im = a_.imshow(img, cmap=cm, vmin=np.percentile(img, 1), vmax=np.percentile(img, 99.7))
        a_.set_title(t)
        a_.set_xlabel("x in ROI (px)"); a_.set_ylabel("y in ROI (px)")
        fig.colorbar(im, ax=a_)
    fig.suptitle(settings_title(run) + "\nSNR = mean(signal - background) / std(background)", fontsize=10)
    fig.tight_layout()
    if save_as:
        _savefig(fig, save_as)
    plt.show()
    return snr_before, snr_after


def _savefig(fig, name):
    _ensure_experiment()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    path = FIG_DIR / f"{name}.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"  figure saved: {path}")
