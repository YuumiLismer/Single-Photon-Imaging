"""Live zoomed view of a beam on the PI-MAX4 with a live profile across the fringes underneath (diode alignment,
interference). Internal trigger, gain 1 by default; same light safety as live_view.py.

    python agent_scripts/live_fringes.py --center 504 450 --half 40 --gate 6400 --acc 2048

Keys (click the window first):
    +  / -     gate width x2 / /2  (3 ns .. 10 us)
    ]  / [     accumulations x2 / /2  (capped at ~0.5 s of gates per readout)
    a          re-centre the zoom on the brightest (smoothed) spot in the full frame
    z  / x     zoom in / out
    f          full frame
    r          profile across rows (vertical stripes) <-> across columns (horizontal stripes)
    s          save the current frame (.npy + .png)
    q          quit (camera goes to gain 1, intensifier off, closed)
"""
import argparse
import pathlib
import sys
import tempfile
import time

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider
from scipy.ndimage import gaussian_filter

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht

SAT_ADU = 60_000
TRIGGER_HZ = 10_000
MAX_ACC = TRIGGER_HZ // 2


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gate", type=float, default=6400.0, help="starting gate width (ns)")
    p.add_argument("--acc", type=int, default=2048, help="starting accumulations per readout")
    p.add_argument("--gain", type=int, default=1)
    p.add_argument("--center", type=int, nargs=2, default=(504, 450), help="zoom centre x y (px)")
    p.add_argument("--half", type=int, default=40, help="zoom half-width (px)")
    p.add_argument("--band", type=int, default=9, help="profile band half-height (px) around the centre")
    p.add_argument("--dump", type=float, default=0, help="every this many s, write the latest frame to <tmp>/live_fringes_latest.npy")
    p.add_argument("--recentre", action="store_true", help="centre on the brightest spot on the first frame")
    p.add_argument("--cmd-file", default=str(pathlib.Path(tempfile.gettempdir()) / "live_fringes_cmd.txt"),
                   help="keys written to this file are applied as if pressed (remote control)")
    a = p.parse_args(argv)
    import camera_control as cc

    st = dict(gate=float(np.clip(a.gate, 3, 10_000)), acc=int(np.clip(a.acc, 1, MAX_ACC)), cx=a.center[0], cy=a.center[1],
              half=a.half, rows=True, dirty=True, view=True, quit=False, save=False, recentre=a.recentre)
    cam, limits = ht.open_camera()
    w_min = max(3.0, float(limits["width_min_ns"]))
    base = {**ht.DEFAULT_SETTINGS, "label": "live_fringes", "gain": a.gain, "trigger_source": "Internal",
            "internal_trigger_hz": TRIGGER_HZ, "gate_delay_ns": 100.0, "photon_mode": "Disabled"}

    def apply():
        try:
            cam.stop_acquisition()
        except Exception:
            pass
        st["gate"] = float(np.clip(st["gate"], w_min, 10_000))
        st["acc"] = int(np.clip(st["acc"], 1, MAX_ACC))
        ht.configure(cam, {**base, "gate_width_ns": st["gate"], "accumulations": st["acc"]}, verbose=False)
        cam.start_acquisition()
        st["dirty"] = False

    def on_key(ev):
        k = ev if isinstance(ev, str) else ev.key
        if k.startswith("half="):                              # remote only: half=150, centre=493,616
            st["half"] = int(k[5:]); st["view"] = True; return
        if k.startswith("centre="):
            st["cx"], st["cy"] = (int(v) for v in k[7:].split(",")); st["view"] = True; return
        if k in ("+", "="):
            st["gate"] *= 2; st["dirty"] = True
        elif k in ("-", "_"):
            st["gate"] /= 2; st["dirty"] = True
        elif k == "]":
            st["acc"] *= 2; st["dirty"] = True
        elif k == "[":
            st["acc"] //= 2; st["dirty"] = True
        elif k == "a":
            st["recentre"] = True
        elif k == "z":
            st["half"] = max(8, int(st["half"] / 1.5)); st["view"] = True
        elif k == "x":
            st["half"] = min(512, int(st["half"] * 1.5)); st["view"] = True
        elif k == "f":
            st["cx"], st["cy"], st["half"] = 512, 512, 512; st["view"] = True
        elif k == "r":
            st["rows"] = not st["rows"]; st["view"] = True
        elif k == "s":
            st["save"] = True
        elif k in ("q", "escape"):
            print(f"quit key {k!r}", flush=True); st["quit"] = True

    plt.ion()
    fig = plt.figure(figsize=(9, 11.5))
    ax = fig.add_axes([0.09, 0.40, 0.86, 0.53]); axp = fig.add_axes([0.09, 0.15, 0.86, 0.19])
    # on-screen controls: box-size slider, zoom/full/find buttons; a click on the image centres the view there
    sld = Slider(fig.add_axes([0.22, 0.075, 0.6, 0.025]), "box size (px)", 16, 1024, valinit=2 * st["half"], valstep=2)
    sld.on_changed(lambda v: st.update(half=max(8, int(v) // 2), view=True))
    buttons = []
    for i, (txt, key) in enumerate([("Zoom in", "z"), ("Zoom out", "x"), ("Full frame", "f"), ("Find beam", "a"), ("Save", "s")]):
        b = Button(fig.add_axes([0.09 + i * 0.175, 0.015, 0.16, 0.045]), txt)
        b.on_clicked(lambda ev, key=key: on_key(key)); buttons.append(b)
    fig.canvas.mpl_connect("button_press_event", lambda ev: ev.inaxes is ax and ev.button == 1 and ev.xdata is not None
                           and on_key(f"centre={int(ev.xdata)},{int(ev.ydata)}"))
    fig.canvas.mpl_connect("key_press_event", on_key)
    try:
        fig.canvas.manager.window.state("zoomed")              # maximised window (TkAgg)
    except Exception:
        pass
    fig.canvas.mpl_connect("close_event", lambda ev: (print("window closed", flush=True), st.update(quit=True)))
    im = ax.imshow(np.zeros((2, 2)), cmap="inferno", origin="upper", interpolation="nearest")
    band_lines = [ax.axhline(0, color="cyan", lw=0.8, ls="--") for _ in range(2)]
    line, = axp.plot([], [], "o-", ms=2.5, color="C3")
    axp.set_ylabel("ADU above median"); axp.grid(alpha=0.3)
    vis_txt = axp.text(0.01, 0.95, "", transform=axp.transAxes, va="top", fontsize=9)
    save_dir = None
    n, t_last = 0, time.time()
    try:
        cmd = pathlib.Path(a.cmd_file)
        while not st["quit"]:
            if cmd.exists():                                    # remote key presses
                try:
                    keys = cmd.read_text().split(); cmd.unlink()
                except OSError:
                    keys = []
                for k in keys:
                    on_key(k)
            try:
                if st["dirty"]:
                    apply()
                if not cc.poll_frame(cam, 0.2):
                    plt.pause(0.01)
                    continue
                frames = cam.read_multiple_images()
            except Exception as e:                              # camera hiccup (acquisition dropped): restart it
                print(f"{time.strftime('%H:%M:%S')} camera error, restarting acquisition: {e}", flush=True)
                st["dirty"] = True
                plt.pause(2)
                continue
            if not frames:
                continue
            img = np.asarray(frames[-1], dtype=float)
            n += 1
            n_sat = int(np.count_nonzero(img >= SAT_ADU))
            if n_sat > 25:
                if st["gate"] <= w_min * 1.01 and st["acc"] <= 1:
                    print("SATURATING at the shortest exposure: intensifier off, stopping.")
                    break
                if st["acc"] > 1:
                    st["acc"] //= 2
                else:
                    st["gate"] /= 2
                st["dirty"] = True
                print(f"saturation ({n_sat} px): exposure reduced", flush=True)
            img = img - np.median(img)
            if a.dump and time.time() - st.get("t_dump", 0) > a.dump:     # for remote monitoring
                np.save(pathlib.Path(tempfile.gettempdir()) / "live_fringes_latest.npy", img); st["t_dump"] = time.time()
            if st["recentre"]:
                sm = gaussian_filter(img, 4)
                st["cy"], st["cx"] = (int(v) for v in np.unravel_index(sm.argmax(), sm.shape))
                st["recentre"] = False; st["view"] = True
                print(f"recentred on ({st['cx']}, {st['cy']}), smoothed peak {sm.max():.0f} ADU", flush=True)
            h, cx, cy = st["half"], st["cx"], st["cy"]
            h = min(h, 512)                                     # keep the whole box on the sensor (slide it, don't crop it)
            x0, y0 = (int(np.clip(c - h, 0, 1024 - 2 * h)) for c in (cx, cy)); x1, y1 = x0 + 2 * h, y0 + 2 * h
            cut = img[y0:y1, x0:x1]
            im.set_data(cut); im.set_extent([x0, x1, y1, y0])
            lo, hi = np.percentile(cut, [1, 99.9]); im.set_clim(lo, max(hi, lo + 1))
            b = a.band
            if st["rows"]:                                      # vertical stripes: profile along x, averaged over rows
                xs, prof = np.arange(x0, x1), img[max(0, cy - b):cy + b + 1, x0:x1].mean(0)
                pos = (cy - b, cy + b)
            else:                                               # horizontal stripes: profile along y
                xs, prof = np.arange(y0, y1), img[y0:y1, max(0, cx - b):cx + b + 1].mean(1)
                pos = (cx - b, cx + b)
            if st["view"]:
                if sld.val != min(1024, 2 * h):
                    sld.set_val(min(1024, 2 * h))
                ax.set_xlim(x0, x1); ax.set_ylim(y1, y0)
                for ln, v in zip(band_lines, pos):
                    if st["rows"]:
                        ln.set_transform(ax.get_yaxis_transform()); ln.set_data([0, 1], [v, v])
                    else:
                        ln.set_transform(ax.get_xaxis_transform()); ln.set_data([v, v], [0, 1])
                axp.set_xlabel("x (px)" if st["rows"] else "y (px)")
                st["view"] = False
            line.set_data(xs, prof)
            axp.set_xlim(xs[0], xs[-1]); axp.set_ylim(min(0, prof.min()) - 5, prof.max() * 1.1 + 5)
            # local fringe contrast near the profile peak: (max - min) / (max + min) of the fine structure within +-8 px
            ip = int(np.argmax(gaussian_filter(prof, 3))); seg = prof[max(0, ip - 8):ip + 9]
            vis = (seg.max() - seg.min()) / max(seg.max() + seg.min(), 1)
            vis_txt.set_text(f"peak {prof.max():.0f} ADU at {xs[np.argmax(prof)]}   local contrast (max-min)/(max+min) within +-8 px: {vis:.2f}")
            now = time.time()
            fps, t_last = 1 / max(now - t_last, 1e-3), now
            ax.set_title(f"LIVE  gain {a.gain}  gate {ht.fmt_ns(st['gate'])} x {st['acc']} acc  |  zoom ({cx},{cy}) +-{h} px  |  "
                         f"max {cut.max():.0f} ADU  sat {n_sat} px  |  {fps:.1f} fps\n"
                         "keys: +/- gate  ]/[ acc  a recentre  z/x zoom  r rows/cols  s save  q quit", fontsize=9)
            if st["save"]:
                if save_dir is None:
                    save_dir = ht.new_experiment("live_fringes") / "figures"
                    save_dir.mkdir(parents=True, exist_ok=True)
                stem = save_dir / f"live_{time.strftime('%H%M%S')}"
                np.save(stem.with_suffix(".npy"), img)
                fig.savefig(stem.with_suffix(".png"), dpi=110)
                print("saved", stem.with_suffix(".png"), flush=True)
                st["save"] = False
            plt.pause(0.01)
    except BaseException as e:
        print("stopped by", repr(e), flush=True); raise
    finally:
        print(f"loop ended after {n} frames", flush=True)
        try:
            cam.stop_acquisition()
        except Exception:
            pass
        ht.shutdown(cam)
        plt.close("all")


if __name__ == "__main__":
    main()
