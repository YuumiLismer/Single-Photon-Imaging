"""Interleaved gate-window scan for the heralded delay (agent, 20 Oct).

    python agent_scripts/window_scan.py --start 27 --stop 2027 --width 100 --acc 50000 \
        --readouts 15 --passes 3 --tag exp3

Why not photon_hunt's coarse step: the beam at the camera is now so weak (a few detected
photons/s, see agent_log E2) that 1 us windows are background-dominated. Narrow windows with
many accumulations per readout keep the trigger time but cut the in-gate background. Each pass
visits every window once in a shuffled order, so slow drifts (pump, trigger rate) hit all
windows alike; every window is scored against all the others pooled, after every pass.
Score = photon-hit pixels (>= --thr ADU above the readout median) in the spot aperture minus the
area-scaled hits at r > 50 px in the same readout (drift-free; see agent_log E3).
"""
import argparse, json, sys, pathlib, datetime
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht
from photon_hunt import merge_runs

p = argparse.ArgumentParser()
p.add_argument("--start", type=float, default=27); p.add_argument("--stop", type=float, default=2027)
p.add_argument("--width", type=float, default=100); p.add_argument("--acc", type=int, default=50_000)
p.add_argument("--readouts", type=int, default=15); p.add_argument("--passes", type=int, default=3)
p.add_argument("--spot", type=int, nargs=2, default=(501, 454)); p.add_argument("--thr", type=float, default=200)
p.add_argument("--gain", type=int, default=9000); p.add_argument("--tag", default="wscan")
p.add_argument("--stop-z", type=float, default=6.0, help="stop early once a window reaches this z")
a = p.parse_args()
spot = tuple(a.spot)
stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
folder = ht.new_experiment(f"agent_{a.tag}") / "figures"   # runs go to the same experiment folder
delays = np.arange(a.start, a.stop, a.width)
base = dict(gain=a.gain, gate_width_ns=a.width, accumulations=a.acc, n_readouts=a.readouts,
            trigger_source="External", trigger_threshold_v=0.5, roi_center=(500, 455), notes=f"agent {a.tag}")
runs = {}
rng = np.random.default_rng(1)


def score(runs):
    H, A, H20, N = {}, {}, {}, {}
    for d, r in runs.items():
        ap = ht.spot_masks(r, spot)[0]
        x0, y0, n = r["meta"]["roi"]
        yy, xx = np.indices((n, n)); big = np.hypot(xx - (spot[0] - x0), yy - (spot[1] - y0)) <= 20
        far = np.hypot(xx - (spot[0] - x0), yy - (spot[1] - y0)) > 50
        hit = ht.above_pedestal(r) >= a.thr
        # net = spot hits minus area-scaled far-field hits of the same readout: removes the ungated,
        # readout-time-proportional background, which follows the trigger-rate drift (exp3: corr 0.86)
        N[d] = hit[:, ap].sum(1) - ap.sum() * hit[:, far].mean(1)
        H[d] = ht.hit_series(r, a.thr, ap); A[d] = ht.spot_series(r, spot); H20[d] = ht.hit_series(r, a.thr, big)
    out = []
    for d in sorted(runs):
        others = lambda S: np.concatenate([S[e] for e in S if e != d])
        out.append(dict(delay=d, n=len(H[d]), net=N[d].mean(), z_net=ht._diff_z(N[d], others(N)),
                        hits=H[d].mean(), z_hits=ht._diff_z(H[d], others(H)),
                        hits20=H20[d].mean(), z_hits20=ht._diff_z(H20[d], others(H20)),
                        adu=A[d].mean(), z_adu=ht._diff_z(A[d], others(A))))
    return out


cam, limits = ht.open_camera()
try:
    ht.cool_down(cam, setpoint_c=-20.0, poll_s=30)
    for k in range(a.passes):
        order = rng.permutation(delays)
        for d in order:
            r = ht.take_run(cam, **{**base, "label": f"{a.tag}_p{k + 1}_{d:g}ns", "gate_delay_ns": float(d)})
            if r is None:
                raise RuntimeError("no readouts (triggers stopped?)")
            if float((r["first"] >= ht.SATURATION_ADU).mean()) > 1e-3:
                raise RuntimeError("saturation > 0.1%: stopping")
            runs[d] = r if d not in runs else merge_runs(runs[d], r)
        sc = score(runs)
        best = max(sc, key=lambda s: s["z_net"])
        print(f"PASS {k + 1}: best {best['delay']:g} ns z_net {best['z_net']:+.1f} z_hits {best['z_hits']:+.1f} "
              f"(hits {best['hits']:.3f}/readout, r20 z {best['z_hits20']:+.1f}, ADU z {best['z_adu']:+.1f})", flush=True)
        for s in sc:
            print(f"  {s['delay']:7g} ns n={s['n']:3d} net {s['net']:+.2f} z {s['z_net']:+5.1f} | hits {s['hits']:.3f} z {s['z_hits']:+5.1f} | r20 {s['hits20']:.2f} "
                  f"z {s['z_hits20']:+5.1f} | ADU {s['adu']:8.0f} z {s['z_adu']:+5.1f}", flush=True)
        (folder / "scores.json").write_text(json.dumps(ht._jsonable(dict(args=vars(a), passes=k + 1, scores=sc)), indent=1))
        if best["z_net"] >= a.stop_z:
            print("STOP: window reached stop-z"); break
    sc = score(runs)
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.2))
    D = [s["delay"] for s in sc]
    for i, (key, lab) in enumerate((("net", "net spot hits (r<=10 minus far field)"), ("hits20", "hits r<=20"), ("adu", "spot ADU (ring-sub)"))):
        ax[i].plot(D, [s[key] for s in sc], "o-"); ax[i].set_xlabel("gate delay (ns)"); ax[i].set_ylabel(lab + " per readout")
        ax[i].set_title(" ".join(f"{s['delay']:g}:{s['z_' + key]:+.1f}" for s in sorted(sc, key=lambda s: -s['z_' + key])[:3]), fontsize=8)
    fig.suptitle(f"{a.tag}: {ht.fmt_ns(a.width)} windows, gain {a.gain}, {a.acc:,} acc x {sc[0]['n']} readouts each, spot {spot}")
    fig.tight_layout(); fig.savefig(folder / "window_scan.png", dpi=110)
    ht.plot_delay_scan([runs[d] for d in sorted(runs)], threshold=a.thr, spot=spot, save_as="delay_scan")
finally:
    ht.shutdown(cam)
print("FIGS", folder)
