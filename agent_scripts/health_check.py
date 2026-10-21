"""Camera health check with the INTERNAL trigger (agent session 3, Tristan: "just see if the camera works").

    python agent_scripts/health_check.py --tag health --passes 2

No herald involved: the gate is opened by the camera's own clock, so every hit is unheralded light
(stray/room light, any signal-arm light) or intensifier/CCD background. Phases:
  S  switch:   External (herald) <-> Internal, gain 1: read-back and readout timing must follow
  A  timing:   gain 1, 3 ns gate, several (rate, accumulations): read-back of Trigger Source /
               Frequency and readout interval vs accumulations / rate (does Internal really take?)
  B  gain:     10 kHz x 10k acc (1 s wall), 10 us gate (0.1 s open) and 3 ns gate, gains 100..9000
  C  gate:     gain 9000, same wall time, gate widths 3 ns .. 50 us; plus equal gate-open time (0.1 s)
               at 1, 10, 100 kHz and 4 s wall time (splits in-gate light from wall-clock background)
  D  closed:   3 ns gate at gain 9000, wall time 0.25-4 s; intensifier off (EM register only); gain 1
Every run: read-back recorded, saturation (> 0.1 % px) aborts, temperature status checked.
Analysis: agent_scripts/health_analyze.py on the experiment folder.
"""
import argparse, json, sys, pathlib, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht

p = argparse.ArgumentParser()
p.add_argument("--tag", default="health"); p.add_argument("--passes", type=int, default=2)
p.add_argument("--phases", default="A,B,C")
a = p.parse_args()
exp = ht.new_experiment(f"agent_{a.tag}"); fig = exp / "figures"
base = dict(trigger_source="Internal", gate_delay_ns=27.0, roi_center=(500, 455), kind="internal",
            photon_mode="Disabled", notes=f"agent {a.tag}", frame_timeout_s=30)
plan = []
if "S" in a.phases:   # switching test: External (APD herald) then Internal, twice; timing must change
    for k in range(2):
        plan.append(("S", dict(gain=1, gate_width_ns=3.0, trigger_source="External", trigger_threshold_v=0.5,
                               accumulations=20000, n_readouts=5, label=f"S_ext_20000acc_{k + 1}")))
        plan.append(("S", dict(gain=1, gate_width_ns=3.0, internal_trigger_hz=1000, accumulations=1000,
                               n_readouts=5, label=f"S_int_1k_1000acc_{k + 1}")))
if "A" in a.phases:
    for hz, acc in [(1000, 500), (1000, 1000), (1000, 2000), (10000, 2000), (10000, 10000), (100000, 100000)]:
        plan.append(("A", dict(gain=1, gate_width_ns=3.0, internal_trigger_hz=hz, accumulations=acc, n_readouts=8,
                               label=f"A_{hz // 1000}k_{acc}acc")))
if "B" in a.phases:
    for g in (100, 1000, 3000, 6000, 9000):
        for w in (3.0, 10000.0):
            plan.append(("B", dict(gain=g, gate_width_ns=w, internal_trigger_hz=10000, accumulations=10000,
                                   n_readouts=15, label=f"B_g{g}_w{int(w)}ns")))
C = [(10000, 10000, w) for w in (3.0, 100.0, 1000.0, 10000.0, 50000.0)] + \
    [(1000, 1000, 100000.0), (1000, 1000, 3.0), (100000, 100000, 1000.0), (100000, 100000, 3.0),
     (1000, 4000, 25000.0), (1000, 4000, 3.0)]
if "C" in a.phases:
    for k in range(a.passes):
        order = C if k % 2 == 0 else C[::-1]
        for hz, acc, w in order:
            plan.append(("C", dict(gain=9000, gate_width_ns=w, internal_trigger_hz=hz, accumulations=acc,
                                   n_readouts=10 if acc == 4000 else 20,
                                   label=f"C_{hz // 1000}k_{acc}acc_w{int(w)}ns_p{k + 1}")))
if "D" in a.phases:   # closed gate (3 ns) vs wall time; intensifier off; gain 1 (what is the closed-gate background?)
    for k in range(a.passes):
        accs = (250, 500, 1000, 2000, 3000, 4000)
        for acc in (accs if k % 2 == 0 else accs[::-1]):
            plan.append(("D", dict(gain=9000, gate_width_ns=3.0, internal_trigger_hz=1000, accumulations=acc,
                                   n_readouts=8, label=f"D_1k_{acc}acc_p{k + 1}")))
    for acc in (1000, 4000):
        plan.append(("D", dict(gain=9000, gate_width_ns=3.0, internal_trigger_hz=1000, accumulations=acc,
                               n_readouts=8, label=f"D_1k_{acc}acc_intensifier_off")))
    plan.append(("D", dict(gain=1, gate_width_ns=3.0, internal_trigger_hz=1000, accumulations=4000, n_readouts=8,
                           label="D_1k_4000acc_gain1")))
    plan.append(("D", dict(gain=9000, gate_width_ns=3.0, internal_trigger_hz=10000, accumulations=40000,
                           n_readouts=8, label="D_10k_40000acc")))
print(f"{len(plan)} runs planned")
_configure = ht.configure
INTENSIFIER_OFF = [False]


def configure_maybe_off(cam, s, **kw):
    """ht.configure, then (for *_intensifier_off runs) switch the intensifier off and check the read-back."""
    _configure(cam, s, **kw)
    if INTENSIFIER_OFF[0]:
        import camera_control as cc
        cc.safe_set_attribute_value(cam, "Enable Intensifier", False); cc.safe_commit_parameters(cam)
        got = cam.get_attribute_value("Enable Intensifier")
        if bool(got):
            raise RuntimeError(f"Enable Intensifier read back {got!r}")
        print("  Enable Intensifier -> False (read back)")


ht.configure = configure_maybe_off
rec = []
cam, limits = ht.open_camera()
try:
    for ph, s in plan:
        INTENSIFIER_OFF[0] = s["label"].endswith("intensifier_off")
        r = ht.take_run(cam, **{**base, **s})
        if r is None:
            raise RuntimeError("no readouts (internal trigger not running?)")
        act, rb = r["meta"]["actual"], r["meta"]["readback"]
        sat = float((r["first"] >= ht.SATURATION_ADU).mean())
        st = r["roi_stack"].astype(float)
        ped = np.median(st.reshape(len(st), -1), axis=1)
        d = dict(phase=ph, label=s["label"], intensifier=rb.get("Enable Intensifier"), file=r["meta"]["file"], want=s, trig=rb.get("Trigger Source"),
                 freq=rb.get("Trigger Frequency"), gain_rb=rb.get("emICCD Gain"), acc_rb=rb.get("Accumulations"),
                 gate_rb=[rb.get("Repetitive Gate delay (ns)"), rb.get("Repetitive Gate width (ns)")],
                 temp=act["temp_c"], temp_status=act["temp_status"], dt=ht.readout_interval_s(r),
                 sat=sat, ped=float(ped.mean()), hits200=float((st - ped[:, None, None] >= 200).sum((1, 2)).mean()))
        rec.append(d)
        print("REC", json.dumps({k: v for k, v in d.items() if k not in ("file", "want")}), flush=True)
        (fig / "records.json").write_text(json.dumps(rec, indent=1))
        if str(d["trig"]) != s.get("trigger_source", "Internal"):
            raise RuntimeError(f"Trigger Source read back {d['trig']!r}")
        if sat > 1e-3:
            raise RuntimeError(f"saturation {sat:.4%} > 0.1 %: stopping")
        if str(d["temp_status"]).lower() != "locked":
            print("  WARNING temperature not locked")
finally:
    ht.shutdown(cam)
print("EXP", exp)
