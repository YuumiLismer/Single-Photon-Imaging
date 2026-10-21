"""External trigger-rate check (agent, session 4): short readouts (<= ~1 s), 0.5 V rising, 50 ohm, DC."""
import sys, json, pathlib, datetime
import matplotlib
matplotlib.use("Agg")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht

exp = ht.new_experiment("agent_rate_check")
base = dict(gain=9000, trigger_source="External", trigger_threshold_v=0.5, roi_center=(500, 455),
            notes="agent session 4 rate check")
cam, limits = ht.open_camera()
rates = []
try:
    ht.cool_down(cam, setpoint_c=-20.0, poll_s=30)
    for acc in (20_000, 30_000):
        rates.append(dict(acc=acc, rate_hz=ht.measure_trigger_rate(cam, base, accumulations=acc, n_readouts=6)))
finally:
    ht.shutdown(cam)
(exp / "figures" / "trigger_rate.json").write_text(json.dumps(ht._jsonable(dict(rates=rates, when=str(datetime.datetime.now())))))
print("RATES", rates)
