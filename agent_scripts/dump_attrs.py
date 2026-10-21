"""Dump every camera attribute (name, value, allowed range/labels) to the experiment folder (agent, session 6).
Read-only: changes nothing on the camera except the usual shutdown safe state."""
import sys, pathlib, json
import matplotlib; matplotlib.use("Agg")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import hoi_troubleshoot as ht
exp = ht.new_experiment("agent_attrs")
cam, lim = ht.open_camera()
out = {}
try:
    attrs = cam.get_all_attributes()
    for n, a in attrs.items():
        try:
            v = cam.get_attribute_value(n)
        except Exception as e:
            v = f"<err {e}>"
        d = dict(value=ht._jsonable(v))
        for k in ("min", "max", "step", "labels", "values", "readable", "writable", "kind"):
            if hasattr(a, k):
                try: d[k] = ht._jsonable(getattr(a, k))
                except Exception: pass
        out[n] = d
    out["_gate_limits"] = ht._jsonable(lim)
finally:
    ht.shutdown(cam)
(exp / "figures" / "attrs.json").write_text(json.dumps(out, indent=1))
for n, d in out.items():
    print(f"{n:40s} {d.get('value') if isinstance(d, dict) else d}   {d.get('labels','') if isinstance(d, dict) else ''}")
