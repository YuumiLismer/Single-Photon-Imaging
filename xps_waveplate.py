"""
xps_waveplate.py -- minimal client for the Newport XPS-D (192.168.0.254, TCP port 5001)
that drives the pump-power half-wave plate (PR50CC rotation stage, positioner S.Pos2 = slot 2).

    import xps_waveplate as xw
    wp = xw.Waveplate()          # connects; does not move anything
    wp.position()                # degrees
    wp.move_to(30.0)             # enables the group if needed, absolute move, waits
    wp.close(disable=True)       # optionally turn the motors off again (as found on 22 Oct)
"""

import socket

HOST, PORT = "192.168.0.254", 5001


class XPSError(RuntimeError):
    pass


class XPS:
    def __init__(self, host=HOST, port=PORT, timeout=30.0):
        self.sock = socket.create_connection((host, port), timeout=timeout)

    def query(self, cmd):
        """Send one XPS API call, return the comma-separated fields after the error code."""
        self.sock.sendall(cmd.encode())
        buf = b""
        while b"EndOfAPI" not in buf:
            chunk = self.sock.recv(4096)
            if not chunk:
                raise XPSError(f"connection closed during {cmd}")
            buf += chunk
        code, _, rest = buf.decode().replace(",EndOfAPI", "").partition(",")
        if int(code) != 0:
            raise XPSError(f"{cmd} -> error {code} {rest}")
        return rest.split(",") if rest else []

    def close(self):
        self.sock.close()


class Waveplate:
    def __init__(self, positioner="S.Pos2", xps=None):
        self.positioner = positioner
        self.group = positioner.split(".")[0]
        self.xps = xps or XPS()
        self.initial_status = self.status()

    def status(self):
        return int(self.xps.query(f"GroupStatusGet({self.group},int *)")[0])

    def position(self):
        return float(self.xps.query(f"GroupPositionCurrentGet({self.positioner},double *)")[0])

    def enable(self):
        s = self.status()
        if 0 <= s <= 10 or s == 42:
            raise XPSError(f"group {self.group} is not initialized/referenced (status {s}); "
                           "home it from the XPS web page first")
        if 20 <= s <= 39:                     # disabled states
            self.xps.query(f"GroupMotionEnable({self.group})")

    def move_to(self, deg):
        """Absolute move in degrees. Blocks until the move is done."""
        self.enable()
        self.xps.query(f"GroupMoveAbsolute({self.positioner},{float(deg):.4f})")
        return self.position()

    def close(self, disable=False):
        if disable:
            self.xps.query(f"GroupMotionDisable({self.group})")
        self.xps.close()
