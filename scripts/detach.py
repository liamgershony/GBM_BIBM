#!/usr/bin/env python3
"""setsid(1) equivalent for macOS: double-fork, os.setsid(), exec the command.

macOS ships no setsid. `nohup` alone only ignores SIGHUP, so a background job
launched from a short-lived shell is still killed when that shell's process group
is signalled -- which is how the first optimal-transport run died silently.

    scripts/detach.py <pidfile> <command...>
"""
import os
import sys

pidfile, cmd = sys.argv[1], sys.argv[2:]
if not cmd:
    sys.exit("usage: detach.py <pidfile> <command...>")

# First fork: the child is guaranteed not to be a process-group leader, so
# setsid() can succeed and put it in a brand-new session with no controlling tty.
if os.fork() > 0:
    os._exit(0)
os.setsid()
pid = os.fork()
if pid > 0:
    with open(pidfile, "w") as fh:
        fh.write(str(pid))
    os._exit(0)
os.execvp(cmd[0], cmd)
