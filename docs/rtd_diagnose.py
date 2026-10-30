"""Run one documentation build step with bounded logs and a timeout."""

import os
import signal
import subprocess
import sys
import time
from pathlib import Path


def main():
    timeout = int(os.environ.get("RTD_STEP_TIMEOUT", "180"))
    command = sys.argv[1:]
    started = time.monotonic()
    print(f"START cwd={Path.cwd()} command={command!r} timeout={timeout}s", flush=True)
    process = subprocess.Popen(command, start_new_session=True)
    try:
        result = process.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        print(f"TIMEOUT after {time.monotonic() - started:.3f}s", flush=True)
        subprocess.run(["ps", "-eo", "pid,ppid,nlwp,pcpu,pmem,etime,args"], check=False)
        for name in ("/sys/fs/cgroup/cpu.stat", "/sys/fs/cgroup/cpu.max", "/proc/loadavg"):
            path = Path(name)
            if path.exists():
                print(f"{name}:\n{path.read_text()}", flush=True)
        os.killpg(process.pid, signal.SIGTERM)
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
        result = 124
    print(f"END elapsed={time.monotonic() - started:.3f}s exit={result}", flush=True)
    return result


if __name__ == "__main__":
    sys.exit(main())
