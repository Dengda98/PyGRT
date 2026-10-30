"""Run compile.sh continuously while Read the Docs waits in short steps."""

import argparse
import os
import signal
import subprocess
import sys
import time
from pathlib import Path


STATE_DIR = Path(".rtd-build")
LOG_FILE = STATE_DIR / "compile.log"
PID_FILE = STATE_DIR / "worker.pid"
RESULT_FILE = STATE_DIR / "exit-code"
OFFSET_FILE = STATE_DIR / "log-offset"
POLL_SECONDS = 45
LOG_BYTES = 6000


def start():
    STATE_DIR.mkdir(exist_ok=True)
    RESULT_FILE.unlink(missing_ok=True)
    OFFSET_FILE.unlink(missing_ok=True)
    with LOG_FILE.open("wb") as output:
        process = subprocess.Popen([sys.executable, "rtd_build.py", "worker"], stdin=subprocess.DEVNULL,
                                   stdout=output, stderr=subprocess.STDOUT, start_new_session=True)
    PID_FILE.write_text(str(process.pid))
    print(f"Started documentation preparation, PID={process.pid}", flush=True)
    return 0


def worker():
    started = time.monotonic()
    print("Running bash compile.sh", flush=True)
    try:
        result = subprocess.run(["bash", "compile.sh"]).returncode
    except OSError as error:
        print(f"Could not start compile.sh: {error}", flush=True)
        result = 1
    if result < 0:
        result = 128 - result
    print(f"Documentation preparation finished: elapsed={time.monotonic() - started:.3f}s, exit={result}", flush=True)
    temporary = RESULT_FILE.with_suffix(".tmp")
    temporary.write_text(str(result))
    temporary.replace(RESULT_FILE)
    return result


def show_log():
    if not LOG_FILE.exists():
        return
    offset = int(OFFSET_FILE.read_text()) if OFFSET_FILE.exists() else 0
    with LOG_FILE.open("rb") as output:
        size = output.seek(0, os.SEEK_END)
        beginning = max(offset, size - LOG_BYTES)
        if beginning > offset:
            print(f"Skipped {beginning - offset} log bytes; showing the latest output", flush=True)
        output.seek(beginning)
        content = output.read(size - beginning)
    if content:
        print(content.decode("utf-8", errors="replace"), end="" if content.endswith(b"\n") else "\n", flush=True)
    OFFSET_FILE.write_text(str(size))


def check(finish=False):
    if not finish:
        deadline = time.monotonic() + POLL_SECONDS
        while not RESULT_FILE.exists() and time.monotonic() < deadline:
            time.sleep(1)
    show_log()
    if RESULT_FILE.exists():
        result = int(RESULT_FILE.read_text())
        print(f"Documentation preparation exit code: {result}", flush=True)
        return result
    if finish:
        # 结束整个后台进程组，避免遗留仍在计算的子进程
        pid = int(PID_FILE.read_text())
        try:
            os.killpg(pid, signal.SIGTERM)
            time.sleep(1)
            os.killpg(pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        print("Documentation preparation exceeded the polling window", flush=True)
        return 124
    print("Documentation preparation is still running", flush=True)
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("start", "worker", "poll", "finish"))
    action = parser.parse_args().action
    if action == "start":
        return start()
    if action == "worker":
        return worker()
    return check(finish=action == "finish")


if __name__ == "__main__":
    sys.exit(main())
