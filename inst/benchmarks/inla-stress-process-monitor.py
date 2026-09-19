"""Read-only process samples for the full-size joint workload."""
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import time

import psutil

root = Path.cwd()
out = root / "artifacts/inla-stress-calibration"
label = "joint40-n97830-t4"
with (out / "full-process-samples.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=["utc", "pid", "parent_pid", "name",
        "cpu_seconds", "threads", "rss_bytes", "private_bytes", "system_available_bytes"])
    writer.writeheader()
    while True:
        found = {}
        for proc in psutil.process_iter(["pid", "ppid", "name", "cmdline"]):
            args = proc.info["cmdline"] or []
            if label in args and any(str(x).endswith("inla-stress-fit.R") for x in args):
                found[proc.pid] = proc
        roots = [p for p in found.values() if p.info["ppid"] not in found]
        processes = {}
        for proc in roots:
            try:
                for child in [proc, *proc.children(recursive=True)]:
                    processes[child.pid] = child
            except psutil.NoSuchProcess:
                pass
        stamp = datetime.now(timezone.utc).isoformat()
        available = psutil.virtual_memory().available
        for proc in processes.values():
            try:
                mem = proc.memory_info()
                cpu = proc.cpu_times()
                writer.writerow(dict(utc=stamp, pid=proc.pid, parent_pid=proc.ppid(),
                    name=proc.name(), cpu_seconds=cpu.user + cpu.system,
                    threads=proc.num_threads(), rss_bytes=mem.rss,
                    private_bytes=getattr(mem, "private", mem.rss),
                    system_available_bytes=available))
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
        stream.flush()
        # Completion journal is updated after each whole workload group.
        groups = json.loads((out / "stress-groups.json").read_text())
        if any(g["group"] == label for g in groups):
            break
        time.sleep(1)
