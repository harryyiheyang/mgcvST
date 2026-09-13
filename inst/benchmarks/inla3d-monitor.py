"""Sample an existing R process and its INLA children; PID is explicit."""
import csv
import sys
import time
from pathlib import Path

import psutil

root = psutil.Process(int(sys.argv[1]))
out = Path("artifacts/inla3d/eb")
with (out / "memory.csv").open("w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["unix_time", "task", "processes", "rss_MiB", "private_MiB"])
    while root.is_running():
        try:
            procs = [root] + root.children(recursive=True)
        except psutil.NoSuchProcess:
            break
        rss = private = count = 0
        for p in procs:
            try:
                mem = p.memory_info()
            except psutil.NoSuchProcess:
                continue
            rss += mem.rss
            private += getattr(mem, "private", 0)
            count += 1
        task = (out / "current-task.txt").read_text().strip()
        writer.writerow([time.time(), task, count, rss / 2**20, private / 2**20])
        f.flush()
        time.sleep(1)
