"""Run owned R jobs, retain failures, and measure process-tree memory."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import psutil

root = Path.cwd()
out = root / "artifacts/inla-stress-calibration"
rscript = Path("C:/Program Files/R/R-4.6.1/bin/Rscript.exe")
mode = sys.argv[1]
jobs = []
if mode == "null":
    for kind in ("pair", "marginal"):
        for mean in ("0.3", "3"):
            for first in range(1, 501, 10):
                folder = out / "null2d" / f"{kind}-{mean}"
                todo = [r for r in range(first, first + 10) if not (folder / f"rep-{r:04d}.rds").exists()]
                if todo:
                    jobs.append(dict(label=f"2d-{kind}-{mean}-{min(todo):04d}",
                        script="inla-score-null.R", args=[kind, mean, str(min(todo)), str(max(todo))],
                        folder=str(folder), first=min(todo), last=max(todo)))
    for kind in ("independent", "joint"):
        for mean in ("0.3", "3"):
            for first in range(1, 501, 10):
                folder = out / "null3d-paired" / f"{kind}-{mean}"
                todo = [r for r in range(first, first + 10) if not (folder / f"rep-{r:04d}.rds").exists()]
                if todo:
                    jobs.append(dict(label=f"3d-{kind}-{mean}-{min(todo):04d}",
                        script="inla-score-null3d.R", args=[kind, mean, str(min(todo)), str(max(todo))],
                        folder=str(folder), first=min(todo), last=max(todo)))
    groups = [("null", 4, jobs)]
elif mode == "stress":
    groups = []
    for workers, threads in ((1, 1), (2, 1), (4, 1), (8, 1), (4, 2)):
        label = f"independent-w{workers}-t{threads}"
        tasks = [dict(label=f"{label}-{r}", script="inla-stress-fit.R",
            args=[label, "independent", "97830", "1", str(threads), str(r)]) for r in range(1, 11)]
        groups.append((label, workers, tasks))
    for n, threads in ((5000, 4), (97830, 4)):
        label = f"joint40-n{n}-t{threads}"
        tasks = [dict(label=f"{label}-{r}", script="inla-stress-fit.R",
            args=[label, "joint", str(n), "40", str(threads), str(r)]) for r in range(1, 11)]
        groups.append((label, 1, tasks))
else:
    raise ValueError(mode)

logdir = out / f"{mode}-logs"
logdir.mkdir(exist_ok=True)
env = os.environ.copy()
env.update(LC_ALL="C", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
records = []
group_records = []
for group, workers, pending in groups:
    active = []
    started = time.monotonic()
    group_peak_rss = group_peak_private = 0
    min_available = psutil.virtual_memory().available
    while pending or active:
        while pending and len(active) < workers:
            job = pending.pop(0)
            log = (logdir / f"{job['label']}.log").open("w")
            cmd = [str(rscript), str(root / "inst/benchmarks" / job["script"]), *job["args"]]
            proc = subprocess.Popen(cmd, cwd=root, env=env, stdout=log, stderr=subprocess.STDOUT,
                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            active.append(dict(job=job, proc=proc, log=log, start=time.monotonic(),
                peak_rss=0, peak_private=0, status="running"))
        rss_total = private_total = 0
        available = psutil.virtual_memory().available
        min_available = min(min_available, available)
        for item in active[:]:
            proc = item["proc"]
            family = []
            try:
                parent = psutil.Process(proc.pid)
                family = parent.children(recursive=True) + [parent]
            except psutil.NoSuchProcess:
                pass
            rss = private = 0
            for child in family:
                try:
                    mem = child.memory_info()
                    rss += mem.rss
                    private += getattr(mem, "private", mem.rss)
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
            item["peak_rss"] = max(item["peak_rss"], rss)
            item["peak_private"] = max(item["peak_private"], private)
            rss_total += rss
            private_total += private
            if private > 24 * 2**30 or available < 6 * 2**30:
                item["status"] = "resource_stop"
                for child in family:
                    try:
                        child.terminate()
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
            code = proc.poll()
            if code is not None:
                item["log"].close()
                job = item["job"]
                record = dict(label=job["label"], group=group, exit_code=code,
                    status=item["status"] if item["status"] != "running" else ("completed" if code == 0 else "failed"),
                    wall_seconds=time.monotonic()-item["start"], peak_rss_bytes=item["peak_rss"],
                    peak_private_bytes=item["peak_private"])
                records.append(record)
                active.remove(item)
                if code != 0 and mode == "null":
                    folder = Path(job["folder"])
                    missing = [r for r in range(job["first"], job["last"] + 1)
                        if not (folder / f"rep-{r:04d}.rds").exists()]
                    if missing:
                        failed = missing.pop(0)
                        folder.mkdir(parents=True, exist_ok=True)
                        (folder / f"rep-{failed:04d}-failure.json").write_text(json.dumps(record, indent=2))
                        if missing:
                            follow = job.copy()
                            follow.update(label=job["label"]+f"-continue-{min(missing)}",
                                first=min(missing), args=job["args"][:2]+[str(min(missing)), str(job["last"])])
                            pending.append(follow)
                if mode == "stress" and record["status"] == "resource_stop":
                    for remaining in pending:
                        records.append(dict(label=remaining["label"], group=group, exit_code=None,
                            status="not_attempted_after_resource_limit", wall_seconds=0,
                            peak_rss_bytes=0, peak_private_bytes=0))
                    pending.clear()
                (out / f"{mode}-jobs.json").write_text(json.dumps(records, indent=2))
        group_peak_rss = max(group_peak_rss, rss_total)
        group_peak_private = max(group_peak_private, private_total)
        (out / f"{mode}-progress.json").write_text(json.dumps(dict(group=group,
            completed_jobs=len(records), active=len(active), pending=len(pending),
            elapsed_seconds=time.monotonic()-started, current_private_bytes=private_total,
            available_bytes=available), indent=2))
        time.sleep(0.5)
    group_records.append(dict(group=group, workers=workers, wall_seconds=time.monotonic()-started,
        peak_sum_rss_bytes=group_peak_rss, peak_sum_private_bytes=group_peak_private,
        min_system_available_bytes=min_available))
    (out / f"{mode}-groups.json").write_text(json.dumps(group_records, indent=2))
print(json.dumps(dict(mode=mode, jobs=len(records), groups=len(group_records))))
