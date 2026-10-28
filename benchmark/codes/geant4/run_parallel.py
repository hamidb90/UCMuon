#!/usr/bin/env python3
"""Run MuonRock on a FILE-mode source as K independent serial jobs, W at a time,
and merge them into one run with the files and event numbering of a single run.

    python3 run_parallel.py --exe build/MuonRock --source muons.txt \\
        --out outputs --name run_stdrock_mono --chunks 32 --workers 8 --seed 20261019

The source lines are dealt round-robin (line i to chunk i mod K), so every
chunk sees the whole population and chunk c's event j is source line j*K + c:
the merged EventID is that source index, as in a single run. Chunk c uses
random seed SEED + c (MuonRock -s). Geant4's environment must be set first
(source geant4.sh).

Writes <out>/<name>_{muons,secondaries,stopped}.csv sorted by EventID,
<name>_summary.txt (Elapsed = CPU time summed over the chunks, the
single-core-equivalent time; the wall time is stated separately) and
<name>_chunks.csv (events, seed, wall and CPU seconds of every chunk).
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

KINDS = ("muons", "secondaries", "stopped")


def split_source(source: Path, k: int, chunk_root: Path) -> tuple[list[Path], int]:
    header, data = [], []
    with open(source) as fh:
        for line in fh:
            (header if line.startswith("#") or not line.strip() else data).append(line)
    paths = []
    for c in range(k):
        d = chunk_root / f"chunk_{c:03d}"
        d.mkdir(parents=True, exist_ok=True)
        p = d / "source.txt"
        with open(p, "w") as fh:
            fh.writelines(header)
            fh.writelines(data[c::k])
        paths.append(p)
    return paths, len(data)


def run_chunks(exe: Path, sources: list[Path], workers: int, seed: int,
               extra: list[str]) -> list[dict]:
    pending = list(enumerate(sources))
    running, done = {}, []
    while pending or running:
        while pending and len(running) < workers:
            c, src = pending.pop(0)
            n = sum(1 for line in open(src) if line.strip() and not line.startswith("#"))
            cmd = [str(exe), "-f", str(src), "-o", str(src.parent / "out"), "-s", str(seed + c),
                   *extra, "-q", str(n)]
            proc = subprocess.Popen(cmd, cwd=src.parent, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            running[proc.pid] = dict(chunk=c, events=n, seed=seed + c, t0=time.perf_counter(), proc=proc)
        pid, status, ru = os.wait4(-1, 0)
        job = running.pop(pid, None)
        if job is None:
            continue
        job["wall_s"] = time.perf_counter() - job.pop("t0")
        job["cpu_s"] = ru.ru_utime + ru.ru_stime
        job.pop("proc")
        if os.waitstatus_to_exitcode(status) != 0:
            sys.exit(f"chunk {job['chunk']} failed (exit {os.waitstatus_to_exitcode(status)})")
        done.append(job)
        print(f"  chunk {job['chunk']:3d}: {job['events']} events, {job['wall_s']:.0f} s wall, "
              f"{job['cpu_s']:.0f} s CPU   ({len(done)}/{len(sources)})", flush=True)
    return sorted(done, key=lambda j: j["chunk"])


def merge(chunk_root: Path, k: int, out: Path, name: str) -> dict:
    counts = {}
    for kind in KINDS:
        rows, header = [], None
        for c in range(k):
            files = glob.glob(str(chunk_root / f"chunk_{c:03d}" / "out" / f"run_*_{kind}.csv"))
            if len(files) != 1:
                sys.exit(f"chunk {c}: expected one {kind} file, found {len(files)}")
            with open(files[0]) as fh:
                h = fh.readline()
                header = header or h
                for seq, line in enumerate(fh):
                    eid, rest = line.split(",", 1)
                    rows.append((int(eid) * k + c, seq, rest))
        rows.sort()
        with open(out / f"{name}_{kind}.csv", "w") as fh:
            fh.write(header)
            fh.writelines(f"{g},{rest}" for g, _, rest in rows)
        counts[kind] = len(rows)
    return counts


def chunk_summaries(chunk_root: Path, k: int) -> list[dict]:
    out = []
    for c in range(k):
        text = open(glob.glob(str(chunk_root / f"chunk_{c:03d}" / "out" / "run_*_summary.txt"))[0]).read()
        get = lambda key: re.search(rf"{key}\s*:\s*([\d.]+)", text)
        out.append({key: float(get(key).group(1)) for key in ("KE min", "KE max") if get(key)})
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--exe", type=Path, required=True)
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--chunks", type=int, default=32)
    ap.add_argument("--workers", type=int, default=max(1, int(0.8 * (os.cpu_count() or 1))))
    ap.add_argument("--seed", type=int, default=20261019)
    ap.add_argument("--keep-chunks", action="store_true")
    ap.add_argument("--xy", type=float, help="slab half-width in cm (MuonRock -xy)")
    ap.add_argument("--incl", type=float, default=0.0,
                    help="beam inclination in degrees (MuonRock -incl): planes at L cos(incl)")
    a = ap.parse_args()

    a.out.mkdir(parents=True, exist_ok=True)
    chunk_root = a.out / f"{a.name}_chunks"
    if chunk_root.exists():
        shutil.rmtree(chunk_root)
    print(f"{a.name}: {a.source} in {a.chunks} chunks, {a.workers} at a time, seeds {a.seed}..{a.seed + a.chunks - 1}")
    sources, n_total = split_source(a.source, a.chunks, chunk_root)
    t0 = time.perf_counter()
    extra = (["-xy", str(a.xy)] if a.xy else []) + (["-incl", str(a.incl)] if a.incl else [])
    jobs = run_chunks(a.exe.resolve(), sources, a.workers, a.seed, extra)
    wall = time.perf_counter() - t0
    counts = merge(chunk_root, a.chunks, a.out, a.name)
    ke = chunk_summaries(chunk_root, a.chunks)
    cpu = sum(j["cpu_s"] for j in jobs)
    with open(a.out / f"{a.name}_chunks.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["chunk", "events", "seed", "wall_s", "cpu_s"])
        w.writeheader()
        for j in jobs:
            w.writerow({**j, "wall_s": f"{j['wall_s']:.2f}", "cpu_s": f"{j['cpu_s']:.2f}"})
    with open(a.out / f"{a.name}_summary.txt", "w") as fh:
        fh.write("=================================================\n"
                 "  MuonRock Run Summary (merged parallel run)\n"
                 "=================================================\n"
                 f"  Run name   : {a.name}\n"
                 f"  Source     : FILE  -> {a.source}\n"
                 f"  Events     : {n_total}\n"
                 f"  Seed       : {a.seed} + chunk index, {a.chunks} chunks\n"
                 f"  Slab half-XY: {a.xy or 2500:g} cm\n"
                 + (f"  Inclination: {a.incl:g} deg (depths are slant lengths)\n" if a.incl else "")
                 + f"  Muon plane crossings : {counts['muons']}\n"
                 f"  Sec hits             : {counts['secondaries']}\n"
                 f"  Stopped muons        : {counts['stopped']}\n"
                 f"  KE min     : {min(s.get('KE min', 1e30) for s in ke):.1f} GeV\n"
                 f"  KE max     : {max(s.get('KE max', 0) for s in ke):.1f} GeV\n"
                 f"  Elapsed    : {cpu:.2f} s   (CPU summed over the chunks)\n"
                 f"  Rate       : {n_total / cpu:.1f} evt/s   (per core)\n"
                 f"  CPU/evt    : {cpu / n_total * 1000:.3f} ms\n"
                 f"  Wall time  : {wall:.2f} s with {a.workers} worker processes\n"
                 "-------------------------------------------------\n")
    print(f"done: {n_total} events, CPU {cpu:.0f} s, wall {wall:.0f} s; merged files in {a.out}")
    if not a.keep_chunks:
        shutil.rmtree(chunk_root)
    return 0


if __name__ == "__main__":
    sys.exit(main())
