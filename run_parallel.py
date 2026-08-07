#!/usr/bin/env sage -python
"""
Parallel driver for the indecomposables database.

Run with::

    sage -python run_parallel.py --degree 3 --disc-max 100000 --workers 20

Design
------
* **Append-only shards.**  Each worker slot owns ``work/deg<N>/shard-<slot>.txt``
  and appends one row per completed field.  Workers never share a file, so no
  locking is needed and a crash costs at most the field in flight.
* **Resume is free.**  On start the driver reads every shard, collects the set
  of finished labels, and skips them.  Re-running the same command after any
  kind of interruption is always the right thing to do.
* **Nothing is silently lost.**  A field that finishes becomes a row in the
  worker's shard; one that does not becomes a line in ``failures.txt``.  Rows
  plus failures equals the scope, so every column in the data file can be an
  invariant of the field rather than carrying a "not computed" NULL.
* **Two timeouts.**  A soft one (SIGALRM inside the worker) which lets the
  worker record the timeout itself and carry on, and a hard one enforced by a
  watchdog in the parent, which kills and replaces a worker that has stopped
  responding.  The soft timeout handles the common case; the hard one exists
  because a signal cannot always interrupt a long-running PARI call.
* **Dynamic scheduling.**  Cost per field varies by orders of magnitude, so
  work is pulled from a queue rather than partitioned up front by line number.

Contract with the computation layer
-----------------------------------
The driver imports the following from ``--compute-module`` (default
``indecomposables.record``) and needs nothing else::

    build_context(coeffs, label)   -> ctx
    compute_record(ctx)            -> record
    failure_line(label, reason)    -> str, no newline
    to_row(record)                 -> str, no newline

Pass ``--self-test`` to run the orchestration against a synthetic workload with
no Sage at all; that mode exercises sharding, resume, timeouts and the watchdog.
"""

from __future__ import annotations

import argparse
import importlib
import multiprocessing as mp
import os
import queue
import random
import signal
import sys
import time
from dataclasses import dataclass
from pathlib import Path

SHARD_GLOB = "shard-*.txt"
POLL = 0.5


# ---------------------------------------------------------------------------
# Work items
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Item:
    label: str
    coeffs: tuple
    disc: int


def parse_input(path: Path, degree: int, disc_min: int, disc_max: int):
    """
    Read ``input/deg<N>.txt``.

    Expected columns ``lmfdb_index|coeffs|monogenic``, optionally with a ``disc``
    column.  Without ``disc`` the driver cannot filter by discriminant without
    constructing every field, which is exactly the waste we want to avoid, so
    that case is refused with a pointer at the fix.
    """
    with path.open(newline="") as f:
        lines = f.read().splitlines()          # tolerates CRLF
    if not lines:
        raise SystemExit(f"{path} is empty")

    header = [c.strip() for c in lines[0].split("|")]
    if "coeffs" not in header:
        raise SystemExit(f"{path}: expected a header line with a 'coeffs' column")
    if "disc" not in header:
        raise SystemExit(
            f"{path}: no 'disc' column.  Regenerate the input table from LMFDB "
            "with the discriminant included -- filtering and sharding both need "
            "it, and recomputing it per field is the thing this driver avoids.")

    idx = {name: i for i, name in enumerate(header)}
    items = []
    for line in lines[1:]:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = [c.strip() for c in line.split("|")]
        disc = abs(int(parts[idx["disc"]]))
        if disc > disc_max:
            break                               # input is sorted by |disc|
        if disc < disc_min:
            continue
        coeffs = tuple(int(c) for c in parts[idx["coeffs"]].split(","))
        index = parts[idx["lmfdb_index"]] if "lmfdb_index" in idx else "1"
        items.append(Item(f"{degree}.{degree}.{disc}.{index}", coeffs, disc))
    return items


def done_labels(work_dir: Path) -> set:
    """Labels already present in any shard.  This is the whole of resume."""
    seen = set()
    for shard in sorted(work_dir.glob(SHARD_GLOB)):
        with shard.open() as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    seen.add(line.split("|", 1)[0])
    return seen


# ---------------------------------------------------------------------------
# Worker
# ---------------------------------------------------------------------------

class _SoftTimeout(Exception):
    pass


def _alarm(signum, frame):
    raise _SoftTimeout()


def worker(slot, work_q, status_q, heartbeat, stop, args):
    """One worker slot.  Owns exactly one shard file, appends and flushes."""
    if args.self_test:
        compute = _SelfTestCompute()
    else:
        mod = importlib.import_module(args.compute_module)
        compute = mod

    base = Path(args.work_dir) / f"degree{args.degree}"
    shard = base / f"shard-{slot:03d}.txt"
    fails = base / "failures.txt"
    signal.signal(signal.SIGALRM, _alarm)
    signal.signal(signal.SIGINT, signal.SIG_IGN)      # parent owns Ctrl-C

    with shard.open("a", buffering=1) as out, fails.open("a", buffering=1) as bad:
        while not stop.is_set():
            try:
                item = work_q.get(timeout=POLL)
            except queue.Empty:
                continue
            if item is None:
                break

            heartbeat[slot] = (item.label, time.time())
            started = time.time()
            status, detail, record = "ok", "", None
            try:
                signal.setitimer(signal.ITIMER_REAL, args.timeout)
                ctx = compute.build_context(item.coeffs, item.label)
                record = compute.compute_record(ctx)
            except _SoftTimeout:
                status, detail = "timeout", f"exceeded {args.timeout}s"
            except Exception as exc:                  # noqa: BLE001 - recorded, not swallowed
                status = "error"
                detail = f"{type(exc).__name__}: {exc}"[:200]
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0)

            target = out if record is not None else bad
            target.write((compute.to_row(record) if record is not None
                          else compute.failure_line(item.label, detail)) + "\n")
            target.flush()
            if args.fsync:
                os.fsync(target.fileno())

            heartbeat[slot] = (None, time.time())
            status_q.put((item.label, status, time.time() - started, detail))


# ---------------------------------------------------------------------------
# Self-test workload (no Sage)
# ---------------------------------------------------------------------------

class _SelfTestCompute:
    """Synthetic workload with realistic cost skew, so the driver is testable."""

    def build_context(self, coeffs, label):
        rng = random.Random(label)
        return {"label": label, "cost": rng.lognormvariate(-2.5, 1.4), "rng": rng}

    def compute_record(self, ctx):
        if ctx["rng"].random() < 0.03:
            raise ValueError("synthetic failure")
        time.sleep(min(ctx["cost"], 30.0))
        n = ctx["rng"].randint(1, 40)
        return {"label": ctx["label"], "status": "ok", "n": n, "detail": ""}

    def failure_line(self, label, reason):
        return f"{label}|{reason}"

    def to_row(self, record):
        return "|".join(str(record[k]) for k in ("label", "status", "n", "detail"))


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--degree", type=int, required=True)
    p.add_argument("--disc-min", type=int, default=1)
    p.add_argument("--disc-max", type=int, required=True)
    p.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    p.add_argument("--input-dir", default="input")
    p.add_argument("--work-dir", default="work")
    p.add_argument("--compute-module", default="indecomposables.record")
    p.add_argument("--timeout", type=float, default=900.0,
                   help="soft per-field timeout in seconds")
    p.add_argument("--hard-timeout-factor", type=float, default=4.0,
                   help="kill a worker stuck this many times past --timeout")
    p.add_argument("--fsync", action="store_true",
                   help="fsync after every row (slower, survives power loss)")
    p.add_argument("--no-resume", action="store_true")
    p.add_argument("--limit", type=int, default=None, help="cap work items (testing)")
    p.add_argument("--shuffle", action="store_true",
                   help="shuffle work order; useful for an unbiased partial run")
    p.add_argument("--seed", type=int, default=20260805)
    p.add_argument("--progress-every", type=float, default=15.0)
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    mp.set_start_method("fork", force=True)

    work_dir = Path(args.work_dir) / f"degree{args.degree}"
    work_dir.mkdir(parents=True, exist_ok=True)

    # Fail fast in the parent: otherwise every worker dies separately on the
    # same import and the only symptom is "0 fields reached".
    if not args.self_test:
        try:
            mod = importlib.import_module(args.compute_module)
        except ImportError as exc:
            raise SystemExit(
                f"cannot import --compute-module {args.compute_module!r}: {exc}\n"
                "install the package with `sage -pip install -e .`, or pass "
                "--self-test to exercise the driver without it")
        missing = [f for f in ("build_context", "compute_record",
                               "failure_line", "to_row")
                   if not hasattr(mod, f)]
        if missing:
            raise SystemExit(
                f"{args.compute_module} is missing {', '.join(missing)}; "
                "see the contract in this file's docstring")

    if args.self_test:
        items = [Item(f"{args.degree}.{args.degree}.{1000 + i}.1", (1, -2, -1), 1000 + i)
                 for i in range(args.limit or 60)]
    else:
        src = Path(args.input_dir) / f"degree{args.degree}.txt"
        items = parse_input(src, args.degree, args.disc_min, args.disc_max)

    total_in_scope = len(items)
    skipped = 0
    if not args.no_resume:
        finished = done_labels(work_dir)
        before = len(items)
        items = [it for it in items if it.label not in finished]
        skipped = before - len(items)

    if args.shuffle:
        random.Random(args.seed).shuffle(items)
    if args.limit:
        items = items[:args.limit]

    print(f"degree {args.degree}, |disc| in [{args.disc_min}, {args.disc_max}]")
    print(f"  in scope: {total_in_scope}   already done: {skipped}   to do: {len(items)}")
    print(f"  workers: {args.workers}   soft timeout: {args.timeout}s   "
          f"hard: {args.timeout * args.hard_timeout_factor}s")
    print(f"  shards:  {work_dir}/shard-NNN.txt")
    if args.dry_run or not items:
        print("nothing to do" if not items else "dry run; exiting")
        return 0

    work_q: mp.Queue = mp.Queue()
    status_q: mp.Queue = mp.Queue()
    for it in items:
        work_q.put(it)
    for _ in range(args.workers):
        work_q.put(None)

    manager = mp.Manager()
    heartbeat = manager.dict({i: (None, time.time()) for i in range(args.workers)})
    stop = mp.Event()

    def spawn(slot):
        p = mp.Process(target=worker, name=f"slot-{slot}",
                       args=(slot, work_q, status_q, heartbeat, stop, args),
                       daemon=True)
        p.start()
        return p

    procs = {slot: spawn(slot) for slot in range(args.workers)}

    interrupted = {"flag": False}

    def on_sigint(signum, frame):
        if interrupted["flag"]:
            print("\nsecond interrupt: exiting now", flush=True)
            raise SystemExit(130)
        interrupted["flag"] = True
        stop.set()
        print("\ninterrupt: letting workers finish the field in flight "
              "(Ctrl-C again to abandon it)", flush=True)

    signal.signal(signal.SIGINT, on_sigint)

    counts = {"ok": 0, "timeout": 0, "error": 0, "killed": 0}
    slowest = []
    t0 = time.time()
    last_report = t0
    hard = args.timeout * args.hard_timeout_factor

    try:
        while any(p.is_alive() for p in procs.values()):
            # drain status
            drained = False
            while True:
                try:
                    label, status, elapsed, detail = status_q.get_nowait()
                except queue.Empty:
                    break
                drained = True
                counts[status] = counts.get(status, 0) + 1
                slowest.append((elapsed, label))
                if status != "ok":
                    print(f"  [{status}] {label}: {detail}", flush=True)

            # watchdog: replace anything wedged past the hard timeout
            now = time.time()
            for slot, proc in list(procs.items()):
                label, since = heartbeat.get(slot, (None, now))
                if label is not None and now - since > hard and proc.is_alive():
                    print(f"  [killed] slot {slot} wedged on {label} "
                          f"({now - since:.0f}s); restarting", flush=True)
                    proc.terminate()
                    proc.join(10)
                    if proc.is_alive():
                        os.kill(proc.pid, signal.SIGKILL)
                        proc.join(5)
                    counts["killed"] += 1
                    _record_kill(work_dir, label, hard)
                    heartbeat[slot] = (None, time.time())
                    if not stop.is_set():
                        procs[slot] = spawn(slot)
                elif not proc.is_alive() and proc.exitcode not in (0, None) \
                        and not stop.is_set():
                    print(f"  [killed] slot {slot} died "
                          f"(exit {proc.exitcode}); restarting", flush=True)
                    procs[slot] = spawn(slot)

            if now - last_report >= args.progress_every:
                _report(counts, len(items), t0, now)
                last_report = now
            if not drained:
                time.sleep(POLL)
    finally:
        stop.set()
        for p in procs.values():
            p.join(timeout=30)
            if p.is_alive():
                p.terminate()
        # Workers can exit with status messages still in the pipe; a row on disk
        # that we failed to count would look like an unfinished field.
        deadline = time.time() + 5
        while time.time() < deadline:
            try:
                label, status, elapsed, detail = status_q.get(timeout=0.2)
            except queue.Empty:
                break
            counts[status] = counts.get(status, 0) + 1
            slowest.append((elapsed, label))

    _report(counts, len(items), t0, time.time(), final=True)
    if slowest:
        slowest.sort(reverse=True)
        print("  slowest fields:")
        for elapsed, label in slowest[:5]:
            print(f"    {elapsed:8.1f}s  {label}")
    done = sum(counts.values())
    if done < len(items):
        print(f"\n{len(items) - done} field(s) not reached; re-run the same "
              "command to resume.")
    return 0


def _record_kill(work_dir: Path, label: str, hard: float):
    """A hard-killed worker cannot write its own line, so the parent does it."""
    with (work_dir / "hard-timeouts.txt").open("a", buffering=1) as f:
        f.write(f"{label}|killed after {hard:.0f}s\n")


def _report(counts, total, t0, now, final=False):
    done = sum(counts.values())
    rate = done / max(now - t0, 1e-9)
    eta = (total - done) / rate if rate > 0 and not final else 0
    bits = "  ".join(f"{k}={v}" for k, v in sorted(counts.items()) if v)
    tag = "final" if final else "progress"
    line = (f"[{tag}] {done}/{total}  {bits or 'ok=0'}  "
            f"{rate * 3600:.0f}/hr  elapsed {now - t0:.0f}s")
    if not final and rate > 0:
        line += f"  eta {eta / 60:.0f}m"
    print(line, flush=True)


if __name__ == "__main__":
    sys.exit(main())
