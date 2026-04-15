#!/usr/bin/env python3

import argparse
import csv
import datetime as dt
import os
import re
import subprocess
import sys
import tempfile
from typing import Dict, List


TIME_PATTERNS = {
    "elapsed_raw": re.compile(r"^\s*Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(.+?)\s*$"),
    "user_seconds": re.compile(r"^\s*User time \(seconds\):\s*([0-9]*\.?[0-9]+)\s*$"),
    "system_seconds": re.compile(r"^\s*System time \(seconds\):\s*([0-9]*\.?[0-9]+)\s*$"),
    "max_rss_kb": re.compile(r"^\s*Maximum resident set size \(kbytes\):\s*([0-9]+)\s*$"),
}


def parse_elapsed_to_seconds(raw: str) -> float:
    # Handles formats like "0:03.12", "12:34.56", "1:02:03.45".
    parts = raw.strip().split(":")
    if len(parts) == 2:
        minutes = int(parts[0])
        seconds = float(parts[1])
        return minutes * 60 + seconds
    if len(parts) == 3:
        hours = int(parts[0])
        minutes = int(parts[1])
        seconds = float(parts[2])
        return hours * 3600 + minutes * 60 + seconds
    raise ValueError(f"Unrecognized elapsed time format: {raw!r}")


def parse_time_v(stderr_text: str) -> Dict[str, float]:
    parsed: Dict[str, float] = {}
    elapsed_raw = None
    for line in stderr_text.splitlines():
        for key, pattern in TIME_PATTERNS.items():
            m = pattern.match(line)
            if not m:
                continue
            if key == "elapsed_raw":
                elapsed_raw = m.group(1)
            elif key == "max_rss_kb":
                parsed[key] = int(m.group(1))
            else:
                parsed[key] = float(m.group(1))
    if elapsed_raw is not None:
        parsed["wall_seconds"] = parse_elapsed_to_seconds(elapsed_raw)
    return parsed


def run_one(
    vg_bin: str,
    mode: str,
    input_path: str,
    threads: int,
    keep_output: bool,
) -> Dict[str, object]:
    if mode == "gfa":
        mode_flag = "-g"
    elif mode == "gfaz":
        mode_flag = "-g"
    else:
        raise ValueError(f"Unsupported mode: {mode}")

    output_path = tempfile.mkstemp(prefix=f"bench_{mode}_t{threads}_", suffix=".vg")[1]
    cmd = ["/usr/bin/time", "-v", vg_bin, "convert", mode_flag, input_path, "-p", "-t", str(threads)]

    try:
        with open(output_path, "wb") as out_f:
            proc = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
        metrics = parse_time_v(proc.stderr or "")
        result: Dict[str, object] = {
            "mode": mode,
            "threads": threads,
            "input_path": input_path,
            "exit_code": proc.returncode,
            "wall_seconds": metrics.get("wall_seconds"),
            "user_seconds": metrics.get("user_seconds"),
            "system_seconds": metrics.get("system_seconds"),
            "max_rss_kb": metrics.get("max_rss_kb"),
            "output_path": output_path if keep_output else "",
        }
        if proc.returncode != 0:
            result["error"] = (proc.stderr or "").strip().replace("\n", " | ")
        else:
            result["error"] = ""
        return result
    finally:
        if not keep_output and os.path.exists(output_path):
            os.remove(output_path)


def parse_threads(raw: str) -> List[int]:
    values = []
    for token in raw.split(","):
        token = token.strip()
        if not token:
            continue
        values.append(int(token))
    if not values:
        raise ValueError("No thread values provided")
    return values


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Benchmark vg convert performance for GFA (-g) vs GFAZ (-z) with /usr/bin/time -v."
    )
    parser.add_argument("--vg", default="./bin/vg", help="Path to vg binary (default: ./bin/vg)")
    parser.add_argument("--gfa", required=True, help="Input GFA path")
    parser.add_argument("--gfaz", required=True, help="Input GFAZ path")
    parser.add_argument(
        "--threads",
        default="1,4,8,16",
        help="Comma-separated thread counts to test (default: 1,4,8,16)",
    )
    parser.add_argument("--csv", required=True, help="Output CSV path")
    parser.add_argument(
        "--keep-output",
        action="store_true",
        help="Keep generated .vg outputs (off by default)",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=1,
        help="Repeats per mode/thread combination (default: 1)",
    )
    args = parser.parse_args()

    vg_bin = os.path.abspath(args.vg)
    gfa_path = os.path.abspath(args.gfa)
    gfaz_path = os.path.abspath(args.gfaz)
    csv_path = os.path.abspath(args.csv)

    if not os.path.exists(vg_bin):
        print(f"error: vg binary not found: {vg_bin}", file=sys.stderr)
        return 2
    if not os.path.exists(gfa_path):
        print(f"error: gfa input not found: {gfa_path}", file=sys.stderr)
        return 2
    if not os.path.exists(gfaz_path):
        print(f"error: gfaz input not found: {gfaz_path}", file=sys.stderr)
        return 2
    if args.repeats < 1:
        print("error: --repeats must be >= 1", file=sys.stderr)
        return 2

    try:
        threads = parse_threads(args.threads)
    except Exception as e:
        print(f"error: invalid --threads value: {e}", file=sys.stderr)
        return 2

    os.makedirs(os.path.dirname(csv_path) or ".", exist_ok=True)

    rows: List[Dict[str, object]] = []
    start_ts = dt.datetime.now().isoformat(timespec="seconds")
    print(f"[{start_ts}] benchmarking started")

    for repeat in range(1, args.repeats + 1):
        for t in threads:
            for mode, path in (("gfa", gfa_path), ("gfaz", gfaz_path)):
                print(f"running mode={mode} threads={t} repeat={repeat}")
                row = run_one(
                    vg_bin=vg_bin,
                    mode=mode,
                    input_path=path,
                    threads=t,
                    keep_output=args.keep_output,
                )
                row["repeat"] = repeat
                row["timestamp"] = dt.datetime.now().isoformat(timespec="seconds")
                rows.append(row)
                if row["exit_code"] != 0:
                    print(
                        f"warning: run failed mode={mode} threads={t} repeat={repeat} "
                        f"exit={row['exit_code']}",
                        file=sys.stderr,
                    )

    fieldnames = [
        "timestamp",
        "repeat",
        "mode",
        "threads",
        "input_path",
        "exit_code",
        "wall_seconds",
        "user_seconds",
        "system_seconds",
        "max_rss_kb",
        "output_path",
        "error",
    ]
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    end_ts = dt.datetime.now().isoformat(timespec="seconds")
    print(f"[{end_ts}] benchmarking complete, wrote {len(rows)} rows to {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
