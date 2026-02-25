#!/usr/bin/env python3
"""
Generic HELIOS job distributor.

Features
--------
- Modes:
    0 : isolated homogeneous      (old run_iso_hom.sh)
    1 : periodic homogeneous      (old run_per_hom.sh)
    2 : isolated layered          (old run_iso_layered.sh)
- Node selection:
    -N nodeA,nodeB,...            explicit list (works anywhere)
    -C NUM                        first NUM nodes from SLURM_NODELIST (Slurm only)
- Optional local 'prepare' step before distributed solve/post.
- Safe SSH:
    * no '-tt', so terminal is not corrupted
    * stdin for ssh is /dev/null (remote jobs can't mess with TTY)
- Even job partition across any number of nodes.
- Layered 'post' can be retried with RETRIES env (default 4).
"""

import argparse
import glob
import os
import shlex
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import List, Sequence

# ------------ Defaults (overridable via env) ----------------

DEFAULT_ACTIVATE = (
    'if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then '
    'source "$HOME/miniconda3/etc/profile.d/conda.sh"; '
    'elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then '
    'source "$HOME/anaconda3/etc/profile.d/conda.sh"; '
    'elif [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then '
    'source "/opt/conda/etc/profile.d/conda.sh"; '
    'else echo "[activate] conda.sh not found" >&2; fi; '
    'conda activate Helios'
)

ACTIVATE = os.environ.get("ACTIVATE", DEFAULT_ACTIVATE)
ROOT_ENV = os.environ.get("ROOT", "")
PY_ENV = os.environ.get("PY", "run_sie.py")

# Safer ssh: no -tt (no forced TTY), overridable with SSH_BIN env
SSH_BIN_ENV = os.environ.get(
    "SSH_BIN",
    "ssh -o BatchMode=yes -o LogLevel=ERROR"
)

MAX_ATTEMPTS = int(os.environ.get("RETRIES", "4"))  # for layered post
HELIOS_CANCEL_FILE = os.environ.get(
    "HELIOS_CANCEL_FILE",
    f"/tmp/helios_cancel_{os.getpid()}"
)


# ------------ Utility helpers ----------------

def build_cmd(args: Sequence[str]) -> str:
    """Shell-escape argv list into a single command string."""
    return " ".join(shlex.quote(a) for a in args)


def resolve_remote_path(p: str, root: str) -> str:
    """Absolute path stays absolute; otherwise root/p or p."""
    if p.startswith("/"):
        return p
    if root:
        return f"{root.rstrip('/')}/{p}"
    return p


def expand_sel(sel: str) -> List[int]:
    """Expand '1,3,5-8' -> [1,3,5,6,7,8]."""
    out: List[int] = []
    if not sel:
        return out
    for tok in sel.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if "-" in tok:
            a_str, b_str = tok.split("-", 1)
            if not (a_str.lstrip("-").isdigit() and b_str.lstrip("-").isdigit()):
                continue
            a, b = int(a_str), int(b_str)
            if a > b:
                a, b = b, a
            out.extend(range(a, b + 1))
        else:
            if tok.lstrip("-").isdigit():
                out.append(int(tok))
    return out


def select_nodes_from_allocation(need: int, ssh_bin: str) -> str:
    """Use SLURM_NODELIST and scontrol show hostnames to pick 'need' nodes."""
    slurm_nodelist = os.environ.get("SLURM_NODELIST", "")
    if not slurm_nodelist:
        print("ERROR: -C/--count requires a Slurm allocation (SLURM_NODELIST not set).",
              file=sys.stderr)
        sys.exit(1)

    try:
        hosts_text = subprocess.check_output(
            ["scontrol", "show", "hostnames", slurm_nodelist],
            text=True
        )
    except Exception as e:
        print(f"ERROR: failed to expand SLURM_NODELIST via scontrol: {e}",
              file=sys.stderr)
        sys.exit(1)

    hosts = [h for h in hosts_text.splitlines() if h.strip()]
    if not hosts:
        print("ERROR: allocation has no hosts", file=sys.stderr)
        sys.exit(1)
    if need > len(hosts):
        print(f"ERROR: requested {need} nodes, but allocation has only {len(hosts)}.",
              file=sys.stderr)
        sys.exit(1)

    ssh_parts = shlex.split(ssh_bin)
    picked: List[str] = []
    for h in hosts[:need]:
        cmd = (
            ssh_parts
            + [
                "-o", "ConnectTimeout=2",
                "-o", "StrictHostKeyChecking=no",
                h,
                "echo ok",
            ]
        )
        try:
            subprocess.run(
                cmd,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                stdin=subprocess.DEVNULL,
                timeout=2,
                check=True,
            )
            picked.append(h)
        except Exception:
            print(f"WARN: {h} not reachable via SSH; continuing...", file=sys.stderr)

    if len(picked) != need:
        print(f"ERROR: only {len(picked)}/{need} allocation nodes reachable via ssh.",
              file=sys.stderr)
        sys.exit(1)

    return ",".join(picked)


# ------------ CLI parsing ----------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generic HELIOS distributed runner (0=iso-hom,1=per-hom,2=iso-layered)."
    )
    parser.add_argument(
        "--mode",
        required=True,
        choices=["0", "1", "2"],
        help="0: isolated homogeneous, 1: periodic homogeneous, 2: isolated layered.",
    )
    parser.add_argument("-s", "--sim", required=True, help="Simulation name (SIM).")
    parser.add_argument(
        "-N", "--nodes",
        help="Comma-separated node list: nodeA,nodeB,... (works on any system)."
    )
    parser.add_argument(
        "-C", "--count",
        type=int,
        dest="node_count",
        help="Use COUNT nodes from SLURM_NODELIST (Slurm only)."
    )
    parser.add_argument(
        "--prepare",
        action="store_true",
        help="Run local 'prepare' once before distributed solve/post.",
    )
    parser.add_argument(
        "-r", "--root",
        help="Project root (remote cwd); defaults to $ROOT or current directory."
    )
    parser.add_argument(
        "-R", "--runner",
        dest="runner",
        help="Path to run_sie.py (absolute or relative to ROOT). Defaults to $PY or 'run_sie.py'.",
    )
    parser.add_argument(
        "--ssh",
        dest="ssh_bin",
        help="SSH command template (default: 'ssh -o BatchMode=yes -o LogLevel=ERROR')."
    )
    parser.add_argument(
        "subcmd",
        nargs="?",
        default="solve",
        help="run_sie.py subcommand: solve/post/... (default: solve).",
    )
    parser.add_argument(
        "run_sie_args",
        nargs=argparse.REMAINDER,
        help="Arguments forwarded to run_sie.py after subcmd.",
    )

    args = parser.parse_args()
    args.mode_num = int(args.mode)
    return args


# ------------ Main logic ----------------

def main() -> int:
    args = parse_args()

    sim = args.sim
    root = args.root if args.root is not None else ROOT_ENV
    runner = args.runner if args.runner is not None else PY_ENV
    ssh_bin = args.ssh_bin if args.ssh_bin is not None else SSH_BIN_ENV
    ssh_parts = shlex.split(ssh_bin)
    mode = args.mode_num

    # Node selection: N explicit, else C via Slurm
    nodes_csv = args.nodes or ""
    if not nodes_csv and args.node_count:
        nodes_csv = select_nodes_from_allocation(args.node_count, ssh_bin)
        print(f"[info] Using allocation nodes: {nodes_csv}")

    if not nodes_csv:
        print("ERROR: provide -N nodeA,nodeB or use -C NUM inside Slurm.",
              file=sys.stderr)
        return 2

    nodes = [n for n in nodes_csv.split(",") if n.strip()]
    if not nodes:
        print("ERROR: empty node list after parsing -N/--nodes.",
              file=sys.stderr)
        return 2

    subcmd = args.subcmd
    rest = list(args.run_sie_args)  # extra args to run_sie.py

    # Resolve remote runner path (run_sie.py)
    remote_py = resolve_remote_path(runner, root)

    # Optional local prepare once
    if args.prepare:
        print(f"[local] prepare -> {runner} prepare {sim} --mode {mode} {' '.join(rest)}")
        prepare_cmd = ["python3", runner, "prepare", sim, "--mode", str(mode)] + rest
        try:
            if root:
                subprocess.run(prepare_cmd, cwd=root, check=True)
            else:
                subprocess.run(prepare_cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"[local] prepare failed with rc={e.returncode}",
                  file=sys.stderr)
            return e.returncode

    # Collect jobs from sim_res/SIM/jobs
    base_root = root if root else os.getcwd()
    res_dir = os.path.join(base_root, "sim_res", sim)
    jdir = os.path.join(res_dir, "jobs")

    if not os.path.isdir(jdir):
        print(f"Error: jobs dir not found: {jdir}", file=sys.stderr)
        return 2

    job_files = sorted(
        glob.glob(os.path.join(jdir, "job.*.xml")),
        key=lambda f: int(os.path.basename(f).split(".")[1]),
    )
    if not job_files:
        print(f"Error: no job.*.xml in {jdir}", file=sys.stderr)
        return 2

    jidx: List[int] = []
    for f in job_files:
        b = os.path.basename(f)
        try:
            i = int(b.split(".")[1])
            jidx.append(i)
        except (IndexError, ValueError):
            continue

    if not jidx:
        print("Error: could not parse any job indices from job.*.xml",
              file=sys.stderr)
        return 2

    # Look for -j / --jobs in REST
    user_sel = ""
    for i, tok in enumerate(rest):
        if tok in ("-j", "--jobs") and i + 1 < len(rest):
            user_sel = rest[i + 1]
            break

    # Build final WANT list
    if user_sel:
        exp = expand_sel(user_sel)
        have = {i: True for i in jidx}
        want = [i for i in exp if i in have]
    else:
        want = list(jidx)

    if not want:
        print("Error: selection produced 0 jobs", file=sys.stderr)
        return 2

    # Strip -j/--jobs from REST (we add -j per node)
    rest_noj: List[str] = []
    skip_next = False
    for tok in rest:
        if skip_next:
            skip_next = False
            continue
        if tok in ("-j", "--jobs"):
            skip_next = True
            continue
        rest_noj.append(tok)

    def has_flag(flag: str) -> bool:
        return flag in rest_noj

    # Preflight: ensure run_sie.py exists on each node
    for node in nodes:
        if root:
            remote_cmd = f"cd {shlex.quote(root)} && test -f {shlex.quote(remote_py)}"
        else:
            remote_cmd = f"test -f {shlex.quote(remote_py)}"
        try:
            subprocess.run(
                ssh_parts + [node, "bash", "-lc", remote_cmd],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                stdin=subprocess.DEVNULL,
                check=True,
            )
        except subprocess.CalledProcessError:
            print(f"[preflight] {node}: runner not found at {remote_py}",
                  file=sys.stderr)
            return 2

    # Partition jobs into contiguous chunks
    nj = len(want)
    nn = len(nodes)
    chunk = (nj + nn - 1) // nn

    slices = []
    for n in range(nn):
        start = n * chunk
        end = min(start + chunk, nj)
        if start >= nj:
            print(f"[{nodes[n]}] no jobs")
            continue
        slice_indices = want[start:end]
        jstr = ",".join(str(i) for i in slice_indices)
        print(f"[dispatch] {nodes[n]} -> {jstr}")
        slices.append((nodes[n], jstr))

    # Shared bookkeeping
    procs_lock = threading.Lock()
    procs: List[subprocess.Popen] = []
    done_counter = 0
    total_slices = len(slices)
    cancel_event = threading.Event()
    cleanup_once = False
    rc_holder: List[int] = []

    def register_proc(p: subprocess.Popen) -> None:
        with procs_lock:
            procs.append(p)

    def cleanup_handler(signum, frame):
        nonlocal cleanup_once
        if cleanup_once:
            return
        cleanup_once = True
        cancel_event.set()
        try:
            Path(HELIOS_CANCEL_FILE).touch()
        except Exception:
            pass
        print("[cleanup] Caught signal; terminating remote sessions...",
              file=sys.stderr)
        with procs_lock:
            current = list(procs)
        for p in current:
            try:
                p.send_signal(signal.SIGINT)
            except Exception:
                pass
        time.sleep(0.2)
        for p in current:
            try:
                p.terminate()
            except Exception:
                pass

    signal.signal(signal.SIGINT, cleanup_handler)
    signal.signal(signal.SIGTERM, cleanup_handler)
    try:
        signal.signal(signal.SIGHUP, cleanup_handler)
    except (AttributeError, ValueError):
        pass

    def progress_done():
        nonlocal done_counter
        done_counter += 1
        print(f"[dist-progress] done {done_counter}")

    def make_remote_body(jstr: str) -> str:
        args_list: List[str] = ["env", "SIE_FORCE_BAR=1", "SIE_NO_SPINNER=1"]
        if mode == 2:
            args_list.append("PYTHONUNBUFFERED=1")
        args_list += [
            "python3",
            remote_py,
            subcmd,
            sim,
            "--mode",
            str(mode),
        ]
        if not has_flag("--overwrite") and not has_flag("--keep"):
            args_list.append("--keep")
        args_list += rest_noj
        args_list += ["-j", jstr]
        cmdline = build_cmd(args_list)
        remote_trap = "trap 'kill -- -$$' INT TERM HUP;"
        if root:
            body = f"{remote_trap} {ACTIVATE} && cd {shlex.quote(root)} && exec {cmdline}"
        else:
            body = f"{remote_trap} {ACTIVATE} && exec {cmdline}"
        return body

    print(f"[dist-progress] total {total_slices}")

    # Layered post gets retries; others are one-shot
    if mode == 2 and subcmd == "post":
        threads = []

        def run_slice_post(node: str, jstr: str):
            attempt = 1
            last_rc = 1
            while attempt <= MAX_ATTEMPTS:
                if cancel_event.is_set():
                    print(f"[dispatch] {node} slice {jstr} cancelled; aborting retries")
                    last_rc = 130
                    break
                body = make_remote_body(jstr)
                print(f"[dispatch] {node} attempt {attempt}/{MAX_ATTEMPTS} cmd: {body}")
                p = subprocess.Popen(
                    ssh_parts + [node, "bash", "-lc", body],
                    stdin=subprocess.DEVNULL,
                )
                register_proc(p)
                try:
                    rc = p.wait()
                except KeyboardInterrupt:
                    rc = 130
                last_rc = rc

                if cancel_event.is_set():
                    print(f"[dispatch] {node} slice {jstr} cancelled during attempt {attempt}")
                    last_rc = 130
                    break

                if rc == 0:
                    print(f"[dispatch] {node} slice {jstr} finished OK on attempt {attempt}")
                    break

                print(f"[dispatch] {node} slice {jstr} failed (rc={rc}) on attempt {attempt}")
                attempt += 1
                if attempt <= MAX_ATTEMPTS:
                    time.sleep(2)
                    print(f"[dispatch] {node} retrying slice {jstr}...")

            if last_rc not in (0, 130):
                rc_holder.append(1)
            progress_done()

        for node, jstr in slices:
            t = threading.Thread(target=run_slice_post, args=(node, jstr), daemon=True)
            threads.append(t)
            t.start()
        for t in threads:
            t.join()

    else:
        # Simple one-shot per slice
        for node, jstr in slices:
            body = make_remote_body(jstr)
            print(f"[dispatch] {node} cmd: {body}")
            p = subprocess.Popen(
                ssh_parts + [node, "bash", "-lc", body],
                stdin=subprocess.DEVNULL,
            )
            register_proc(p)

        with procs_lock:
            local_procs = list(procs)
        for p in local_procs:
            try:
                rc = p.wait()
            except KeyboardInterrupt:
                rc = 130
            if rc not in (0, 130):
                rc_holder.append(1)
            progress_done()

    return 1 if rc_holder else 0


if __name__ == "__main__":
    sys.exit(main())
