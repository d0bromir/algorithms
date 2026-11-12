#!/usr/bin/env python3
"""certa-qc against fastp, then certa map: wall and CPU time.

The data are those of tests/qc_equivalence.py, larger: a 2 Mbp simulated
genome and 400 000 pairs of 150 bases (gzip input) with adapter read-through,
poly-G tails and poor ends. Each configuration runs once per round, and the
rounds interleave the configurations. It prints every run and the medians of
wall time and CPU time (user + system, of the process and its threads). Every
run's SAM must be identical.

  two-step:  fastp -w W -o c1 -O c2 (gzip or plain)  then  certa map -t T c1 c2
  one pass:  certa-qc -t T --fastp -w W

Usage:
  qc_speed.py --certa BUILD/certa --certa-qc BUILD/certa-qc --fastp FASTP/fastp
              [--workdir DIR] [--pairs 400000] [--rounds 3]
"""
import argparse, hashlib, os, resource, statistics, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tests"))
import qc_equivalence  # noqa: E402  (the simulated data)

# (name, kind, fastp -w, certa -t, intermediate file extension)
CONFIGS = [
    ("two-step gzip", "two", 4, 4, ".fq.gz"),
    ("two-step plain", "two", 4, 4, ".fq"),
    ("certa-qc", "one", 4, 4, None),
    ("certa-qc", "one", 1, 4, None),
    ("certa-qc", "one", 1, 3, None),
    ("two-step plain", "two", 1, 1, ".fq"),
    ("certa-qc", "one", 1, 1, None),
]


def timed(cmd):
    r0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    t = time.time()
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
    return time.time() - t, (r1.ru_utime - r0.ru_utime) + (r1.ru_stime - r0.ru_stime)


def run(a, cfg):
    name, kind, w, t, ext = cfg
    o = os.path.join(a.workdir, f"run_{kind}_{w}_{t}_{(ext or '').strip('.').replace('.', '_')}")
    os.makedirs(o, exist_ok=True)
    idx = os.path.join(a.workdir, "genome.cidx")
    inputs = ["-i", os.path.join(a.workdir, "in1.fq.gz"), "-I", os.path.join(a.workdir, "in2.fq.gz")]
    reports = ["-j", o + "/fastp.json", "-h", o + "/fastp.html"]
    mapped = ["-t", str(t), "-o", o + "/out.sam", "-u", o + "/u1.fq", "-U", o + "/u2.fq"]
    steps = ""
    if kind == "one":
        wall, cpu = timed([a.certa_qc, idx] + mapped + ["--fastp"] + inputs + ["-w", str(w)] + reports)
    else:
        c = [o + "/c1" + ext, o + "/c2" + ext]
        w1, c1 = timed([a.fastp] + inputs + ["-o", c[0], "-O", c[1], "-w", str(w)] + reports)
        w2, c2 = timed([a.certa, "map", idx] + c + mapped)
        wall, cpu = w1 + w2, c1 + c2
        steps = f"  (fastp wall {w1:.2f} s cpu {c1:.2f} s; certa map wall {w2:.2f} s cpu {c2:.2f} s)"
    with open(o + "/out.sam", "rb") as f:
        return wall, cpu, hashlib.md5(f.read()).hexdigest(), steps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--certa", required=True)
    ap.add_argument("--certa-qc", dest="certa_qc", required=True)
    ap.add_argument("--fastp", required=True)
    ap.add_argument("--workdir", default="qc_speed_work")
    ap.add_argument("--pairs", type=int, default=400000)
    ap.add_argument("--rounds", type=int, default=3)
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)
    qc_equivalence.make_data(a.workdir, a.certa, pairs=a.pairs, seed=8)
    print(f"{a.pairs} pairs, {os.cpu_count()} CPUs, {a.rounds} rounds", flush=True)
    res = {cfg: [] for cfg in CONFIGS}
    for r in range(a.rounds):
        for cfg in CONFIGS:
            res[cfg].append(run(a, cfg))
            wall, cpu, _, steps = res[cfg][-1]
            print(f"round {r + 1}  {cfg[0]:15} -w {cfg[2]} -t {cfg[3]}  wall {wall:6.2f} s  cpu {cpu:6.2f} s{steps}",
                  flush=True)
    sams = {x[2] for runs in res.values() for x in runs}
    print(f"distinct SAM files: {len(sams)}")
    print(f"{'configuration':15} {'-w':>3} {'-t':>3} {'wall (median)':>14} {'cpu (median)':>13} {'wall range':>14}")
    for cfg, runs in res.items():
        walls, cpus = [x[0] for x in runs], [x[1] for x in runs]
        print(f"{cfg[0]:15} {cfg[2]:3} {cfg[3]:3} {statistics.median(walls):12.2f} s {statistics.median(cpus):11.2f} s"
              f" {min(walls):6.2f}-{max(walls):.2f} s")
    return 0 if len(sams) == 1 else 1


if __name__ == "__main__":
    sys.exit(main())
