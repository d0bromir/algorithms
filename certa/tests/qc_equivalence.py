#!/usr/bin/env python3
"""certa-qc against the two-step pipeline it replaces.

For each case it runs

  two-step:  fastp [opts] -o c1.fq [-O c2.fq]  then  certa map ref.cidx c1.fq [c2.fq] [map opts]
  one pass:  certa-qc ref.cidx [map opts] --fastp [opts]

and compares, byte for byte: the SAM, the uncertified FASTQ files, and the
files fastp writes for the reads it drops (--failed_out, --unpaired1/2, after
decompression). It also compares fastp's JSON report (without its volatile
keys) and certa's --stats JSON (without "seconds").

The data are simulated: a 2 Mbp genome and 40 000 pairs of 150 bases with
NovaSeq X qualities. Some reads get adapter read-through, poly-G tails, poor
3' ends or poor qualities throughout, so fastp trims and drops reads.

Usage:
  qc_equivalence.py --certa BUILD/certa --certa-qc BUILD/certa-qc --fastp FASTP/fastp
                    [--workdir DIR] [--jobs 2] [--only TEXT]
Exit status 0 if every case is identical, 1 otherwise.
"""
import argparse, gzip, hashlib, json, os, random, shutil, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
SIMULATE = os.path.join(HERE, "..", "tools", "simulate.py")
ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
FASTP_VOLATILE = {"command", "fastp_version", "time_consuming"}

# fastp option sets: name -> (options, paired only). "{D}" is the case's directory.
FASTP_SETS = {
    "default": ([], False),
    "cut_and_length": (["--cut_front", "--cut_tail", "-q", "20", "-l", "50"], False),
    "polyg_dedup_failed": (["--trim_poly_g", "--dedup", "--failed_out", "{D}/failed.fq.gz"], False),
    "front_tail_umi": (["-f", "3", "-t", "2", "--umi", "--umi_loc", "read1", "--umi_len", "6"], False),
    "detect_pe_correction": (["--detect_adapter_for_pe", "--correction"], True),
    "unpaired_failed": (["--unpaired1", "{D}/U1.fq", "--unpaired2", "{D}/U2.fq", "--failed_out", "{D}/failed.fq"], True),
}
# certa map modes: name -> options. Interval mode writes every read to SAM (no -u/-U).
CERTA_MODES = {
    "default": [],
    "intervals": ["--intervals"],
    "full_cert": ["--full-cert"],
    "bio": ["--bio"],
}
# (fastp -w, certa -t). The plain input is large enough for fastp to use all
# four workers (fastp lowers -w for small inputs).
THREADS = [(1, 1), (4, 4)]
# fastp's --dedup keeps the first copy of a duplicate it processes. With more
# than one worker that depends on the order of the packs, so fastp's own output
# changes from run to run (fastp-gpu at 73a3009 too). These sets run at -w 1.
ONE_WORKER_SETS = {"polyg_dedup_failed"}


def run(cmd, log):
    with open(log, "wb") as f:
        p = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=f)
    return p.returncode


def content_md5(path):
    opener = gzip.open if path.endswith(".gz") else open
    h = hashlib.md5()
    with opener(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def make_data(wd, certa):
    """Simulated genome, index and reads with fastp work in them."""
    g = os.path.join(wd, "genome.fa")
    if not os.path.exists(os.path.join(wd, "genome.cidx")):
        subprocess.run([sys.executable, SIMULATE, "genome", g, "--length", "2000000", "--contigs", "3",
                        "--seed", "5"], check=True, stdout=subprocess.DEVNULL)
        subprocess.run([sys.executable, SIMULATE, "reads", g, os.path.join(wd, "sim1.fq"), "--n", "40000",
                        "--length", "150", "--quality", "novaseqx", "--mates", os.path.join(wd, "sim2.fq"),
                        "--insert", "300,60", "--seed", "6"], check=True, stdout=subprocess.DEVNULL)
        subprocess.run([certa, "index", g, "-o", os.path.join(wd, "genome.cidx"), "-t", "4"], check=True,
                       stderr=subprocess.DEVNULL)
    rng = random.Random(17)

    def records(p):
        lines = open(p).read().split("\n")
        return [lines[i:i + 4] for i in range(0, len(lines) - 3, 4)]

    r1, r2 = records(os.path.join(wd, "sim1.fq")), records(os.path.join(wd, "sim2.fq"))
    for a, b in zip(r1, r2):
        u, n = rng.random(), len(a[1])
        if u < 0.10:                                   # adapter read-through
            cut = rng.randint(20, n - 10)
            a[1] = (a[1][:cut] + ADAPTER_R1 + "A" * n)[:n]
            b[1] = (b[1][:cut] + ADAPTER_R2 + "A" * n)[:n]
        elif u < 0.15:                                 # poly-G tail, no signal
            k = rng.randint(12, 40)
            a[1], a[3] = a[1][:n - k] + "G" * k, a[3][:n - k] + "#" * k
        elif u < 0.22:                                 # poor 3' end
            k = rng.randint(10, 60)
            b[3] = b[3][:n - k] + "".join(rng.choice("#$%&'") for _ in range(k))
        elif u < 0.24:                                 # poor throughout: fails the quality filter
            a[3] = "".join(rng.choice("#$%&'()") for _ in range(n))
    text = lambda rs: "".join("\n".join(x) + "\n" for x in rs)
    with open(os.path.join(wd, "in1.fq"), "w") as f: f.write(text(r1))
    with open(os.path.join(wd, "in2.fq"), "w") as f: f.write(text(r2))
    with gzip.open(os.path.join(wd, "in1.fq.gz"), "wt", compresslevel=4) as f: f.write(text(r1))
    with gzip.open(os.path.join(wd, "in2.fq.gz"), "wt", compresslevel=4) as f: f.write(text(r2))


def run_case(a, case):
    name, paired, gz, fset, mode, (w, t) = case
    d = os.path.join(a.workdir, name.replace("/", "_"))
    shutil.rmtree(d, ignore_errors=True)
    two, one = os.path.join(d, "two"), os.path.join(d, "one")
    os.makedirs(two); os.makedirs(one)
    ext = ".fq.gz" if gz else ".fq"
    inputs = ["-i", os.path.join(a.workdir, "in1" + ext)]
    if paired: inputs += ["-I", os.path.join(a.workdir, "in2" + ext)]
    idx = os.path.join(a.workdir, "genome.cidx")
    fopts = FASTP_SETS[fset][0]
    mopts = CERTA_MODES[mode]

    def outputs(dd):
        if mode == "intervals": return ["-o", dd + "/out.sam", "--stats", dd + "/map.json"]
        o = ["-o", dd + "/out.sam", "-u", dd + "/u1.fq", "--stats", dd + "/map.json"]
        return o + (["-U", dd + "/u2.fq"] if paired else [])

    # two-step: fastp writes the trimmed reads, certa map reads them back
    f2 = [x.replace("{D}", two) for x in fopts]
    trimmed = ["-o", two + "/c1.fq"] + (["-O", two + "/c2.fq"] if paired else [])
    rc = run([a.fastp] + inputs + trimmed + f2 + ["-w", str(w), "-j", two + "/fastp.json", "-h", two + "/fastp.html"],
             two + "/fastp.err")
    if rc == 0:
        reads = [two + "/c1.fq"] + ([two + "/c2.fq"] if paired else [])
        rc = run([a.certa, "map", idx] + reads + mopts + ["-t", str(t)] + outputs(two), two + "/map.err")
    if rc != 0:
        return name, "ERROR", "two-step failed: " + open(two + "/fastp.err").read()[-200:]
    # one pass
    f1 = [x.replace("{D}", one) for x in fopts]
    rc = run([a.certa_qc, idx] + mopts + ["-t", str(t)] + outputs(one) +
             ["--fastp"] + inputs + f1 + ["-w", str(w), "-j", one + "/fastp.json", "-h", one + "/fastp.html"],
             one + "/err")
    if rc != 0:
        return name, "FAIL", "certa-qc failed: " + open(one + "/err").read()[-300:].replace("\n", " ")

    problems = []
    for f in sorted(os.listdir(two)):
        if f in ("c1.fq", "c2.fq", "fastp.html", "fastp.json", "map.json") or f.endswith(".err"):
            continue
        if not os.path.exists(os.path.join(one, f)) or content_md5(os.path.join(two, f)) != content_md5(os.path.join(one, f)):
            problems.append(f)
    for f, volatile in (("fastp.json", FASTP_VOLATILE), ("map.json", {"seconds"})):
        x, y = json.load(open(os.path.join(two, f))), json.load(open(os.path.join(one, f)))
        for k in volatile:
            x.pop(k, None); y.pop(k, None)
        if x != y:
            problems.append(f + ":" + ",".join(sorted(k for k in set(x) | set(y) if x.get(k) != y.get(k)))[:80])
    sam_records = sum(1 for line in open(os.path.join(one, "out.sam")) if not line.startswith("@"))
    if not problems:
        shutil.rmtree(d, ignore_errors=True)
    return name, ("FAIL" if problems else "PASS"), (" ".join(problems) if problems else f"{sam_records} SAM records")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--certa", required=True); ap.add_argument("--certa-qc", dest="certa_qc", required=True)
    ap.add_argument("--fastp", required=True); ap.add_argument("--workdir"); ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--only")
    a = ap.parse_args()
    a.workdir = a.workdir or tempfile.mkdtemp(prefix="qc_equiv_")
    os.makedirs(a.workdir, exist_ok=True)
    make_data(a.workdir, a.certa)
    cases = []
    for paired in (False, True):
        for fset, (_, pe_only) in FASTP_SETS.items():
            if pe_only and not paired: continue
            for mode in CERTA_MODES:
                for w, t in THREADS:
                    if fset in ONE_WORKER_SETS and w > 1: continue
                    gz = fset == "default" and w > 1          # gzip input once per mode
                    name = f"{'PE' if paired else 'SE'}/{fset}/{mode}/w{w}t{t}{'/gz' if gz else ''}"
                    if a.only and a.only not in name: continue
                    cases.append((name, paired, gz, fset, mode, (w, t)))
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        rows = sorted(ex.map(lambda c: run_case(a, c), cases))
    for name, status, detail in rows:
        print(f"{status:5} {name}  {detail}")
    fails = sum(1 for r in rows if r[1] != "PASS")
    print(f"cases={len(rows)} pass={len(rows) - fails} fail={fails}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
