"""speed.csv for the dashboard from the final protocol (one session per host, 64 cores)."""
import csv
import pathlib
import sys

RES = pathlib.Path(sys.argv[1])
OUT = pathlib.Path(sys.argv[2])


def rows(p):
    out, head = [], None
    for ln in (RES / p).read_text(encoding='utf-8').splitlines():
        if not ln.strip() or ln.startswith('#'):
            continue
        f = ln.split('\t')
        if head is None:
            head = f
        elif len(f) == len(head):
            out.append(dict(zip(head, f)))
    return out


NAMES = {'bwa-mem2': 'BWA-MEM2', 'minibwa': 'minibwa', 'strobealign': 'strobealign', 'minimap2': 'minimap2 -x sr',
         'bowtie2': 'Bowtie 2', 'certa-cpu': 'CERTA CPU + minibwa', 'certa-cpu-iv': 'CERTA CPU interval mode',
         'certa-gpu': 'CERTA GPU + minibwa', 'certa-gpu-iv': 'CERTA GPU interval mode',
         'certa-cpu-tier': 'CERTA CPU certified tier only', 'certa-gpu-tier': 'CERTA GPU certified tier only'}
HOSTS = {'galaxy': ('hg002_2026-10/final_galaxy.tsv', 'ARM64, 64 of 128 cores'),
         'a2': ('hg002_2026-10/final_a2.tsv', 'x86-64, 64 cores, memory interleaved over 4 NUMA nodes')}
out = []
for host, (path, hw) in HOSTS.items():
    by = {}
    for r in rows(path):
        by.setdefault((r['mode'], r['pipeline'], r['config']), []).append(r)
    for (mode, pipe, cfg), rs in by.items():
        rs = sorted(rs, key=lambda r: float(r['wall_s']))
        r = rs[len(rs) // 2]
        gpu = pipe.startswith('certa-gpu')
        out.append(dict(host=host, mode=mode, pipeline=NAMES[pipe], config=cfg,
                        hardware=hw + (' + A100' if gpu else ''), reads=20001073,
                        wall_s=float(r['wall_s']), cpu_s=round(float(r['user_s']) + float(r['sys_s']), 2),
                        runs=len(rs), complete='no' if pipe.endswith('tier') else 'yes'))
with open(OUT, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=list(out[0]))
    w.writeheader()
    w.writerows(out)
for r in out:
    print(r['host'], r['mode'], r['pipeline'], r['config'], r['wall_s'], r['cpu_s'], r['runs'])
