"""The dashboard's slow-mapper datasets from bench/results/hg002_2026-10:
slow_speed.csv (timings, medians), slow_reliability.csv (weaknesses.py
counts) and slow_hifi.csv (HiFi certificates on each aligner's alignments).

  python make_slow_csv.py RESULTS_DIR OUT_DIR
"""
import csv
import pathlib
import re
import statistics
import sys

RES = pathlib.Path(sys.argv[1]) / 'hg002_2026-10'
OUT = pathlib.Path(sys.argv[2])  # reads: reads, or pairs (paired)


def tsv(name):
    lines = [ln for ln in (RES / name).read_text(encoding='utf-8').splitlines() if ln and not ln.startswith('#')]
    head = lines[0].split('\t')
    return [dict(zip(head, ln.split('\t'))) for ln in lines[1:]]


def write(name, rows):
    with open(OUT / name, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(name, len(rows), 'rows')


TOOLS = {'certa-cpu': 'CERTA', 'certa-gpu': 'CERTA + A100', 'yara': 'Yara 1.0.5', 'razers3': 'RazerS 3 3.5.12',
         'minimap2': 'minimap2 2.31', 'lr_ceiling': 'lr_ceiling', 'winnowmap': 'Winnowmap 2.03'}
SETTINGS = {
    'certa': 'certificate radius <= 5 edits',
    'yara': 'every alignment with <= 4 % errors (6 edits at 150 bp)',
    'razers3': '100 % of matches with >= 96 % identity (6 edits at 150 bp)',
    'minimap2': 'HiFi preset, with CIGARs',
    'lr_ceiling': 'proves the optimum of each read',
    'winnowmap': 'HiFi preset, repetitive 15-mers from meryl; indexes the FASTA in each run',
}
cert = {(r['data'], r['config']): round(float(r['certified_pct']) / 100, 4) for r in tsv('slow_certified.tsv')}
index_s = float(re.search(r'winnowmap: ([0-9.]+) s of its run',
                          (RES / 'slow_hifi_and_logs.txt').read_text(encoding='utf-8')).group(1))
groups = {}
for r in tsv('slow_galaxy.tsv'):
    groups.setdefault((r['data'], r['tool'], r['config']), []).append(r)
speed = []
for (data, tool, cfg), rs in groups.items():
    med = lambda k: statistics.median(float(x[k]) for x in rs)
    config = cfg.replace(' --intervals', '').replace('default', '(default)') if tool.startswith('certa') else cfg
    reads = {'single-end': 2000000, 'paired': 2000000, 'single-end-200000': 200000, 'paired-200000': 200000,
             'hifi': 242534}[data]
    role = 'mapper'
    if tool == 'lr_ceiling':
        role = 'CERTA HiFi pipeline' if 'minimap2' in cfg else 'check of Winnowmap'
    elif tool == 'minimap2':
        role = 'CERTA HiFi pipeline'
    speed.append(dict(
        data=data, reads=reads, tool=TOOLS[tool], config=config,
        settings=SETTINGS['certa' if tool.startswith('certa') else tool], role=role,
        wall_s=round(med('wall_s'), 2), cpu_s=round(med('user_s') + med('sys_s'), 1),
        rss_gb=round(max(int(x['max_rss_kb']) for x in rs) / 1e6, 1), runs=len(rs),
        status='ok' if rs[0]['exit'] == '0' else f"aborted (exit {rs[0]['exit']})",
        certified='' if not tool.startswith('certa') else cert[(data, cfg.replace(' --intervals', ''))],
        index_s=index_s if tool == 'winnowmap' else ''))
write('slow_speed.csv', speed)

# weaknesses.py output: a "==" title per comparison, then its two tables.
rel, title = [], None
NAMES = {'Yara': 'Yara 1.0.5', 'RazerS3:all': 'RazerS 3 (best record)', 'RazerS3-primary': 'RazerS 3 (flagged primary)'}
text = (RES / 'slow_reliability.txt').read_text(encoding='utf-8')
for block in text.split('\n== ')[1:]:
    title = block.splitlines()[0]
    data = 'single-end' if title.startswith('single-end, 2,000,000') else 'paired' if title.startswith('paired') else 'single-end-200000'
    ties = int(re.search(r'proven exact multi-copy ties (\d+)', block).group(1))
    m = re.search(r'end-to-end \(S0 S1 SR\) (\d+), local \(SL PR\) (\d+)', block)
    n_e2e, n_local = int(m.group(1)), int(m.group(2))
    main = {ln.split()[0]: ln for ln in block.splitlines() if ln.split() and ln.split()[0] in NAMES and 'of ties' in ln}
    tier = {ln.split()[0]: ln for ln in block.splitlines() if ln.split() and ln.split()[0] in NAMES and 'of ties' not in ln}
    for key, name in NAMES.items():
        if key not in main:
            continue
        n = [int(x) for x in re.findall(r'(?<![.(\d])(\d+)(?![.\d%])', tier[key].split(None, 1)[1])]
        over = int(re.findall(r'(\d+) \([0-9.]+% of ties\)', main[key])[0])
        rel.append(dict(data=data, tool=name, end_to_end_certified=n_e2e, e2e_missed_better=n[0],
                        e2e_more_edits=n[1], e2e_unmapped=n[2], local_certified=n_local, local_missed_better=n[3],
                        local_unmapped=n[4], proven_ties=ties, overconfident_ties=over))
write('slow_reliability.csv', rel)

hifi = []
logs = (RES / 'slow_hifi_and_logs.txt').read_text(encoding='utf-8')
for name, log, tag in (('minimap2 2.31', 'lr_mm2.log', 'minimap2'), ('Winnowmap 2.03', 'lr_wm.log', 'winnowmap')):
    sec = logs.split(f'== {log}')[1].split('==')[0]
    proven, above = re.search(r'OPT proven for (\d+) reads .*above the best known alignment for (\d+)', sec).groups()
    aligned = re.search(rf'^{tag}: \d+ primary alignments, (\d+) reads with one', logs, re.M).group(1)
    hifi.append(dict(aligner=name, reads=242534, reads_aligned=int(aligned), proven=int(proven),
                     proven_above_aligner=int(above)))
write('slow_hifi.csv', hifi)
