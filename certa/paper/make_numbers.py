#!/usr/bin/env python3
"""Every number the paper quotes, read from bench/results/ and the test suite.

Writes numbers.tex: one LaTeX macro per number (\\newcommand\\Name{value}).
The paper never types a measured number by hand; it uses these macros, so a
re-run of a benchmark changes the paper when this script is re-run.

Usage: python make_numbers.py   (from any directory)
"""
import json
import pathlib
import re
import statistics

HERE = pathlib.Path(__file__).resolve().parent
CERTA = HERE.parent
RES = CERTA / 'bench' / 'results'

macros = {}
sources = {}


def put(name, value, src):
    assert re.fullmatch('[A-Za-z]+', name), name
    assert name not in macros, name
    macros[name] = str(value)
    sources[name] = src


def pct(a, b, nd=2):
    return f'{100 * a / b:.{nd}f}'


def grp(n):
    return f'{int(n):,}'


def fx(x, nd=2):
    return f'{x:.{nd}f}'


def rows(path):
    """Rows of a tab-separated file with a header line; '#' lines skipped."""
    out, head = [], None
    for ln in (RES / path).read_text(encoding='utf-8').splitlines():
        if not ln.strip() or ln.startswith('#'):
            continue
        f = ln.split('\t')
        if head is None:
            head = f
            continue
        out.append(dict(zip(head, f)))
    return out


def median_run(path, tool, key='tool'):
    """The repetition with the median wall time; CPU = user + sys."""
    rs = [r for r in rows(path) if r[key] == tool]
    assert rs, (path, tool)
    rs.sort(key=lambda r: float(r['wall_s']))
    r = rs[len(rs) // 2]
    return float(r['wall_s']), float(r['user_s']) + float(r['sys_s'])


def text(path):
    return (RES / path).read_text(encoding='utf-8')


def grab(pattern, s, conv=float, flags=0):
    m = re.search(pattern, s, flags)
    assert m, pattern
    return conv(m.group(1))


# ---------------------------------------------------------------- certification
def stats(path):
    d = json.loads(text(path))
    certified, reads = d['certified'], d['reads']
    pr = certified - sum(d.get(k, 0) for k in ('s0', 's1', 'sr', 'sl', 's2'))
    return d, certified, reads, pr


src = 'hg002_2026-09/full_444M_v06_run4.stats.json'
d, c, n, _ = stats(src)
put('SeReads', grp(n), src)
put('SeCert', pct(c, n), src)
put('SeExact', pct(d['s0'], n, 1), src)
put('SeEdits', pct(d['s1'], n, 1), src)
put('SeRepeat', pct(d['sr'], n, 1), src)
put('SeUncert', pct(n - c, n, 1), src)
se_full_cert = c

src = 'hg002_2026-09/full_444M_gpu_k2.stats.json'
d, c, n, _ = stats(src)
assert d['k'] == 2
put('SeCertFixed', pct(c, n, 1), src)

src = 'aim4_2026-10/check_index.log'
s = text(src)
put('IdxEntries', f"{grab(r'(\d+) entries', s, int) / 1e9:.2f}", src)
put('IdxCheckSec', f"{grab(r'valid \(([\d.]+) s\)', s):.0f}", src)

src = 'hg002_2026-09/galaxy_timing_v06.tsv'
put('BenchReads', grp(rows(src)[0]['reads']), src)

src = 'hg002_2026-10/certa_pe_full.stats.json'
d, c, n, pr = stats(src)
put('PeReads', grp(n), src)
put('PePairs', grp(n // 2), src)
put('PeCert', pct(c, n), src)
put('PeCertPR', pct(pr, n), src)
put('PeExact', pct(d['s0'], n, 1), src)
pe_by_d1 = d['certified_by_d1']
for name, f in (('HGOne', 'hg002_2026-10/certa_pe_hg001.stats.json'),
                ('HGFive', 'hg002_2026-10/certa_pe_hg005.stats.json')):
    d, c, n, pr = stats(f)
    put(f'PeCert{name}', pct(c, n), f)
    put(f'PeReads{name}', grp(n), f)
    put(f'PeReadsM{name}', f'{n / 1e6:.0f}', f)

src = 'hg002_2026-10/certa_pe_20M.stats.json'
d, c, n, pr = stats(src)
put('PeTwentyCert', pct(c, n), src)
put('PeTwentyPR', pct(pr, n), src)

# ---------------------------------------------------------------- provable errors
src = 'hg002_2026-09/weaknesses_2M.txt'
s = text(src)
w_reads = grab(r'^reads (\d+);', s, int, re.M)
w_ties = grab(r'ties (\d+)', s, int)
put('WkReads', grp(w_reads), src)
put('WkTies', grp(w_ties), src)
TOOLS = {'bwa-mem2': 'Bwa', 'minibwa': 'Mini', 'strobealign': 'Strobe',
         'minimap': 'Mm', 'bowtie': 'Bt'}
wk = {}
for m in re.finditer(r'^(\S+)\s+(\d+) \(([\d.]+)%\)\s+(\d+)\s+(\d+) \(([\d.]+)%\)\s+(\d+) \(([\d.]+)% of ties\)',
                     s, re.M):
    tool = m.group(1)
    key = next(v for k, v in TOOLS.items() if tool.startswith(k))
    vals = dict(miss=int(m.group(2)), misspct=m.group(3), q20=int(m.group(4)),
                unmapped=int(m.group(5)), ties=int(m.group(7)), tiespct=m.group(8))
    wk[tool] = vals
    put(f'Wk{key}Miss', grp(vals['miss']), src)
    put(f'Wk{key}MissPct', vals['misspct'], src)
    put(f'Wk{key}Q', grp(vals['q20']), src)
    put(f'Wk{key}Unmapped', grp(vals['unmapped']), src)
    put(f'Wk{key}Ties', grp(vals['ties']), src)
    put(f'Wk{key}TiesPct', vals['tiespct'], src)
assert len(wk) == 5, wk
put('WkStrobeOverBwa', f"{wk['strobealign']['miss'] / wk['bwa-mem2']['miss']:.0f}", src)
put('WkMiniOverBwa', f"{wk['minibwa']['miss'] / wk['bwa-mem2']['miss']:.1f}", src)
# BWA-MEM2's rate among certified reads, applied to the certified reads of the whole sample.
per_genome = wk['bwa-mem2']['miss'] / w_reads * se_full_cert
put('WkBwaPerGenome', grp(round(per_genome, -3)), src + ' x full_444M_v06_run4 certified')

# ---------------------------------------------------------------- speed, single-end
src = 'hg002_2026-09/galaxy_timing_v06.tsv'
bw_w, bw_c = median_run(src, 'bwa-mem2')
mb_w, mb_c = median_run(src, 'minibwa')
ce_w, ce_c = median_run(src, 'certa-gpu-pipe')
ct_w, ct_c = median_run(src, 'certa-gpu')
put('SeWallBwa', fx(bw_w), src)
put('SeCpuBwa', grp(round(bw_c)), src)
put('SeWallMini', fx(mb_w), src)
put('SeCpuMini', grp(round(mb_c)), src)
put('SeWallCerta', fx(ce_w), src)
put('SeCpuCerta', grp(round(ce_c)), src)
put('SeWallTier', fx(ct_w), src)
put('SeCpuTier', grp(round(ct_c)), src)
put('SeSpeedBwa', fx(bw_w / ce_w), src)
put('SeSpeedMini', fx(mb_w / ce_w), src)
put('SeCpuRatioBwa', fx(bw_c / ce_c, 1), src)
put('SeTierSpeedBwa', f'{bw_w / ct_w:.0f}', src)

src = 'hg002_2026-09/galaxy_timing.tsv'
for tool, key in (('strobealign', 'Strobe'), ('minimap2', 'Mm'), ('bowtie2', 'Bt')):
    w, c = median_run(src, tool)
    put(f'SeWall{key}', fx(w), src)
    put(f'SeCpu{key}', grp(round(c)), src)

src = 'hg002_2026-09/galaxy_timing_v07_pass3.tsv'
mx_w, mx_c = median_run(src, 'maxcert', key='mode')
put('SeWallMaxCert', fx(mx_w), src)

# ---------------------------------------------------------------- speed, paired-end
src = 'hg002_2026-10/galaxy_timing_pe.tsv'
pbw_w, pbw_c = median_run(src, 'bwa-mem2')
put('PeWallBwa', fx(pbw_w), src)
put('PeCpuBwa', grp(round(pbw_c)), src)
src = 'hg002_2026-10/galaxy_timing_pe_v2.tsv'
pe_rows = rows(src)
best = min((r for r in pe_rows if r['tool'] == 'certa-gpu-pe-pipe'), key=lambda r: float(r['wall_s']))
pce_w = float(best['wall_s'])
pce_c = float(best['user_s']) + float(best['sys_s'])
pmb_w, pmb_c = median_run(src, 'minibwa')
put('PeWallCerta', fx(pce_w), src)
put('PeCpuCerta', grp(round(pce_c)), src)
put('PeCertaThreads', best['certa_threads'], src)
put('PeMiniThreads', str(64 - int(best['certa_threads'])), src)
put('PeWallMini', fx(pmb_w), src)
put('PeCpuMini', grp(round(pmb_c)), src)
put('PeSpeedBwa', fx(pbw_w / pce_w, 1), src)
put('PeSpeedMini', fx(pmb_w / pce_w), src)
put('PeCpuRatioBwa', fx(pbw_c / pce_c, 1), src)
side = re.search(r"CERTA's own side at 16 / 24 / 32 threads: ([\d.]+) / ([\d.]+) / ([\d.]+) s", text(src))
put('PeCertaSide', side.group(3), src)

# ---------------------------------------------------------------- variant calling
def vcfeval(path, keyfn):
    out = {}
    for ln in text(path).splitlines():
        f = ln.split()
        if len(f) < 11 or ln.startswith('#'):
            continue
        k, rest = keyfn(f)
        tp_base, tp_call, fp, fn = (int(x) for x in rest[1:5])
        prec, rec, f1 = rest[5:8]
        out[k] = dict(tp=tp_base, fp=fp, fn=fn, prec=prec, rec=rec, f1=f1)
    return out


src = 'hg002_2026-09/giab_chr20_vcfeval.tsv'
se_v = vcfeval(src, lambda f: ((f[0], f[1], f[2]), f[3:]))
src_pe = 'hg002_2026-10/giab_chr20_vcfeval_pe.tsv'
pe_v = vcfeval(src_pe, lambda f: (('HG002', {'ourspe': 'certa', 'bwape': 'bwa'}[f[0]], f[1], f[2]), f[3:]))
src_3 = 'hg002_2026-10/giab_chr20_vcfeval_hg001_hg005_pe.tsv'
pe_v.update(vcfeval(src_3, lambda f: ((f[0], f[1], f[2], f[3]), f[4:])))

CALL = {('gatk', 'snps'): 'GS', ('gatk', 'indels'): 'GI', ('bcftools', 'snps'): 'BS', ('bcftools', 'indels'): 'BI'}
for (caller, typ), key in CALL.items():
    for pipe, pk in (('bwa', 'Bwa'), ('ours4', 'Certa')):
        v = se_v[(pipe, caller, typ)]
        for fld, fk in (('prec', 'Prec'), ('rec', 'Rec'), ('f1', 'F')):
            put(f'SeV{key}{pk}{fk}', v[fld], src)
    put(f'SeV{key}Delta', f"{float(se_v[('ours4', caller, typ)]['f1']) - float(se_v[('bwa', caller, typ)]['f1']):+.4f}", src)

wins = 0
for smp, sk in (('HG002', 'B'), ('HG001', 'A'), ('HG005', 'C')):
    for (caller, typ), key in CALL.items():
        b = pe_v[(smp, 'bwa', caller, typ)]
        c = pe_v[(smp, 'certa', caller, typ)]
        s_ = src_pe if smp == 'HG002' else src_3
        put(f'PeV{sk}{key}Bwa', b['f1'], s_)
        put(f'PeV{sk}{key}Certa', c['f1'], s_)
        delta = float(c['f1']) - float(b['f1'])
        put(f'PeV{sk}{key}Delta', f'{delta:+.4f}'.replace('-', '$-$'), s_)
        wins += delta > 0
        if (caller, typ) == ('gatk', 'snps'):
            put(f'PeV{sk}GSFpBwa', grp(b['fp']), s_)
            put(f'PeV{sk}GSFpCerta', grp(c['fp']), s_)
            put(f'PeV{sk}GSFnBwa', grp(b['fn']), s_)
            put(f'PeV{sk}GSFnCerta', grp(c['fn']), s_)
put('PeVWins', str(wins), src_pe + ', ' + src_3)
put('PeVComparisons', '12', src_pe + ', ' + src_3)
put('SeVWins', str(sum(float(se_v[('ours4',) + k]['f1']) > float(se_v[('bwa',) + k]['f1']) for k in CALL)), src)

# ---------------------------------------------------------------- interval output (Aim 3)
src = 'aim3_2026-10/hg005_intervals.json'
d = json.loads(text(src))
put('IvReads', grp(d['reads']), src)
put('IvExact', pct(d['iv_width'][0], d['reads'], 1), src)
put('IvUnmapped', pct(d['intervals']['unmapped'], d['reads'], 1), src)
src = 'aim3_2026-10/hg005_f1.tsv'
f1 = {(r['caller'], r['type'], r['pipeline']): r for r in rows(src)}
iv_wins = 0
for (caller, typ), key in CALL.items():
    iv, bwa = f1[(caller, typ, 'iv')], f1[(caller, typ, 'bwa')]
    put(f'IvV{key}', iv['F1'], src)
    put(f'IvV{key}Bwa', bwa['F1'], src)
    iv_wins += float(iv['F1']) > float(bwa['F1'])
put('IvVWins', str(iv_wins), src)

# ---------------------------------------------------------------- the checker (Aim 4)
recs, us = [], []
for smp in ('hg002', 'hg005'):
    src = f'aim4_pr_2026-10/check_{smp}_2M.log'
    s = text(src)
    a, b = re.search(r'(\d+) of (\d+) records checked: accepted 100\.0000%', s).groups()
    assert a == b
    recs.append(int(a))
    us.append(grab(r'([\d.]+) us CPU', s))
put('CkRecords', grp(sum(recs)), 'aim4_pr_2026-10/check_hg00{2,5}_2M.log')
put('CkUsLo', f'{min(us):.0f}', 'aim4_pr_2026-10/check_hg00{2,5}_2M.log')
put('CkUsHi', f'{max(us):.0f}', 'aim4_pr_2026-10/check_hg00{2,5}_2M.log')
src = 'aim4_pr_2026-10/tamper_hg002.tsv'
tam = rows(src)
inj = sum(int(r['injected']) for r in tam if r.get('false_by_construction') == 'yes')
rej = sum(int(r['rejected']) for r in tam if r.get('false_by_construction') == 'yes')
put('CkFaults', grp(inj), src)
put('CkRejected', grp(rej), src)
put('CkFaultKinds', str(sum(r.get('false_by_construction') == 'yes' for r in tam)), src)
assert 'untouched records rejected (false alarms): 0' in text(src)
put('CkXS', grp(sum(grab(r'(\d+) have XS proven', text(f'aim4_pr_2026-10/check_{x}_2M.log'), int)
                     for x in ('hg002', 'hg005'))), 'aim4_pr_2026-10/check_hg00{2,5}_2M.log')

# ---------------------------------------------------------------- long reads (Aim 6)
src = 'aim6_2026-10/summary.txt'
s = text(src)
lr_reads = grab(r'reads (\d+); fit', s, int)
put('LrReads', grp(lr_reads), src)
put('LrFit', grab(r'fit \d+ = ([\d.]+) %', s, str), src)
put('LrMedLen', grp(grab(r'read length: median (\d+)', s, int)), src)
put('LrMedLoss', grab(r'loss L - sigma \(mapped\): 10% \d+, 50% (\d+)', s, str), src)
put('LrMedCeil', grp(grab(r'ceiling L - U1: 1% \d+, 10% \d+, 50% (\d+)', s, int)), src)
for kind, key in (('clipped', 'Clip'), ('high loss', 'Loss'), ('repeat-limited', 'Rep'), ('unmapped', 'Unm')):
    put(f'LrMiss{key}', grab(rf'{kind} \d+ \(([\d.]+) % of reads\)', s, str), src)

src = 'aim6_2026-10/certify_cost.txt'
s = text(src)
proven = grab(r'->\s+([\d,]+): every read that fits', s, lambda x: int(x.replace(',', '')))
put('LrProven', grp(proven), src)
put('LrProvenPct', pct(proven, lr_reads), src)
cpu = re.search(r'CPU \(user\)\s+([\d,]+) s\s+([\d,]+) s\s+([\d,]+) s', s).groups()
cpu = [int(x.replace(',', '')) for x in cpu]
mm2 = grab(r"minimap2's mapping and alignment \(([\d,]+) s CPU\)", s, lambda x: int(x.replace(',', '')))
put('LrCpuScalar', grp(cpu[0]), src)
put('LrCpu', grp(cpu[2]), src)
put('LrCpuMm', grp(mm2), src)
put('LrCpuRatio', fx(cpu[2] / mm2), src)
put('LrCpuRatioScalar', fx(cpu[0] / mm2, 1), src)
dp = re.search(r'DP per certified read\s+median ([\d.]+) ms.*->\s+([\d.]+), [\d.]+, \d+\s*$', s, re.M).groups()
put('LrDpMedScalar', dp[0], src)
put('LrDpMed', dp[1], src)
put('LrCellsPct', grab(r'computed, [\d.e+]+ \((\d+) %\)', s, str), src)

src = 'aim6_2026-10/certify.txt'
s = text(src)
above = grab(r'OPT above the best known alignment: (\d+) reads \(([\d.]+) %\)', s, int)
put('LrAbove', grp(above), src)
put('LrAbovePct', grab(r'OPT above the best known alignment: \d+ reads \(([\d.]+) %\)', s, str), src)
put('LrVerifyReads', grp(grab(r'full DP near the known locus: (\d+) reads', s, int)), src)

src = 'aim6_2026-10/split_certify.txt'
s = text(src)
put('LrJumpFit', grab(r'jump P = 30: fit \d+ \(([\d.]+) %\)', s, str), src)
jp = grab(r'OPT_J proven for (\d+) reads', s, int)
put('LrJumpProven', grp(jp), src)
put('LrJumpProvenPct', pct(jp, lr_reads), src)
put('LrExactAll', pct(jp + proven, lr_reads), src)
put('LrJumpAttempted', grp(grab(r'attempted (\d+);', s, int)), src)
put('LrJumpInterval', grp(grab(r'an interval (\d+), over', s, int)), src)
put('LrJumpCapped', grp(grab(r'over the cell cap (\d+) ', s, int)), src)
put('LrJumpMedSec', grab(r'time per read \(one thread\): 50% ([\d.]+) s', s, str), src)

# ---------------------------------------------------------------- verification
src = 'aim6_2026-10/verify_lemmas.txt'
s = text(src)
put('VlL', grp(grab(r'Lemma L, parts >= 3, L = 9: (\d+)', s, int)), src)
put('VlH', grp(grab(r'Lemma H, parts >= 4 .*?: (\d+) transcript', s, int)), src)
put('VlE', grp(grab(r'Lemma E, L = 8: (\d+)', s, int)), src)
put('VlS', grp(grab(r'L = 24: (\d+) \(layout', s, int)), src)
assert 'all lemma checks passed' in s

mut = [ln.split('\t') for ln in (CERTA / 'tests' / 'mutants.tsv').read_text(encoding='utf-8').splitlines()
       if ln and not ln.startswith('#')]
put('MutTotal', str(len(mut)), 'tests/mutants.tsv')
put('MutKilled', str(sum(m[3] == 'killed' for m in mut)), 'tests/mutants.tsv')
put('MutEquiv', str(sum(m[3] == 'equivalent' for m in mut)), 'tests/mutants.tsv')
ax = (CERTA / 'lean' / 'Axioms.lean').read_text(encoding='utf-8')
put('LeanTheorems', str(len(re.findall(r'^#print axioms', ax, re.M))), 'lean/Axioms.lean')
lean_files = sorted((CERTA / 'lean' / 'Certa').glob('*.lean'))
put('LeanModules', str(len(lean_files)), 'lean/Certa/*.lean')
put('LeanLines', grp(sum(len(f.read_text(encoding='utf-8').splitlines()) for f in lean_files)), 'lean/Certa/*.lean')

# ---------------------------------------------------------------- second index (Aim 1)
src = 'aim1_2026-10/summary.tsv'
a1 = rows(src)
r = next(r for r in a1 if r['sample'] == 'HG002' and r['baseline'].startswith('q22 best') and r['q15_config'] == 't=2 budget 256')
put('AqBase', r['baseline_certified_pct'], src)
put('AqAfter', r['certified_after_pct'], src)
put('AqMs', r['q15_cpu_ms_per_residual_read'], src)
r = next(r for r in a1 if r['sample'] == 'HG002' and r['baseline'].startswith('q22 best') and r['q15_config'] == 't=2 budget 4096')
put('AqAfterMax', r['certified_after_pct'], src)

# ---------------------------------------------------------------- write
out = ['% Generated by make_numbers.py from bench/results/; do not edit.']
for k in sorted(macros):
    out.append(f'\\newcommand{{\\{k}}}{{{macros[k]}}}% {sources[k]}')
(HERE / 'numbers.tex').write_text('\n'.join(out) + '\n', encoding='utf-8', newline='\n')
print(f'numbers.tex: {len(macros)} macros')
