#!/usr/bin/env python3
"""The paper's data figures, drawn from bench/results/ (PDF, into figures/).

Palette: the dataviz reference instance, light mode (print). Categorical
slots 1-4 in fixed order (validated: adjacent CVD dE >= 9.1, normal-vision
>= 22.9; slots 3-4 sit below 3:1 on white, so every figure also carries a
legend or direct labels, and the values are in the tables).

Usage: python make_figures.py   (needs matplotlib; run make_numbers.py too)
"""
import json
import pathlib
import re

import matplotlib
import matplotlib.ticker

matplotlib.use('pdf')
import matplotlib.pyplot as plt  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RES = HERE.parent / 'bench' / 'results'
OUT = HERE / 'figures'
OUT.mkdir(exist_ok=True)

BLUE, ORANGE, AQUA, YELLOW = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
GRAY = '#a3a29c'        # neutral marks ("not certified", reference bars)
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
MM = 1 / 25.4
COL = 86 * MM           # Bioinformatics single column

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7, 'axes.titlesize': 7.5, 'axes.labelsize': 7, 'xtick.labelsize': 6.5,
    'ytick.labelsize': 6.5, 'legend.fontsize': 6.5, 'axes.edgecolor': INK2,
    'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK,
    'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'pdf.fonttype': 42,
    'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
})


def rows(path):
    out, head = [], None
    for ln in (RES / path).read_text(encoding='utf-8').splitlines():
        if not ln.strip() or ln.startswith('#'):
            continue
        f = ln.split('\t')
        if head is None:
            head = f
            continue
        if len(f) == len(head):
            out.append(dict(zip(head, f)))
    return out


def median_wall(path, tool):
    rs = sorted((r for r in rows(path) if r['tool'] == tool), key=lambda r: float(r['wall_s']))
    return float(rs[len(rs) // 2]['wall_s'])


def grid(ax, axis):
    ax.grid(axis=axis, color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)


def hbar(ax, y, left, width, color, h=0.56):
    """A bar of height h, square at its baseline and rounded at the data end
    when it ends a stack; touching segments are separated by a surface gap."""
    ax.barh(y, width, left=left, height=h, color=color, edgecolor='white', linewidth=0.8)


# ---------------------------------------------------------------- Fig. 2
def fig_errors_vs_speed():
    s = (RES / 'hg002_2026-09/weaknesses_2M.txt').read_text(encoding='utf-8')
    cert = int(re.search(r'^reads (\d+);', s, re.M).group(1))
    miss = {m.group(1): int(m.group(2)) for m in re.finditer(r'^(\S+)\s+(\d+) \(', s, re.M)}
    wall = {
        'bwa-mem2': median_wall('hg002_2026-09/galaxy_timing_v06.tsv', 'bwa-mem2'),
        'minibwa': median_wall('hg002_2026-09/galaxy_timing_v06.tsv', 'minibwa'),
        'strobealign': median_wall('hg002_2026-09/galaxy_timing.tsv', 'strobealign'),
        'minimap2': median_wall('hg002_2026-09/galaxy_timing.tsv', 'minimap2'),
        'bowtie2': median_wall('hg002_2026-09/galaxy_timing.tsv', 'bowtie2'),
    }
    certa = median_wall('hg002_2026-09/galaxy_timing_v06.tsv', 'certa-gpu-pipe')
    names = {'bwa-mem2': 'BWA-MEM2', 'minibwa': 'minibwa', 'strobealign': 'strobealign',
             'minimap2': 'minimap2 -x sr', 'bowtie2': 'Bowtie 2'}
    fig, ax = plt.subplots(figsize=(COL, 54 * MM))
    grid(ax, 'both')
    off = {'bwa-mem2': (5, 6), 'minibwa': (5, 7), 'strobealign': (6, -2),
           'minimap2': (6, -2), 'bowtie2': (-4, 7)}
    for t, w in wall.items():
        y = miss[t] / cert * 1e6
        ax.scatter(w, y, s=26, color=INK2, edgecolor='white', linewidth=1.0, zorder=3)
        ha = 'right' if off[t][0] < 0 else 'left'
        ax.annotate(names[t], (w, y), xytext=off[t], textcoords='offset points',
                    ha=ha, va='center', fontsize=6.5, color=INK)
    ax.scatter(certa, 0, s=30, color=BLUE, edgecolor='white', linewidth=1.0, zorder=4, clip_on=False)
    ax.annotate('CERTA: 0 by construction', (certa, 0), xytext=(7, 0),
                textcoords='offset points', ha='left', va='center', fontsize=6.5, color=INK)
    ax.set_xscale('log')
    ax.set_xticks([20, 30, 50, 100, 200])
    ax.set_xticklabels(['20', '30', '50', '100', '200'])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlim(18, 260)
    ax.set_ylim(-150, 2900)
    ax.set_xlabel('Wall time, 20 M reads, 64 threads (s, log scale)')
    ax.set_ylabel('Provably misplaced reads\nper million certified reads')
    fig.tight_layout(pad=0.3)
    fig.savefig(OUT / 'errors_vs_speed.pdf', metadata={'CreationDate': None})
    plt.close(fig)


# ---------------------------------------------------------------- Fig. 3
def fig_certification():
    sets = [
        ('HG002 single-end', 'hg002_2026-09/full_444M_v06_run4.stats.json'),
        ('HG002 paired', 'hg002_2026-10/certa_pe_full.stats.json'),
        ('HG001 paired', 'hg002_2026-10/certa_pe_hg001.stats.json'),
        ('HG005 paired', 'hg002_2026-10/certa_pe_hg005.stats.json'),
    ]
    segs = [('Exact (S0)', BLUE), ('1 to R edits (S1)', ORANGE), ('Certified repeat (SR)', AQUA),
            ('Pair certificate (PR)', YELLOW), ('Not certified', GRAY)]
    fig, ax = plt.subplots(figsize=(COL, 40 * MM))
    for i, (label, path) in enumerate(sets):
        d = json.loads((RES / path).read_text(encoding='utf-8'))
        n = d['reads']
        pr = d['certified'] - sum(d.get(k, 0) for k in ('s0', 's1', 'sr', 'sl', 's2'))
        vals = [d['s0'], d['s1'], d['sr'], pr, n - d['certified']]
        left = 0.0
        y = len(sets) - 1 - i
        for (name, color), v in zip(segs, vals):
            w = 100 * v / n
            hbar(ax, y, left, w, color)
            left += w
        ax.text(100.8, y, f"{100 * d['certified'] / n:.2f} %", va='center', ha='left',
                fontsize=6.5, color=INK)
    ax.set_yticks(range(len(sets)))
    ax.set_yticklabels([s[0] for s in reversed(sets)])
    ax.tick_params(axis='y', length=0)
    ax.spines['left'].set_visible(False)
    ax.set_xlim(0, 100)
    ax.set_xlabel('Share of reads (%)')
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, c in segs]
    ax.legend(handles, [s for s, _ in segs], ncol=3, loc='lower left', bbox_to_anchor=(-0.32, 1.0),
              handlelength=1.0, handleheight=0.8, columnspacing=1.0, borderaxespad=0.2)
    fig.tight_layout(pad=0.3)
    fig.subplots_adjust(right=0.86)
    fig.savefig(OUT / 'certification.pdf', metadata={'CreationDate': None})
    plt.close(fig)


# ---------------------------------------------------------------- Fig. 4
def fig_hifi_cost():
    s = (RES / 'aim6_2026-10/certify_cost.txt').read_text(encoding='utf-8')
    cpu = [int(x.replace(',', '')) for x in
           re.search(r'CPU \(user\)\s+([\d,]+) s\s+([\d,]+) s\s+([\d,]+) s', s).groups()]
    mm2 = int(re.search(r"minimap2's mapping and alignment \(([\d,]+) s CPU\)", s).group(1).replace(',', ''))
    nreads = int(re.search(r'the full sample \(([\d,]+) reads\)', s).group(1).replace(',', ''))
    bars = [('Scalar DP', cpu[0], BLUE), ('Vectorized DP row', cpu[1], BLUE),
            ('+ batched lookups,\n   DP pruned to live cells', cpu[2], BLUE),
            ('minimap2 -x map-hifi\n(map and align)', mm2, GRAY)]
    fig, ax = plt.subplots(figsize=(COL, 40 * MM))
    grid(ax, 'x')
    for i, (label, v, c) in enumerate(bars):
        y = len(bars) - 1 - i
        ax.barh(y, v, height=0.56, color=c)
        ax.text(v + 150, y, f'{v:,} s', va='center', ha='left', fontsize=6.5, color=INK)
    ax.set_yticks(range(len(bars)))
    ax.set_yticklabels([b[0] for b in reversed(bars)])
    ax.tick_params(axis='y', length=0)
    ax.spines['left'].set_visible(False)
    ax.set_xlim(0, 14800)
    ax.set_xticks([0, 4000, 8000, 12000])
    ax.set_xticklabels(['0', '4,000', '8,000', '12,000'])
    ax.set_xlabel(f'CPU time for {nreads:,} HiFi reads (s)')
    fig.tight_layout(pad=0.3)
    fig.savefig(OUT / 'hifi_cost.pdf', metadata={'CreationDate': None})
    plt.close(fig)


fig_errors_vs_speed()
fig_certification()
fig_hifi_cost()
print('figures:', ', '.join(sorted(p.name for p in OUT.glob('*.pdf'))))
