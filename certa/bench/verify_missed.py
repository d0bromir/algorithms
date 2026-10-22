"""For reads a tool leaves unmapped though CERTA proves an end-to-end
alignment with <= 5 edits: fetch the reference under CERTA's alignment
(samtools faidx) and recompute the semi-global edit distance of the whole
read there; it must equal CERTA's XE. Also: does the read hold an N?
  verify_missed.py certa.sam tool.sam ref.fa samtools [max reads]"""
import subprocess
import sys

CIG = __import__('re').compile(r'(\d+)([MIDNSHP=X])')
comp = str.maketrans('ACGTN', 'TGCAN')

cert = {}
with open(sys.argv[1]) as f:
    for line in f:
        if line[0] == '@':
            continue
        t = line.rstrip('\n').split('\t')
        flag = int(t[1])
        if flag & 0x904:
            continue
        tags = {x[:2]: x[5:] for x in t[11:]}
        if tags.get('XT') in ('S0', 'S1') and int(tags.get('XE', 99)) <= 5:
            cert[t[0]] = (t[2], int(t[3]), t[5], t[9], flag & 16, int(tags['XE']))
mapped = set()
with open(sys.argv[2]) as f:
    for line in f:
        if line[0] == '@':
            continue
        t = line.split('\t', 3)
        if not int(t[1]) & 4:
            mapped.add(t[0])
miss = [k for k in cert if k not in mapped][: int(sys.argv[5]) if len(sys.argv) > 5 else 20]


def semiglobal(read, ref):  # read fully aligned, reference ends free
    prev = [0] * (len(ref) + 1)
    for i, c in enumerate(read, 1):
        cur = [i] + [0] * len(ref)
        for j, r in enumerate(ref, 1):
            cur[j] = min(prev[j - 1] + (0 if c == r and c != 'N' else 1), prev[j] + 1, cur[j - 1] + 1)
        prev = cur
    return min(prev)


ok = bad = withn = 0
lens = []
for k in miss:
    chrom, pos, cigar, seq, rev, xe = cert[k]
    ops = CIG.findall(cigar)
    span = sum(int(n) for n, o in ops if o in 'MD=XN')
    cl = int(ops[0][0]) if ops[0][1] == 'S' else 0  # soft clips: the end-to-end alignment reaches past them
    cr = int(ops[-1][0]) if ops[-1][1] == 'S' else 0
    a, b = max(1, pos - cl - 20), pos + span + cr + 20
    out = subprocess.run([sys.argv[4], 'faidx', sys.argv[3], f'{chrom}:{a}-{b}'], capture_output=True, text=True).stdout
    ref = ''.join(out.split('\n')[1:]).upper()
    d = semiglobal(seq, ref)  # SAM SEQ is already on the forward strand
    withn += 'N' in seq
    ok += d == xe
    bad += d != xe
    within = xe <= int(0.04 * len(seq))  # the tools' 4 % of the read length
    lens.append((len(seq), within))
    print(f'{k}\t{chrom}:{pos}\tXE {xe}\trecomputed {d}\tlength {len(seq)}\twithin 4 %: {within}\tN: {"N" in seq}')
print(f'checked {ok + bad}: recomputed distance = XE for {ok}, differs for {bad}; reads with N: {withn}; '
      f'shorter than 150: {sum(1 for L, w in lens if L < 150)}; XE within 4 % of the length: '
      f'{sum(1 for L, w in lens if w)}')
