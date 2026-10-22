"""Of the reads a tool leaves unmapped though CERTA proves an end-to-end
alignment: how many have a proven edit count within the tool's error rate
(4 % of the read length, floored), i.e. the misses the tool's own promise
covers.  within_rate.py certa.sam tool.sam"""
import sys

cert = {}
with open(sys.argv[1]) as f:
    for line in f:
        if line[0] == '@':
            continue
        t = line.rstrip('\n').split('\t')
        flag = int(t[1])
        if flag & 0x904:
            continue
        key = t[0] + ('/1' if flag & 64 else '/2' if flag & 128 else '')
        tags = {x[:2]: x[5:] for x in t[11:]}
        if tags.get('XT') in ('S0', 'S1', 'SR'):
            cert[key] = (int(tags['XE']), len(t[9]))
mapped = set()
with open(sys.argv[2]) as f:
    for line in f:
        if line[0] == '@':
            continue
        t = line.split('\t', 3)
        flag = int(t[1])
        if not flag & 4:
            mapped.add(t[0] + ('/1' if flag & 64 else '/2' if flag & 128 else ''))
within = beyond = 0
for k, (xe, L) in cert.items():
    if k in mapped:
        continue
    if xe <= int(0.04 * L):
        within += 1
    else:
        beyond += 1
n_within = sum(1 for xe, L in cert.values() if xe <= int(0.04 * L))
print(f'end-to-end certified {len(cert)}, of which with a proven edit count within 4 % of the length {n_within}; '
      f'unmapped by the tool: {within} within its error rate, {beyond} beyond it')
