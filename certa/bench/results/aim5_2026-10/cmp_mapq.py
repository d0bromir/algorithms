import collections
import sys


def load(path):
    out = {}
    with open(path) as f:
        for line in f:
            if line[0] == "@":
                continue
            t = line.split("\t")
            tier = "other"
            for x in t[11:]:
                if x.startswith("XT:Z:"):
                    tier = x[5:].strip()
            out[(t[0], t[1])] = (int(t[4]), tier, t[2], t[3], t[5])
    return out


old, new = load(sys.argv[1]), load(sys.argv[2])
ch, dn, up, moved = collections.Counter(), collections.Counter(), collections.Counter(), collections.Counter()
for k, (q, tier, *pos) in old.items():
    if k not in new:
        continue
    q2, tier2, *pos2 = new[k]
    if pos != pos2:
        moved[tier] += 1
    if q != q2:
        ch[tier] += 1
        dn[tier] += q >= 20 and q2 < 20
        up[tier] += q < 20 and q2 >= 20
print("records: old %d new %d, only in new %d, only in old %d" % (
    len(old), len(new), len(set(new) - set(old)), len(set(old) - set(new))))
for tier in sorted(set(ch) | set(moved)):
    print("%s: MAPQ changed %d (fell below 20: %d, rose to >= 20: %d); placement changed %d" % (
        tier, ch[tier], dn[tier], up[tier], moved[tier]))
