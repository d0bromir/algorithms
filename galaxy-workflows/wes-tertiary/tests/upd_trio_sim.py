"""Uniparental disomy and parent-of-origin validation of the clinical report on a real trio (README section 7).

    python3 tests/upd_trio_sim.py tools/wes_clinical_report/wes_clinical_report.py NAME SEX child.vcf         child_thresholds.bed.gz mother.vcf father.vcf [child_roh.txt] [--split tune|holdout|all]

child.vcf is the VCF workflow's annotated VCF, child_thresholds the mosdepth per-target thresholds, the parents'
VCFs their normalised calls, child_roh the bcftools roh output. Standard library only.

Real trio: the best segment per chromosome and parent (false-positive background). Simulations replace the
proband's genotypes on one chromosome (or a segment) with genotypes drawn from the parents' calls:
  iso:  homozygous for one allele of the donor parent at every site (unphased: a random allele per site)
  het:  the donor parent's genotype
  tri:  donor's two alleles + one allele of the other parent (allele fraction from reads, DP 40, binomial)
Sites outside the proband's >= 20x targets are dropped (an exome would not call them). Odd chromosomes and X tune,
even chromosomes are held out."""
import copy
import importlib.util
import random
import resource
import sys
import time

spec = importlib.util.spec_from_file_location("r", sys.argv[1])
R = importlib.util.module_from_spec(spec)
spec.loader.exec_module(R)
name, sex, cvcf, tcov, mvcf, fvcf = sys.argv[2:8]
roh_path = sys.argv[8] if len(sys.argv) > 8 and not sys.argv[8].startswith("--") else None
split = sys.argv[sys.argv.index("--split") + 1] if "--split" in sys.argv else "all"
t0 = time.time()
_s, _c, variants, prefix = R.parse_vcf(cvcf)
covered = R.Covered(R.load_thresholds(tcov))
trio = R.Trio(mvcf, fvcf, None, None, prefix, "none")
roh = R.load_roh(roh_path)
t_load = time.time() - t0
rng = random.Random(7)

t1 = time.time()
segs, rates = R.upd_scan(variants, trio, sex, covered, (), roh)
t_scan = time.time() - t1
print(f"{name} ({sex}): load {t_load:.1f} s, scan {t_scan:.1f} s, peak "
      f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2:.2f} GiB; error rates "
      + ", ".join(f"{w} {r:.4f}" for w, r in rates.items()))
sites = R.upd_sites(variants, trio, sex, covered)
bg = []
for (c, who), s in sorted(sites.items()):
    lod, i, j = R.best_segment([x for _p, x in s], rates[who])
    bg.append((lod, c, who, sum(x for _p, x in s[i:j + 1]) if j >= i else 0, len(s)))
bg.sort(reverse=True)
print("  real trio, called segments:", [(d["chrom"], d["from"], round(d["lod"], 1)) for d in segs] or "none")
print("  background, top 5 segment LODs:", [(c, w, round(lod, 1), ex, n) for lod, c, w, ex, n in bg[:5]])
print(f"  informative sites per chromosome and parent: median "
      f"{sorted(len(s) for s in sites.values())[len(sites) // 2]}")

AFS = {(v.chrom, v.pos, v.ref, v.alt): v.af for v in variants}


def mk(c, p, r, a, zyg, af, ab):
    v = R.Variant()
    v.chrom, v.pos, v.ref, v.alt, v.zyg, v.af, v.qual, v.dp, v.ab = c, p, r, a, zyg, af, 60.0, 40, ab
    v.anns, v.gene = [], ""
    return v


def alleles(call):
    """parent genotype at an ALT site as (alt copies of 2)"""
    return 2 if call and call[0] == "hom" else (1 if call else 0)


def simulate(kind, donor, chrom, start=0, end=10 ** 10):
    other = R.OTHER[donor]
    keep = [v for v in variants if not (v.chrom == chrom and start <= v.pos <= end)]
    keys = set()
    for w in ("mother", "father"):
        keys |= {k for k in trio.vcf[w] if k[0] == chrom and start <= k[1] <= end}
    keys |= {(v.chrom, v.pos, v.ref, v.alt) for v in variants if v.chrom == chrom and start <= v.pos <= end}
    for k in sorted(keys):
        if len(k[2]) != 1 or len(k[3]) != 1 or (k[0], k[1]) not in covered:
            continue
        d, o = alleles(trio.vcf[donor].get(k)), alleles(trio.vcf[other].get(k))
        if kind == "iso":
            a = 2 * (rng.random() < d / 2)
            n = 2
        elif kind == "het":
            a, n = d, 2
        else:  # tri: both donor alleles plus one random allele of the other parent
            a, n = d + (rng.random() < o / 2), 3
        if a == 0:
            continue
        reads = sum(rng.random() < a / n for _ in range(40))
        zyg = "hom" if a == n else "het"
        keep.append(mk(*k, zyg, AFS.get(k, 0.05), reads / 40))
    return keep


def half(c):
    return "tune" if c == "X" or int(c) % 2 else "holdout"


events = [("15", None, None), ("7", None, None), ("14", None, None), ("6", None, None), ("11", None, None),
          ("20", None, None), ("16", None, None), ("1", None, None), ("2", None, None), ("22", None, None),
          ("11", 0, 21_000_000), ("15", 22_000_000, 30_000_000), ("6", 140_000_000, 171_000_000),
          ("14", 95_000_000, 107_000_000), ("20", 50_000_000, 65_000_000)]
print("  simulated UPD (chromosome, segment, donor, kind): called / LOD / kind / sites")
for c, s, e in events:
    if split != "all" and half(c) != split:
        continue
    for donor in ("mother", "father"):
        for kind in ("iso", "het"):
            sim = simulate(kind, donor, c, s or 0, e or 10 ** 10)
            got, _ = R.upd_scan(sim, trio, sex, covered, (), roh)
            hit = [d for d in got if d["chrom"] == c and d["from"] == donor]
            extra = [d for d in got if not (d["chrom"] == c and d["from"] == donor)]
            seg = f"{(s or 0) / 1e6:.0f}-{e / 1e6:.0f} Mb" if e else "whole"
            h = hit[0] if hit else None
            print(f"    chr{c:3} {seg:12} {donor:6} {kind}: "
                  + (f"FOUND LOD {h['lod']:.0f}, {h['kind']}, {h['excluded']}/{h['sites']} sites, "
                     f"{h['start'] / 1e6:.1f}-{h['end'] / 1e6:.1f} Mb" if h else "missed")
                  + (f"; other calls {[(d['chrom'], d['from']) for d in extra]}" if extra else ""))
print("  simulated trisomy: extra copy from / maternal share / sites")
for c in ("21", "18", "13", "16", "22", "15"):
    if split != "all" and half(c) != split:
        continue
    for donor in ("mother", "father"):
        sim = simulate("tri", donor, c)
        print(f"    chr{c:3} {donor:6}: " + R.gain_origin_text(R.gain_origin(sim, trio, c)))
    print(f"    chr{c:3} none  : " + R.gain_origin_text(R.gain_origin(variants, trio, c)))
print("  X parents (real):", R.x_parents_text(R.x_parents(variants, trio, covered)))
