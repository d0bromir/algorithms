#include "map_config.h"

#include <algorithm>
#include <cstdio>
#include <stdexcept>

namespace certa {

MapConfig MapConfig::from_args(const Args& a) {
  if (a.pos.size() != 2 && a.pos.size() != 3)
    throw std::runtime_error(
        "usage: certa map ref.cidx reads.fq[.gz] [mates.fq[.gz]] [-k 5] [--k1 2] [--budget 16] [--budget2 256] "
        "[--budget3 0] [--support3 2] [--local 0] [--full-cert] [--adaptive-parts 0|1|2] [--vbudget 4096] [--index2 second.cidx] [--mapq-limit 8] [-t N] [--gpu] [--device 0] [-o out.sam] "
        "[-u uncertified.fq] [-U uncertified_mates.fq] [--fallback-mates pair|single -S single_mates.fq] [--max-dist 850] [--stats s.json] [--batch N] [--io-threads N] [--serial-stages] "
        "[--bio [--bio-snv 0.001] [--bio-titv 2.1] [--bio-indel 0.0001] [--bio-seq-indel 0.00001] [--bio-gap-ext 0.5] "
        "[--bio-str-min 4] [--bio-str-slope 3] [--bio-str-cap 24] [--bio-gap-cap 50] [--bio-delta 60] [--bio-max-band 1024] "
        "[--bio-budget 4096] [--bio-insert 400,100] [--bio-improper 40] [--bio-rescue-loci 8]]");
  MapConfig c;
  c.index_path = a.pos[0];
  c.reads_path = a.pos[1];
  // Paired-end: a second FASTQ with the mates, in the same order.
  if (a.pos.size() == 3) c.mates_path = a.pos[2];

  // Pass 1 runs every read with the small budget and without S2 (uniform,
  // fast GPU warps). Pass 2 re-runs only the reads pass 1 could not certify,
  // as a compacted batch with --budget2 and S2, so the expensive work stays
  // among hard reads. Each re-run read gets exactly its single-pass result.
  Params& p = c.p1;
  p.k = a.geti("-k", 5);
  p.budget = a.geti("--budget", a.geti("--cap", 16));  // --cap: the v0.1 name
  if (p.k < 0 || p.k > KMAX) throw std::runtime_error("-k must be in [0, 5]");
  if (p.budget < 1 || p.budget > GPU_BUDGET_MAX) throw std::runtime_error("--budget must be in [1, 256]");
  // MAPQ only: the second-best estimate beyond the certified radius, up to
  // this many edits (0 = off).
  p.mapq_limit = a.geti("--mapq-limit", 8);
  // S2 is off by default: on synthetic truth most S2 placements with MAPQ >= 20
  // were wrong, and on HG002 S2 reads were enriched at false-positive calls.
  p.s2_limit = a.geti("--s2", 0);
  if (p.s2_limit < 0 || p.s2_limit > S2_MAX) throw std::runtime_error("--s2 must be in [0, 10]");
  c.out.s2_limit = p.s2_limit;
  Params& p2 = c.p2;
  p2 = p;
  p2.budget = a.geti("--budget2", 256);
  if (p2.budget != 0 && (p2.budget < p.budget || p2.budget > GPU_BUDGET_MAX))
    throw std::runtime_error("--budget2 must be 0 or in [budget, 256]");
  if (p2.budget == 0) p2.budget = p.budget;
  p.s2_limit = 0;
  // Pass 1 also uses a smaller radius (--k1): wide bands for R up to k are
  // only needed by the few reads pass 1 cannot certify.
  p.k = std::min(p2.k, a.geti("--k1", 2));
  if (p.k < 0) throw std::runtime_error("--k1 must be >= 0");

  // --full-cert: the most certified reads the certificates allow, at a cost
  // in speed: pass 3 with budget 4096, pass V (adaptive parts) and pass L
  // with t = 2 unless given.
  c.full_cert = a.has("--full-cert");
  // Pass 3 (CPU, opt-in): larger budgets for reads whose parts are all
  // repetitive. It certifies ~0.9% more reads but costs more CPU than the
  // fallback mapper spends on them, so it is off by default; use
  // --budget3 4096 (or --full-cert) for maximum certification.
  Params& p3 = c.p3;
  p3 = p2;
  p3.budget = a.geti("--budget3", c.full_cert ? HOST_BUDGET_MAX : 0);
  if (p3.budget != 0 && (p3.budget < p2.budget || p3.budget > HOST_BUDGET_MAX))
    throw std::runtime_error("--budget3 must be 0 or in [budget2, 4096]");
  if (p3.budget == 0) p3.budget = p2.budget;
  // Pass 3 verifies only clusters hit by >= t distinct parts (radius |S| - t),
  // which skips the flood of single-part hits from repetitive parts.
  p3.min_support = a.geti("--support3", 2);
  if (p3.min_support < 1 || p3.min_support > PMAX) throw std::runtime_error("--support3 must be in [1, 16]");
  // Pass L (CPU, opt-in): certified local alignment (tier SL) for reads still
  // uncertified. --local t verifies chains hit by >= t parts; 2 is the
  // cheap setting (t = 1 costs ~10x more on HG002 for +0.16 points).
  Params& pl = c.pl;
  pl = p2;
  pl.budget = std::max(p2.budget, p3.budget);
  pl.local = a.geti("--local", c.full_cert ? 2 : 0);
  pl.local_only = 1;
  if (pl.local < 0 || pl.local > PMAX) throw std::runtime_error("--local must be in [0, 16]");

  PairParams& pp = c.pp;
  pp.se = p2;
  pp.se.budget = std::max(p2.budget, p3.budget);
  pp.se.min_support = 1;
  pp.max_dist = a.geti("--max-dist", 850);
  if (pp.max_dist < 0) throw std::runtime_error("--max-dist must be >= 0");
  c.out.max_dist = pp.max_dist;
  c.out.unc_tags = a.has("--unc-tags");

  // Interval output: no read goes to the fallback; tier R (CPU) for every
  // read still uncertified, pairs placed by certa/interval.h.
  c.intervals = a.has("--intervals");
  RepeatParams& rp = c.rp;
  rp.tau = a.geti("--iv-tau", 256);
  rp.cmax = a.geti("--iv-cmax", 4096);
  rp.samples = a.geti("--iv-samples", 8);
  rp.sub_margin = a.geti("--iv-margin", 20);
  rp.fallback_samples = a.geti("--iv-fallback", 8);
  rp.lb_chains = a.geti("--iv-chains", 8);
  if (rp.tau < 0 || rp.cmax < 1 || rp.samples < 1 || rp.sub_margin < 0 || rp.fallback_samples < 0 ||
      rp.lb_chains < 0)
    throw std::runtime_error("bad --iv-* parameters");
  const std::string mapq_mode = a.get("--iv-mapq", "point");
  if (mapq_mode != "point" && mapq_mode != "cons") throw std::runtime_error("--iv-mapq must be point or cons");
  c.out.mapq_cons = mapq_mode == "cons";
  c.out.witness = a.has("--witness");
  if (c.out.witness && !c.intervals) throw std::runtime_error("--witness needs --intervals");

  // --index2: a second index of the same reference (e.g. q = 15: more, shorter
  // parts) on which the reads still uncertified after every pass get every
  // pass again (CPU); a certificate found there is a certificate.
  c.index2_path = a.get("--index2", "");
  c.serial_stages = a.has("--serial-stages");
  if (c.index2() && c.out.witness)
    throw std::runtime_error("--witness with --index2: the witnesses would refer to two indexes");
  if (c.index2() && c.serial_stages)
    throw std::runtime_error("--index2 needs the fused stages (no --serial-stages)");

  // Pass V (CPU; --full-cert unless --witness): parts chosen per read on the
  // main index, the rarest disjoint q-mers (1), then also q-mers searched with
  // up to one edit (2), for the reads passes 1-3 leave uncertified; a
  // certificate found there is taken (core.h, certify_vparts).
  Params& pv = c.pv;
  pv = p2;
  pv.vparts = a.geti("--adaptive-parts", c.full_cert && !c.out.witness ? 2 : 0);
  if (pv.vparts < 0 || pv.vparts > 2) throw std::runtime_error("--adaptive-parts must be 0, 1 or 2");
  if (pv.vparts > 0 && c.out.witness)
    throw std::runtime_error("--witness with --adaptive-parts: no witness describes adaptive parts");
  pv.budget = a.geti("--vbudget", HOST_BUDGET_MAX);
  if (pv.budget < 1 || pv.budget > HOST_BUDGET_MAX) throw std::runtime_error("--vbudget must be in [1, 4096]");
  pv.k = VMAX;
  pv.min_support = 1;
  pv.s2_limit = 0;
  pv.local = pv.local_only = 0;

  // Outputs.
  if (c.intervals && (a.has("-u") || a.has("-U") || a.has("-S")))
    throw std::runtime_error("--intervals writes every read to SAM: no -u/-U/-S");
  c.sam_path = a.get("-o", "");
  c.fq_path = a.get("-u", "");
  c.fq2_path = a.get("-U", "");
  c.stats_path = a.get("--stats", "");
  // --fallback-mates single: the uncertified mate of a half-certified pair goes
  // single-end to -S (default "pair": both mates to -u/-U, renamed name~k).
  const std::string fb_mode = a.get("--fallback-mates", "pair");
  if (fb_mode != "pair" && fb_mode != "single") throw std::runtime_error("--fallback-mates must be pair or single");
  if (fb_mode == "single") {
    if (!a.has("-S")) throw std::runtime_error("--fallback-mates single needs -S single_mates.fq");
    c.single_path = a.get("-S", "");
    c.out.single_mates = true;
  }
  if (c.paired() && a.has("-u") && !a.has("-U")) throw std::runtime_error("paired input: -u needs -U for the mates");

  c.gpu = a.has("--gpu");
  c.device = a.geti("--device", 0);
  c.threads = a.geti("-t", default_threads());
  // Threads that parse and encode the next batch while the current one maps.
  c.io_threads = a.geti("--io-threads", std::max(1, std::min(16, c.threads / 4)));
  c.batch_size = static_cast<size_t>(a.geti("--batch", c.gpu ? 1000000 : 200000));

  // --bio: one CPU pass per read that certifies its maximum-likelihood
  // placement (certa/bio.h, THEORY.md section 12) under the model below;
  // pairs are placed jointly. The other passes are not run.
  c.bio = a.has("--bio");
  if (c.bio) {
    if (c.gpu) throw std::runtime_error("--bio runs on the CPU: drop --gpu");
    if (c.intervals) throw std::runtime_error("--bio does not combine with --intervals (its bounds are BWA-MEM scores)");
    if (c.index2()) throw std::runtime_error("--bio does not combine with --index2");
    if (c.serial_stages) throw std::runtime_error("--bio does not combine with --serial-stages");
  }
  BioParams& bp = c.bp;
  bp.snv = a.getd("--bio-snv", bp.snv);
  bp.titv = a.getd("--bio-titv", bp.titv);
  bp.indel = a.getd("--bio-indel", bp.indel);
  bp.seq_indel = a.getd("--bio-seq-indel", bp.seq_indel);
  bp.gap_ext = a.getd("--bio-gap-ext", bp.gap_ext);
  bp.str_min = a.geti("--bio-str-min", bp.str_min);
  bp.str_slope = a.geti("--bio-str-slope", bp.str_slope);
  bp.str_cap = a.geti("--bio-str-cap", bp.str_cap);
  bp.gap_cap = a.geti("--bio-gap-cap", bp.gap_cap);
  bp.delta = a.geti("--bio-delta", bp.delta);
  bp.max_band = a.geti("--bio-max-band", bp.max_band);
  bp.improper = a.geti("--bio-improper", bp.improper);
  bp.rescue_loci = a.geti("--bio-rescue-loci", bp.rescue_loci);
  const std::string ins = a.get("--bio-insert", "");
  if (!ins.empty()) {
    const size_t comma = ins.find(',');
    if (comma == std::string::npos) throw std::runtime_error("--bio-insert takes MEAN,SD");
    bp.insert_mean = std::atof(ins.substr(0, comma).c_str());
    bp.insert_sd = std::atof(ins.substr(comma + 1).c_str());
  }
  if (!(bp.snv > 0 && bp.snv < 0.5) || !(bp.titv > 0) || !(bp.indel > 0 && bp.indel < 0.5) ||
      !(bp.seq_indel >= 0 && bp.seq_indel < 0.5) || !(bp.gap_ext > 0 && bp.gap_ext < 1))
    throw std::runtime_error("--bio-snv/--bio-indel must be in (0, 0.5), --bio-seq-indel in [0, 0.5), "
                             "--bio-gap-ext in (0, 1), --bio-titv > 0");
  if (bp.str_min < 1 || bp.str_slope < 0 || bp.str_cap < 0 || bp.gap_cap < 0 || bp.gap_cap > 200 || bp.delta < 0 ||
      bp.max_band < 1 || bp.improper < 0 || bp.rescue_loci < 0 || !(bp.insert_sd > 0) || !(bp.insert_mean > 0))
    throw std::runtime_error("bad --bio-* parameters");
  Params& pb = c.pb;
  pb = c.p2;
  pb.k = KMAX;
  pb.budget = a.geti("--bio-budget", HOST_BUDGET_MAX);
  if (pb.budget < 1 || pb.budget > HOST_BUDGET_MAX) throw std::runtime_error("--bio-budget must be in [1, 4096]");
  pb.min_support = 1;
  pb.s2_limit = 0;
  pb.mapq_limit = 0;
  pb.local = pb.local_only = 0;
  pb.vparts = 0;
  return c;
}

void MapConfig::fit_index(int s) {
  if (bio) return;  // pass V is not run
  if (pass_v() && s != 1) {
    std::fprintf(stderr, "[map] pass V needs an index of every position (s = 1); off\n");
    pv.vparts = 0;
  }
}

std::string MapConfig::command_line() const {
  if (bio) {
    char b[512];
    std::snprintf(b, sizeof b,
                  "certa map --bio --bio-budget %d --bio-snv %g --bio-titv %g --bio-indel %g --bio-seq-indel %g "
                  "--bio-gap-ext %g --bio-str-min %d --bio-str-slope %d --bio-str-cap %d --bio-gap-cap %d "
                  "--bio-delta %d --bio-max-band %d --bio-insert %g,%g --bio-improper %d --bio-rescue-loci %d --max-dist %d",
                  pb.budget, bp.snv, bp.titv, bp.indel, bp.seq_indel, bp.gap_ext, bp.str_min, bp.str_slope,
                  bp.str_cap, bp.gap_cap, bp.delta, bp.max_band, bp.insert_mean, bp.insert_sd, bp.improper,
                  bp.rescue_loci, pp.max_dist);
    return b;
  }
  // The certification passes beyond the defaults, when on (the default
  // command line is unchanged).
  std::string extra;
  if (pass3()) extra += " --budget3 " + std::to_string(p3.budget) + " --support3 " + std::to_string(p3.min_support);
  if (pass_l()) extra += " --local " + std::to_string(pl.local);
  if (pass_v()) extra += " --adaptive-parts " + std::to_string(pv.vparts) + " --vbudget " + std::to_string(pv.budget);
  if (index2()) extra += " --index2 " + index2_path;
  return "certa map -k " + std::to_string(p2.k) + " --k1 " + std::to_string(p1.k) + " --budget " +
         std::to_string(p1.budget) + " --budget2 " + std::to_string(p2.budget) + " --s2 " +
         std::to_string(p2.s2_limit) + (gpu ? " --gpu" : "") + extra;
}

}  // namespace certa
