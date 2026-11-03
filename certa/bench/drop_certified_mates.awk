# Filters the fallback mapper's SAM in paired mode. CERTA sends both mates of
# a pair to the fallback even when one mate is already certified, renaming the
# pair "name~1" or "name~2" (the certified mate). This drops the fallback's
# records of the certified mate and restores the original name. With
# -v div=N it also divides the MAPQ and MQ:i of the records it keeps by N,
# rounding down: on the reads CERTA does not certify, minibwa's MAPQ is about
# four times the Phred of its error rate (bench/VALIDATION_RESULTS.md, "MAPQ
# of the uncertified mates"). POSIX awk.
#   minibwa map idx unc_1.fq unc_2.fq | awk -f bench/drop_certified_mates.awk -v div=4
BEGIN { FS = OFS = "\t" }
/^@/ { print; next }
{
  n = length($1)
  if (n > 2 && substr($1, n - 1, 1) == "~") {
    k = substr($1, n, 1)
    mate = (int($2 / 64) % 2) ? "1" : ((int($2 / 128) % 2) ? "2" : "")
    if (mate == k) next
    $1 = substr($1, 1, n - 2)
  }
  if (div > 1) {
    $5 = int($5 / div)
    for (i = 12; i <= NF; i++) if ($i ~ /^MQ:i:/) $i = "MQ:i:" int(substr($i, 6) / div)
  }
  print
}
