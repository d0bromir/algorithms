# Filters the fallback mapper's SAM in paired mode. CERTA sends both mates of
# a pair to the fallback even when one mate is already certified, renaming the
# pair "name~1" or "name~2" (the certified mate). This drops the fallback's
# records of the certified mate and restores the original name. POSIX awk.
#   minibwa map idx unc_1.fq unc_2.fq | awk -f bench/drop_certified_mates.awk
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
  print
}
