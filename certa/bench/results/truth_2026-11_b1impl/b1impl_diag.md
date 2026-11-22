# B1-impl diagnostics (500 k pairs, ~900 k commit_lift calls)

Temporary counters in `commit_lift` (since removed), exit-reason breakdown:

```
call=100000 early=12109  ct=0 here=88799  lo=0 align=0 score=245  commit=0
call=200000 early=23987  ct=0 here=177912 lo=0 align=0 score=483  commit=0
call=300000 early=36284  ct=0 here=266527 lo=0 align=0 score=746  commit=0
call=400000 early=48405  ct=0 here=355330 lo=0 align=0 score=1002 commit=0
call=500000 early=60542  ct=0 here=444095 lo=0 align=0 score=1263 commit=0
call=600000 early=72963  ct=0 here=532540 lo=0 align=0 score=1498 commit=0
call=700000 early=85374  ct=0 here=621058 lo=0 align=0 score=1716 commit=1
call=800000 early=97326  ct=0 here=710107 lo=0 align=0 score=1915 commit=1
call=900000 early=109518 ct=0 here=798992 lo=0 align=0 score=2147 commit=2
```

- `early`: lift state != 1, MAPQ < 40, L > LMAX, or mate not certified (~12%).
- `here`: lift within 150 bp of the GRCh38 placement, same contig/strand (~89%) —
  the chain projects the CHM13 copy back to the collapsed GRCh38 locus.
- `lo`/`align`: window invalid / affine_align failed (0).
- `score`: re-aligned elsewhere but score <= GRCh38 placement, deferred (~0.24%).
- `commit`: 2 in 900 k.

`pair_interval` confirmed `t2t_commit=1` and a non-null lift pointer (the commit
loop runs); the near-zero commit count is the biological/chain limitation, not a
wiring bug.
