# A5 fit: phantom parameter sweep on HG002 (4M single-end subsample)

Ranked by calibration error (Sum|misplaced-expected| over MAPQ bands);
pi_unique=0 throughout (locus prior, A4). N=98%/95% are misplaced per
million among the N most confident mates with a truth.

| delta | pi | calib err | N=98% /M | N=95% /M | MAPQ>=40 err | mates |
|---|---|---|---|---|---|---|
| 0.006 | 0.003 | 1,756 | 130.8 | 54.9 | 190 | 3,201,304 |
| 0.003 | 0.003 | 1,766 | 135.1 | 54.9 | 190 | 3,201,304 |
| 0.003 | 0.01 | 1,773 | 134.5 | 54.9 | 190 | 3,201,304 |
| 0.006 | 0.01 | 1,804 | 129.0 | 54.9 | 190 | 3,201,304 |
| 0.01 | 0.003 | 1,903 | 135.7 | 54.9 | 190 | 3,201,304 |
| 0.01 | 0.01 | 2,039 | 134.5 | 54.9 | 190 | 3,201,304 |
| 0.01 | 0.03 | 2,956 | 128.7 | 54.9 | 190 | 3,201,304 |
| 0.006 | 0.03 | 3,281 | 145.7 | 54.9 | 190 | 3,201,304 |
| 0.003 | 0.03 | 3,360 | 146.6 | 54.9 | 190 | 3,201,304 |

**Frozen: delta=0.006, pi=0.003, pi_unique=0** (calib err 1,756).
