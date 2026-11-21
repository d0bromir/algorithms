#!/usr/bin/env python3
"""Unit tests for truth_ceiling.py (A1, the offline gap-graded MAPQ ceiling).

  python3 -m unittest test_truth_ceiling -v   (from certa/bench)
"""
import io
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import truth_ceiling as tc  # noqa: E402


def sam(qname, flag, tags=""):
    return f"{qname}\t{flag}\tchr1\t100\t60\t150M\t*\t0\t0\t{'A'*150}\t{'I'*150}" + (("\t" + tags) if tags else "")


class TestPhred(unittest.TestCase):
    def test_bounds_and_monotone(self):
        self.assertEqual(tc.phred(0, 0), 3)  # round(-10*log10(0.5/1)) = 3
        for n in (10, 1000, 100000):
            qs = [tc.phred(e, n) for e in range(0, n + 1, max(1, n // 20))]
            self.assertTrue(all(0 <= q <= 60 for q in qs))
            self.assertTrue(all(a >= b for a, b in zip(qs, qs[1:])))  # more errors -> lower MAPQ

    def test_no_errors_high(self):
        self.assertGreater(tc.phred(0, 100000), 40)
        self.assertLess(tc.phred(50000, 100000), 5)


class TestPava(unittest.TestCase):
    def test_monotone_non_increasing(self):
        rng = np.random.default_rng(0)
        for _ in range(20):
            counts = rng.integers(1, 1000, 50)
            errors = rng.integers(0, counts + 1)
            cp, ep = tc.pava_decreasing(counts, errors)
            rate = ep / cp
            self.assertTrue(np.all(rate[:-1] >= rate[1:] - 1e-12))

    def test_preserves_weighted_total(self):
        counts = np.array([10, 20, 30, 40])
        errors = np.array([1, 15, 2, 20])  # violators everywhere
        cp, ep = tc.pava_decreasing(counts, errors)
        rate = ep / cp  # block rate replicated per bin
        # PAVA conserves the count-weighted error total: sum(counts * rate) = sum(errors)
        self.assertAlmostEqual(float(np.sum(counts * rate)), float(errors.sum()))
        self.assertTrue(np.all(rate[:-1] >= rate[1:] - 1e-12))

    def test_already_monotone_unchanged(self):
        counts = np.array([100, 100, 100])
        errors = np.array([50, 20, 1])
        cp, ep = tc.pava_decreasing(counts, errors)
        np.testing.assert_array_equal(cp, counts)
        np.testing.assert_array_equal(ep, errors)


class TestQmap(unittest.TestCase):
    def test_non_decreasing_in_gap(self):
        counts = np.zeros(tc.GAPMAX + 1, np.int64)
        errors = np.zeros(tc.GAPMAX + 1, np.int64)
        counts[2], errors[2] = 1000, 200   # 20% error at small gap
        counts[100], errors[100] = 1000, 1  # ~0.1% at large gap
        q = tc.build_qmap(counts, errors)
        self.assertLess(q[2], q[100])
        self.assertTrue(np.all(q[:-1] <= q[1:]))  # larger gap -> >= MAPQ
        self.assertTrue(np.all(q <= 60))

    def test_empty(self):
        q = tc.build_qmap(np.zeros(tc.GAPMAX + 1, np.int64), np.zeros(tc.GAPMAX + 1, np.int64))
        self.assertTrue(np.all(q == 60))


class TestAtN(unittest.TestCase):
    def test_pro_rata_boundary(self):
        mq = np.array([60, 60, 40, 40, 40])
        bad = np.array([1, 0, 1, 0, 0], bool)
        # n=3: two above (60s, 1 bad) + 1 of the three 40s (1/3 bad) = 1 + 1/3
        self.assertAlmostEqual(tc.at_n(mq, bad, 3), 1 + 1 / 3)

    def test_all(self):
        mq = np.array([10, 20, 30])
        bad = np.array([1, 1, 0], bool)
        self.assertEqual(tc.at_n(mq, bad, 5), 2.0)

    def test_top_only(self):
        mq = np.array([60, 50, 40])
        bad = np.array([0, 1, 1], bool)
        self.assertEqual(tc.at_n(mq, bad, 1), 0.0)


class TestReadSamGap(unittest.TestCase):
    def test_extraction(self):
        lines = [
            sam("a0", 0, "AS:i:100\tXS:i:90"),    # gap 10
            sam("a0", 128, "AS:i:100"),           # no XS -> GAPMAX
            sam("a1", 4),                          # unmapped -> -1
            sam("a1", 128, "AS:i:80\tXS:i:80"),   # gap 0
            sam("a0", 256, "AS:i:50\tXS:i:50"),   # secondary, skipped
            "@HD\tVN:1.6",                         # header, skipped
        ]
        gap = tc.read_sam_gap(io.StringIO("\n".join(lines) + "\n"), 2)
        self.assertEqual(gap[0, 0], 10)
        self.assertEqual(gap[0, 1], tc.GAPMAX)
        self.assertEqual(gap[1, 0], -1)
        self.assertEqual(gap[1, 1], 0)

    def test_negative_gap_clipped(self):
        gap = tc.read_sam_gap(io.StringIO(sam("a0", 0, "AS:i:50\tXS:i:90") + "\n"), 1)
        self.assertEqual(gap[0, 0], 0)


class TestRerankCeiling(unittest.TestCase):
    def test_small_gap_errors_demoted(self):
        counts = np.zeros(tc.GAPMAX + 1, np.int64)
        errors = np.zeros(tc.GAPMAX + 1, np.int64)
        counts[2], errors[2] = 2, 2     # gap 2: all misplaced -> MAPQ ~1
        counts[100], errors[100] = 2, 0  # gap 100: none misplaced -> MAPQ ~8
        qmap = tc.build_qmap(counts, errors)
        mq = np.array([60, 60, 60, 60])
        bad = np.array([1, 1, 0, 0], bool)
        gap = np.array([2, 2, 100, 100])
        gq = np.where(gap < 0, 60, qmap[np.clip(gap, 0, tc.GAPMAX)])
        rerank = np.minimum(mq, gq)
        self.assertTrue(np.all(rerank <= mq))
        # baseline: all MAPQ 60 tie, top 2 carries 1.0 misplaced pro rata;
        # re-ranked: the two correct (gap 100) outrank the two errors (gap 2) -> 0.
        self.assertAlmostEqual(tc.at_n(mq, bad, 2), 1.0)
        self.assertAlmostEqual(tc.at_n(rerank, bad, 2), 0.0)

    def test_no_competitor_no_penalty(self):
        qmap = tc.build_qmap(np.zeros(tc.GAPMAX + 1, np.int64), np.zeros(tc.GAPMAX + 1, np.int64))
        mq = np.array([60, 30])
        gap = np.array([-1, -1])  # competitor did not map these mates
        gq = np.where(gap < 0, 60, qmap[np.clip(gap, 0, tc.GAPMAX)])
        np.testing.assert_array_equal(np.minimum(mq, gq), mq)


if __name__ == "__main__":
    unittest.main()
