# Copyright 2018 by Ivan Antonov. All rights reserved.

from math import log10
import unittest

from run_triplexator import compute_pvalue


class Test(unittest.TestCase):
    def test_compute_pvalue(self):
        rna_len = 1595   # MEG3 length
        dna_len = 400
        prm_str = '-l 10 -e 10 -fr off'

        # True estimates were computed by R (see 181211.F1000.poisson/)
        # lambda <- predict(lambda_lm, data.frame(DNA_len = 400))
        # ppois(1, lambda, lower.tail = FALSE)
        # ppois(10, lambda, lower.tail = FALSE)
        # ppois(100, lambda, lower.tail = FALSE)
        true_log_pvalue = {
            1: -log10(0.1762126),
            10: -log10(5.949311e-10),
            100: -log10(3.539248e-173)}

        for num_tpx in [1, 10, 100]:
            pvalue = compute_pvalue(rna_len, dna_len, num_tpx, prm_str)
            self.assertTrue(0 < pvalue < 1)

            log_diff = abs(-log10(pvalue) - true_log_pvalue[num_tpx])
            self.assertTrue(log_diff < 0.01)

        # Make sure the function returns -1 for unsupported params
        self.assertEqual(compute_pvalue(100, dna_len, num_tpx, prm_str), -1)
        self.assertEqual(compute_pvalue(rna_len, dna_len, num_tpx, ''), -1)

