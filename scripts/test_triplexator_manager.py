# Copyright 2018 by Ivan Antonov. All rights reserved.

from math import log10
import unittest

from triplexator_manager import compute_pvalue


class Test(unittest.TestCase):
    def test_compute_pvalue(self):
        rna_len = 1595   # MEG3 length
        dna_len = 400
        prm_str = '-l 10 -e 10 -fr off'

        # True estimates were computed by R (see 190110.F1000.poisson/)
        # params <- data.frame(RNA_len = 1595, DNA_len = 400)
        # lambda <- predict(lambda_lm, params)
        # dpois(10, lambda) + ppois(10, lambda, lower.tail = FALSE)
        true_log_pvalue = {
            0: -log10(1),
            1: -log10(0.3359263),
            10: -log10(2.511582e-11),
            100: -log10(1.160769e-197)}

        for num_tpx in true_log_pvalue.keys():
            pvalue = compute_pvalue(rna_len, dna_len, num_tpx, prm_str)
            self.assertTrue(0 <= pvalue <= 1)

            log_diff = abs(-log10(pvalue) - true_log_pvalue[num_tpx])
            self.assertTrue(log_diff < 0.01)

        # Make sure the function returns -1 for unsupported params
        self.assertEqual(compute_pvalue(rna_len, dna_len, num_tpx, ''), -1)

