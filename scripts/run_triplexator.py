#!/usr/bin/env python3

# Copyright 2018 by Ivan Antonov. All rights reserved.

"""Run Triplexator for given RNA and DNA sequences and output all predictions."""

from pprint import pprint

import argparse
import logging
import os
import re
import shutil
import subprocess          # For triplexator run
import tempfile

from Bio import SeqIO
from scipy.stats import poisson


def main(args):
    if args.quiet:
        logging.getLogger().setLevel(logging.WARNING)
    else:
        logging.getLogger().setLevel(logging.INFO)

    # Create tmp folder for triplexator output files
    tmp_dir = tempfile.mkdtemp(prefix='__TPX.', dir=os.getcwd())
    prm_with_od = args.prm + ' -od ' + tmp_dir
    all_scores = run_triplexator(args.rna_fn, args.dna_fn, prm_with_od)
    summary = read_tpx_summary(os.path.join(tmp_dir, 'triplex_search.summary'))
    if not args.keep:
        shutil.rmtree(tmp_dir)    # remove tmp folder

    output(args.rna_fn, args.dna_fn, all_scores, summary, args.prm)

    # Move it to a separate script?
    if args.tpx_cov is not None:
        raise NotImplementedError("The function needs to be updated")
        save_tpx_coverage(all_tpx, summary, rna_dict, dna_dict, args)
    if args.all_tpx is not None:
        raise NotImplementedError("The function needs to be updated")
        save_all_tpx_info(all_tpx, args)

def output(rna_fn, dna_fn, all_scores, summary, prm_str):
    "Prints triplexator predictions for ALL(!) RNA-DNA pairs."
    head_arr = ['RNA', 'DNA', 'RNA_len', 'DNA_len',
                'Num_scores', 'Max_score', 'Sum_score',
                't_pot', 'Num_tpx_abs', 'pvalue']
    print('\t'.join(head_arr))
    for rna in SeqIO.parse(rna_fn, "fasta"):
        for dna in SeqIO.parse(dna_fn, "fasta"):
            to_print = {'RNA': rna.id, 'DNA': dna.id,
                        'RNA_len': len(rna.seq), 'DNA_len': len(dna.seq),
                        'Num_scores': 0, 'Max_score': 0, 'Sum_score': 0,
                        'Num_tpx_abs': 0, 't_pot': 0, 'pvalue': 1}
            info = summary.get(rna.id, {}).get(dna.id, None)
            if info is None:
                # No predictions for this RNA-DNA pair
                print('\t'.join([str(to_print[k]) for k in head_arr]))
                continue

            scores = all_scores.get(rna.id, {}).get(dna.id, [])
            to_print['Num_scores'] = len(scores)
            to_print['Max_score'] = max(scores)
            to_print['Sum_score'] = sum(scores)
            to_print['Num_tpx_abs'] = int(info['Total (abs)'])
            to_print['t_pot'] = info['Total (rel)']
            to_print['pvalue'] = '%.1e' % compute_pvalue(
                len(rna.seq), len(dna.seq), to_print['Num_tpx_abs'], prm_str)
            print('\t'.join([str(to_print[k]) for k in head_arr]))

def compute_pvalue(rna_len, dna_len, num_tpx, prm_str):
    if prm_str != '-l 10 -e 10 -fr off':
        logging.warning("Can't compute pvalue because Triplexator was run with "
                        "unsupported set of params '%s'" % prm_str)
        return -1

    # R: coef(lambda_lm)
    p_lambda = -0.6884666667 + 0.00053712*rna_len + 0.0006028061*dna_len

    # P(X = x) (i.e. PMF) is added because we want to compute P(X >= x)
    return poisson.pmf(num_tpx, p_lambda) + poisson.sf(num_tpx, p_lambda)

def run_triplexator(rna_fn, dna_fn, prm):
    # https://stackoverflow.com/a/18244485/310453
    cmd_str = 'triplexator ' + prm + ' -ss ' + rna_fn + ' -ds ' + dna_fn
    out_arr = subprocess.getoutput(cmd_str).splitlines()

    # get the header row & remove the comment char
    head_str = out_arr.pop(0).lstrip('# ')
    head_arr = head_str.split('\t')

    all_scores = {}
    for line in out_arr:
        vals = line.split('\t')
        if( len(vals) != len(head_arr) ):
            logging.warning('Something is wrong: len(vals) != len(head_arr)')
            continue
        tpx = dict(zip(head_arr, vals))
        score = int(tpx['Score'])

        rna_id, dna_id = tpx['Sequence-ID'], tpx['Duplex-ID']
        all_scores.setdefault(rna_id, {}).setdefault(dna_id, [])
        all_scores[rna_id][dna_id].append(score)

    return all_scores

def read_tpx_summary(fn):
    with open(fn) as f:
        # Get the header fields as list
        head_arr = f.readline().lstrip('# ').split('\t')
        summary = {}
        for line in f:
            vals = re.split('\s+', line)
            tpx = dict(zip(head_arr, vals))
            rna_id, dna_id = tpx['Sequence-ID'], tpx['Duplex-ID']
            if summary.get(rna_id, {}).get(dna_id, None) is not None:
                logging.warning("The key '%s' is duplicated in file '%s'" %
                                (key, fn))
                continue
            summary.setdefault(rna_id, {})[dna_id] = tpx
    return summary

def save_all_tpx_info(all_pairs, args):
    all_keys = ['Sequence-ID', 'TFO start', 'TFO end', 'Duplex-ID', 'TTS start', 'TTS end',
                'Score', 'Error-rate', 'Errors', 'Motif', 'Strand', 'Orientation', 'Guanine-rate']
    # 'TFO start'  =>  'TFO_start'
    rename_d = {'Sequence-ID': 'RNA', 'Duplex-ID': 'DNA'}
    colnames = [rename_d.get(k, re.compile(r'[ \-]+').sub('_', k)) for k in all_keys]

    out_fn = os.path.join(args.init_cwd, args.all_tpx)
    with open(out_fn, 'w') as f:
        f.write('\t'.join(colnames) + '\n')
        for pair_tpx in all_pairs.values():
            for tpx in pair_tpx:
                to_print = [tpx[k] for k in all_keys]
                f.write('\t'.join(to_print) + '\n')

def save_tpx_coverage(all_tpx, summary, rna_dict, dna_dict, args):
    for key in summary:
        fn_key = re.compile(r'[^\w\d]+').sub('_', key)   # '...chr18:657...'   => ''...chr18_657...'
        logging.info("Generating coverage file '%s'..." % fn_key)

        rna_out_fn = os.path.join(args.init_cwd, args.tpx_cov + fn_key + '.RNA.txt')
        _save_tpx_coverage_for_pair_and_type(key, 'RNA', summary, all_tpx, rna_dict, rna_out_fn)

        dna_out_fn = os.path.join(args.init_cwd, args.tpx_cov + fn_key + '.DNA.txt')
        _save_tpx_coverage_for_pair_and_type(key, 'DNA', summary, all_tpx, dna_dict, dna_out_fn)

def _save_tpx_coverage_for_pair_and_type(key, seq_type, summary, all_tpx, seq_dict, out_fn):
    s_info = summary[key]
    seq_name = s_info['Sequence-ID'] if seq_type == 'RNA' else s_info['Duplex-ID']

    # create .chrom.sizes files for RNA and DNA
    sizes_fn = seq_type + '.chrom.sizes'
    with open(sizes_fn, 'w') as f:
        f.write("%s\t%i\n" % (seq_name, len(seq_dict[seq_name])))

    bed_fn = seq_type + '.tpx.bed'
    with open(bed_fn, 'w') as f:
        if seq_type == 'RNA':
            for tpx in sorted(all_tpx[key], key=lambda k: int(k['TFO start'])):
                f.write("%s\t%s\t%s\n" % (seq_name, tpx['TFO start'], tpx['TFO end']))
        elif seq_type == 'DNA':
            for tpx in sorted(all_tpx[key], key=lambda k: int(k['TTS start'])):
                f.write("%s\t%s\t%s\n" % (seq_name, tpx['TTS start'], tpx['TTS end']))

    subprocess.run('bedtools genomecov -d -i ' + bed_fn + ' -g ' + sizes_fn + ' > ' + out_fn, shell=True)

def parse_args():
    parser = argparse.ArgumentParser(
        description='Run the Triplexator and process its output')
    parser.add_argument('rna_fn', metavar='RNA.fna',
                        help='corresponds to the -ss Triplexator option')
    parser.add_argument('dna_fn', metavar='DNA.fna',
                        help='corresponds to the -ds Triplexator option')
    parser.add_argument('--prm', metavar='STR', default='-l 10 -e 10 -fr off',
                        help="additional Triplexator params "
                       "(default is '-l 10 -e 10 -fr off')")
    parser.add_argument('--tpx_cov', metavar='PREFIX',
                        help='compute triplex coverage info for both RNA and DNA')
    parser.add_argument('--all_tpx', metavar='FN.txt',
                        help='save info about all the predicted triplexes')
    parser.add_argument('--keep', action='store_true',
                        help="do not remove tmp folder")
    parser.add_argument('-q', '--quiet', action='store_true',
                        help="do not write info messages to stderr")
    return parser.parse_args()

if __name__ == '__main__':
    main(parse_args())

