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

DEFAULT_PRM = '-l 10 -e 10 -fr off'


def main_pvalue(args):
    """Prints p-value for each RNA-DNA pair."""
    all_scores = run_triplexator(
        args.rna_fn, args.dna_fn, args.prm, args.tmp_dir, value_key='Score')
    summary = read_tpx_summary(os.path.join(args.tmp_dir, 'triplex_search.summary'))
    output(args.rna_fn, args.dna_fn, all_scores, summary, args.prm)

def main_coverage(args):
    """Prints output to STDOUT."""
    all_dna = list(SeqIO.parse(args.dna_fn, "fasta"))
    if len(all_dna) != 1:
        raise ValueError(
            "File '%s' contains '%s' sequences, but it is expected to "
            "contain 1 DNA only!" % (args.dna_fn, len(all_dna)))

    # Create tmp folder for triplexator output files
    all_tpx = run_triplexator(
        args.rna_fn, args.dna_fn, args.prm, args.tmp_dir)

    poly_ga_dict = _find_poly_purine_regions(
        all_dna[0].seq, args.strand)

    merged_bed_fn = _create_merged_bed(
        all_tpx, all_dna[0].id, args.strand, args.tmp_dir)

    _print_tpx_coverage(
        all_dna[0], merged_bed_fn, poly_ga_dict, args.strand, args.tmp_dir)

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

def run_triplexator(rna_fn, dna_fn, prm, tmp_dir, value_key=None):
    """Returns a dict of dicts: RNA_ID -> DNA_ID -> LIST."""
    # https://stackoverflow.com/a/18244485/310453
    prm_with_od = prm + ' -od ' + tmp_dir
    cmd_str = ' '.join([
        'triplexator', prm_with_od, '-ss', rna_fn, '-ds', dna_fn])
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
        tpx['Score'] = int(tpx['Score'])

        rna_id, dna_id = tpx['Sequence-ID'], tpx['Duplex-ID']
        all_scores.setdefault(rna_id, {}).setdefault(dna_id, [])

        if value_key is None:
            # The value is a list of dicts
            value = tpx
        else:
            # The value is a list of simple objects, e.g. ints
            value = tpx[value_key]
        all_scores[rna_id][dna_id].append(value)

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

def _find_poly_purine_regions(dna_seq, strand):
    """Retruns a dict where keys are the dna_seq coordiantes (0-based) and
    the values are 1's indicating that these positions correspond to the 
    poly-purine elements.
    """
    dna_seq = dna_seq.upper()
    if strand == -1:
        dna_seq = dna_seq.complement()

    poly_ga_re = re.compile(r'[GA]{10,}')
    poly_ga_dict = {}
    for mo in re.finditer(poly_ga_re, str(dna_seq)):
        for coord in range(mo.start(), mo.end()):
            poly_ga_dict[coord] = 1

    return poly_ga_dict

def _create_merged_bed(all_tpx, dna_id, strand, tmp_dir):
    # We want: 'chr19:13618723-13619244   13   23   ENST00000637775.1   9   +'
    bed_keys = ['Duplex-ID', 'TTS start', 'TTS end', 'Sequence-ID', 'Score',
                'Strand']
    strand_char = '+' if strand == 1 else '-'

    merged_bed_fn = os.path.join(tmp_dir, 'all_tpx_merged.bed')
    out_f = open(merged_bed_fn, 'w')
    for rna_id in sorted(all_tpx.keys()):
        bed_txt = ''
        for tpx in all_tpx[rna_id][dna_id]:
            if tpx['Strand'] != strand_char:
                continue
            vals = [str(tpx[k]) for k in bed_keys]
            bed_txt += "\t".join(vals) + "\n"

        # Merge triplex coordinates from the same RNA
        merge_cmd = 'sort -k1,1 -k2,2n | bedtools merge -s'
        merged_bed_txt = subprocess.run(
            merge_cmd, shell=True, stdout=subprocess.PIPE,
            input=bed_txt, universal_newlines=True
        ).stdout

        # Add RNA name and append it to file
        for line in merged_bed_txt.splitlines():
            # 'chr19:13618723-13619244  28  55'  =>  '...  ENST00000479263.5  +'
            out_f.write("\t".join([line, rna_id, strand_char]) + "\n")
    out_f.close()

    return merged_bed_fn

def _print_tpx_coverage(dna_record, merged_bed_fn, poly_ga_dict,
                        strand, tmp_dir):
    # create .chrom.sizes file (just 1 line)
    sizes_fn = os.path.join(tmp_dir, 'target_seq.chrom.sizes')
    with open(sizes_fn, 'w') as f:
        f.write("%s\t%i\n" % (dna_record.id, len(dna_record.seq)))

    # Sort the .bed file
    sorted_bed_fn = merged_bed_fn + '.sorted'
    sort_cmd = "sort -k1,1 -k2,2n %s  >  %s" % (
        merged_bed_fn, sorted_bed_fn)
    subprocess.run(sort_cmd, shell=True)

    coverage_cmd = 'bedtools genomecov -d -i %s -g %s' % (
        sorted_bed_fn, sizes_fn)
    coverage_txt = subprocess.run(
        coverage_cmd, shell=True, stdout=subprocess.PIPE,
        universal_newlines=True
    ).stdout

    # print header
    head_arr = ['name', 'strand', 'coord', 'nt', 'poly_ga', 'num_rna']
    print("\t".join(head_arr))

    dna_seq = dna_record.seq.upper()
    if strand == -1:
        dna_seq = dna_seq.complement()

    for line in coverage_txt.splitlines():
        # line = 'chr19:13618723-13619244   442   1'
        name, coord, num_rna = line.split()
        coord = int(coord)
        to_print = {'name': name, 'coord': coord, 'num_rna': num_rna,
                    'strand': '+' if strand == 1 else '-',
                    'nt': dna_seq[coord-1],
                    'poly_ga': poly_ga_dict.get(coord-1, 0)}
        print('\t'.join([str(to_print[k]) for k in head_arr]))

def _add_common_arguments(parser):
    parser.add_argument('rna_fn', metavar='RNA.fna',
                        help='corresponds to the -ss Triplexator option')
    parser.add_argument('dna_fn', metavar='DNA.fna',
                        help='corresponds to the -ds Triplexator option')
    parser.add_argument('--prm', metavar='STR', default=DEFAULT_PRM,
                        help="additional Triplexator params "
                       "(default is '%s')" % DEFAULT_PRM)
    parser.add_argument('--keep', action='store_true',
                        help="do not remove tmp folder")
    parser.add_argument('-q', '--quiet', action='store_true',
                        help="do not write info messages to stderr")

def parse_args():
    """Parse command line arguments using ArgumentParser
    (https://docs.python.org/3.5/library/argparse.html#sub-commands)."""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(metavar='<subcommand>')
    subparsers.required = True   # https://stackoverflow.com/a/18283730/310453

    # pvalue command
    parser_pvalue = subparsers.add_parser(
        'pvalue', help="run Triplexator and compute p-value for each "
        "RNA-DNA pair")
    parser_pvalue.set_defaults(func=main_pvalue)
    _add_common_arguments(parser_pvalue)

    # coverage command
    parser_coverage = subparsers.add_parser(
        'coverage', help="For the given DNA sequence (1 only!) compute "
        "its coverage by triplexes from different RNAs")
    parser_coverage.set_defaults(func=main_coverage)
    _add_common_arguments(parser_coverage)
    parser_coverage.add_argument(
        'strand', type=int, choices=[1, -1], metavar='STRAND',
        help='DNA strand to be considered')

    # parse the args
    args = parser.parse_args()

    if args.quiet:
        logging.getLogger().setLevel(logging.WARNING)
    else:
        logging.getLogger().setLevel(logging.INFO)

    # Create tmp folder for triplexator output files
    args.tmp_dir = tempfile.mkdtemp(prefix='__TPX.', dir=os.getcwd())
    args.func(args)   # call whatever function was selected
    if not args.keep:
        shutil.rmtree(args.tmp_dir)    # remove tmp folder

if __name__ == '__main__':
    parse_args()

