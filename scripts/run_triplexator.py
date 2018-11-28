#!/usr/bin/env python3

# $Id: run_triplexator.py 2943 2018-10-22 16:27:18Z antonov $

###
# Ivan Antonov (vanya.antonov@gmail.com)
#

import sys
import os
import re
import argparse
from pprint import pprint
from warnings import warn
from math import log10

import logging    # https://www.youtube.com/watch?v=-RcDmGNSuvU
logging.basicConfig(level=logging.INFO)
#logging.basicConfig(level=logging.DEBUG)

import tempfile, shutil    # To create tmp dir: https://stackoverflow.com/a/3223615/310453
import subprocess          # For triplexator run

from mylib.baseutilbio import fasta2dict

###
# SUBROUTINES
def run(args):
    all_tpx = run_triplexator(args)
    summary = read_tpx_summary('triplex_search.summary', args)    
    
    rna_dict = fasta2dict(args.rna_fn)
    dna_dict = fasta2dict(args.dna_fn)
    output(all_tpx, summary, rna_dict, dna_dict, args)
    if args.tpx_cov is not None:
        save_tpx_coverage(all_tpx, summary, rna_dict, dna_dict, args)
    if args.all_tpx is not None:
        save_all_tpx_info(all_tpx, args)

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

def output(all_tpx, summary, rna_dict, dna_dict, args):
    head_arr = ['RNA', 'DNA', 'RNA_len', 'DNA_len', 'Num_tpx', 'Max_score', 'Sum_score',  't_pot', 'log_tpot']
    print('\t'.join(head_arr))
    for tpx in sorted(summary.values(), reverse=True, key=lambda k: k['t_pot']):
        rna_seq = rna_dict.get(tpx['Sequence-ID'], '')
        dna_seq = dna_dict.get(tpx['Duplex-ID'],   '')
        sites = all_tpx[ tpx['key'] ]
        scores = [int(d['Score']) for d in sites]
        to_print = [
            tpx['Sequence-ID'], tpx['Duplex-ID'], str(len(rna_seq)), str(len(dna_seq)),
            str(len(sites)), str(max(scores)), str(sum(scores)), str(tpx['t_pot']),
            '%.4f' % tpx['log_tpot'], 
        ]
        print('\t'.join(to_print))
    
    # Output pairs without triplexes at the end
    for rna_name in rna_dict.keys():
        for dna_name in dna_dict.keys():
            key = rna_name + args.delim + dna_name
            if(key in summary.keys()): continue
            rna_seq  = rna_dict[rna_name]
            dna_seq  = dna_dict[dna_name]
            to_print = [rna_name, dna_name, str(len(rna_seq)), str(len(dna_seq)), '0','0','0','0','0']
            print('\t'.join(to_print))

def run_triplexator(args):
    # https://stackoverflow.com/a/18244485/310453
    cmd_str = 'triplexator ' + args.prm + ' -ss ' + args.rna_fn + ' -ds ' + args.dna_fn
    out_arr = subprocess.getoutput(cmd_str).splitlines()
    
    head_str = out_arr.pop(0)          # get the header row
    head_str = head_str.lstrip('# ')   # remove the comment line
    head_arr = head_str.split('\t')
    
    all_tpx = {}
    for line in out_arr:
        vals = line.split('\t')
        if( len(vals) != len(head_arr) ):
            warn('Something is wrong: len(vals) != len(head_arr)')
            continue
        tpx  = dict( zip(head_arr, vals) )
        key  = tpx['Sequence-ID'] + args.delim + tpx['Duplex-ID']
        if(key not in all_tpx):
            all_tpx[key] = []
        all_tpx[key].append(tpx)
    return all_tpx

def read_tpx_summary(fn, args):
    with open(fn) as f:
        head_str = f.readline()          # Skip the header: https://stackoverflow.com/a/4796785/310453
        head_str = head_str.lstrip('# ')   # remove the comment line
        head_arr = head_str.split('\t')
        
        summary = {}
        for line in f:
            vals = re.split('\s+', line)
            tpx  = dict( zip(head_arr, vals) )
            tpx['t_pot'] = float( tpx['Total (rel)'] )
            tpx['log_tpot'] = log10(tpx['t_pot'] * 10**6 + 1)

            key  = tpx['Sequence-ID'] + args.delim + tpx['Duplex-ID']
            tpx['key'] = key
            
            if(key in summary):
                warn("The key '%s' is duplicated" % key)
            else:
                summary[key] = tpx
    return summary

###
# Parse command line arguments: https://stackoverflow.com/a/30493366/310453
parser = argparse.ArgumentParser(description='Run the Triplexator and process its output')
parser.add_argument('rna_fn',   metavar='RNA.fna',    help='corresponds to the -ss Triplexator option')
parser.add_argument('dna_fn',   metavar='DNA.fna',    help='corresponds to the -ds Triplexator option')
parser.add_argument('--prm',    metavar='STR',        help='specify additional Triplexator params', default='')
parser.add_argument('--tpx_cov', metavar='PREFIX', help='compute triplex coverage info  for both RNA and DNA')
parser.add_argument('--all_tpx', metavar='FN.txt', help='save info about all the predicted triplexes')
parser.add_argument('--silent', action ='store_true', help='')
args = parser.parse_args()

###
#my $START_TIME = time;

args.rna_fn = os.path.abspath(args.rna_fn)
args.dna_fn = os.path.abspath(args.dna_fn)
args.delim = '_vs_'
args.init_cwd = os.getcwd()
args.tmp_dir = tempfile.mkdtemp(prefix='__TPX.', dir=os.getcwd())

os.chdir(args.tmp_dir)
run( args )
shutil.rmtree(args.tmp_dir)    # remove tmp folder

#warn "\nElapsed time: ".(time-$START_TIME)." sec\n" if !$SILENT;
###

