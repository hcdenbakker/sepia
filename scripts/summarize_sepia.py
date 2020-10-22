#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  2 14:00:25 2019

@author: henkdenbakker
given a taxonomy file & a reads file:
    1. reconcile multiple hits with ancestral node
    2. remove ids below a certain threshold
    3. create krona file
"""
import sys
from sys import argv
from collections import defaultdict
from argparse import ArgumentParser
import numpy as np
import gzip
import io

def parse_args():
    "Parse the input arguments, use '-h' for help."
    parser = ArgumentParser(description='Reconcile - read classification processing tool using reconcilaiation of taxonomy')
    # inputs
    parser.add_argument(
        '-i','--input_file', type=str, required=True, default= 'None', metavar = 'input_file',
        help='Read classification file produced by ColorID')
    parser.add_argument(
        '-p','--prefix', type=str, required=True, default='None', metavar = 'prefix',
        help='Prefix for output files.')
    parser.add_argument(
        '-r','--ratio', type=str, required=False, default= 0.0, metavar = 'ratio',
        help='Ratio of hits to non-hits to use for filtering.')
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()

def reconcile(reads_file, out_prefix, threshold):
    #create dict with species:lineage
    kmer_ratios = defaultdict(list)
    kronafile = open(out_prefix + '.krona.txt', 'w')
    plus_file = open(out_prefix + '.plus.txt', 'w')
    summary = {}
    with open(reads_file, 'r') as reads:
        for line in reads:
            line_split = line.split('\t')
            classification = line_split[1]
            if int(line_split[3]) == 0:
                tax_lineage = line_split[1]
            else:
                if float(line_split[2])/float(line_split[3]) < threshold:
                    tax_lineage = 'below_threshold'
                else:    
                    tax_lineage = line_split[1]
            if tax_lineage in summary:
                summary[tax_lineage] += 1
            else:
                summary[tax_lineage] = 1
            if line_split[3] == '0':
                kmer_ratios[tax_lineage].append(0.0)
            else:
                try:
                    kmer_ratios[tax_lineage].append(float(line_split[2])/float(line_split[3]))
                except ZeroDivisionError:
                    kmer_ratios[tax_lineage].append(0.0)
    for f in summary:    
        kronafile.write(str(float(summary[f])) + '\t' + '\t'.join(f.split(';' )) + '\n')
        average = round(sum(kmer_ratios[f])/len(kmer_ratios[f]), 3)
        array = np.histogram(kmer_ratios[f], bins = 50)
        tensor = list(array[0]/sum(array[0]))
        tensor_string = [str(a) for a in tensor]
        plus_file.write(f + '\t' + str(summary[f]) + '\t' + str(average) + '\t' + '\t'.join(tensor_string) + '\n')
    plus_file.close()
    kronafile.close()

def reconcile_gz(reads_file, out_prefix, threshold):
    #create dict with species:lineage
    kmer_ratios = defaultdict(list)
    kronafile = open(out_prefix + '.krona.txt', 'w')
    plus_file = open(out_prefix + '.plus.txt', 'w')
    summary = {}
    reads = io.BufferedReader(gzip.open(reads_file))
    for line in reads:
        line=line.decode()
        line_split = line.split('\t')
        classification = line_split[1]
        if int(line_split[3]) == 0:
            tax_lineage = line_split[1]
        else:
            if float(line_split[2])/float(line_split[3]) < threshold:
                tax_lineage = 'below_threshold'
            else:
                tax_lineage = line_split[1]
        if tax_lineage in summary:
            summary[tax_lineage] += 1
        else:
            summary[tax_lineage] = 1
        if line_split[3] == '0':
            kmer_ratios[tax_lineage].append(0.0)
        else:
            try:
                kmer_ratios[tax_lineage].append(float(line_split[2])/float(line_split[3]))
            except ZeroDivisionError:
                kmer_ratios[tax_lineage].append(0.0)
    for f in summary:
        kronafile.write(str(float(summary[f])) + '\t' + '\t'.join(f.split(';' )) + '\n')
        average = round(sum(kmer_ratios[f])/len(kmer_ratios[f]), 3)
        array = np.histogram(kmer_ratios[f], bins = 50)
        tensor = list(array[0]/sum(array[0]))
        tensor_string = [str(a) for a in tensor]
        plus_file.write(f + '\t' + str(summary[f]) + '\t' + str(average) + '\t' + '\t'.join(tensor_string) + '\n')
    plus_file.close()
    kronafile.close()

def main():
    
    args = parse_args()
    reads_file = args.input_file
    out_prefix = args.prefix
    ratio = args.ratio
    if reads_file.endswith('gz'):
        reconcile_gz(reads_file, out_prefix,float(ratio))
    else:
        reconcile(reads_file, out_prefix,float(ratio))


if __name__ == '__main__':
    main()
