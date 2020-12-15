#!/usr/bin/env python3

import os
import sys
import argparse

parser = argparse.ArgumentParser(
    description='Summarize hmmscan domtblout output file')
parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                    default=sys.stdin,
                    help='hmmscan domtblout file')
parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                    default=sys.stdout,
                    help='summarized output')
parser.add_argument('-e', '--Evalue_threshold', type=float,
                    default=0.01,
                    help='Evalue threshold')
parser.add_argument('-i', '--iEvalue_threshold', type=float,
                    default=0.01,
                    help='iEvalue threshold')
parser.add_argument('-c', '--coverage_threshold', type=float,
                    default=0.6,
                    help='coverage threshold')
parser.add_argument('-sc', '--soft_coverage_threshold', type=float,
                    help='soft coverage threshold')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
args = parser.parse_args()

if (args.soft_coverage_threshold == None) or (args.soft_coverage_threshold > args.coverage_threshold):
    args.soft_coverage_threshold = args.coverage_threshold

header = '\t'.join(('Query name', 'Relaxase MOB family', 'Profile HMM', 'Coverage', 'Start', 'End', 'Evalue', 'i-Evalue'))
# args.outfile.write(header + '\n')

hit_found = False
partial_found = False
prev_line_query_name = ''
for line in args.infile:
    if line[0] == '#':
        continue
    items = line.strip().split()
    if len(items) != 23:
        print(os.path.basename(sys.argv[0]) + ':',
              'error: bad domtblout format', line, file=sys.stderr)
        sys.exit(1)

    line_target_family = items[0].strip().split('_')[1][:4]
    line_target_name = items[0]
    line_query_name = items[3]
    line_evalue = items[6]
    line_ievalue = items[12]
    line_tlen = int(items[2])
    line_hmm_from = int(items[15])
    line_hmm_to = int(items[16])
    line_hlen = line_hmm_to - line_hmm_from + 1
    line_coverage = line_hlen / line_tlen
    line_ali_from = items[17]
    line_ali_to = items[18]

    if prev_line_query_name != line_query_name:
        if hit_found or partial_found:
            output = '{0}\t{1}\t{2}\t{3:.2f}\t{4}\t{5}\t{6}\t{7}'.format(hit_query_name, hit_target_family, hit_target_name, hit_coverage, hit_ali_from, hit_ali_to, hit_evalue, hit_ievalue)
            args.outfile.write(output + '\n')
            hit_found = False
            partial_found = False
        prev_line_query_name = line_query_name

    if (float(line_ievalue) <= args.iEvalue_threshold) and (float(line_evalue) <= args.Evalue_threshold):
        if (line_coverage >= args.coverage_threshold):
            if hit_found:
                if float(line_ievalue) < float(hit_ievalue):
                    hit_target_family = line_target_family
                    hit_target_name = line_target_name
                    hit_query_name = line_query_name
                    hit_evalue = line_evalue
                    hit_ievalue = line_ievalue
                    hit_coverage = line_coverage
                    hit_ali_from = line_ali_from
                    hit_ali_to = line_ali_to
            else:
                hit_found = True
                partial_found = False

                hit_target_family = line_target_family
                hit_target_name = line_target_name
                hit_query_name = line_query_name
                hit_evalue = line_evalue
                hit_ievalue = line_ievalue
                hit_coverage = line_coverage
                hit_ali_from = line_ali_from
                hit_ali_to = line_ali_to
        elif (line_coverage >= args.soft_coverage_threshold) and not hit_found:
            if partial_found:
                if float(line_ievalue) < float(hit_ievalue):
                    hit_target_family = line_target_family
                    hit_target_name = line_target_name
                    hit_query_name = line_query_name
                    hit_evalue = line_evalue
                    hit_ievalue = line_ievalue
                    hit_coverage = line_coverage
                    hit_ali_from = line_ali_from
                    hit_ali_to = line_ali_to
            else:
                partial_found = True

                hit_target_family = line_target_family
                hit_target_name = line_target_name
                hit_query_name = line_query_name
                hit_evalue = line_evalue
                hit_ievalue = line_ievalue
                hit_coverage = line_coverage
                hit_ali_from = line_ali_from
                hit_ali_to = line_ali_to

if hit_found or partial_found:
    output = '{0}\t{1}\t{2}\t{3:.2f}\t{4}\t{5}\t{6}\t{7}'.format(hit_query_name, hit_target_family, hit_target_name, hit_coverage, hit_ali_from, hit_ali_to, hit_evalue, hit_ievalue)
    args.outfile.write(output + '\n')
