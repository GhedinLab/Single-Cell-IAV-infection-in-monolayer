import argparse
import os
import re

parser = argparse.ArgumentParser('Extract the patterns in cigar string for each IAV sam file in a directory.')

parser.add_argument(
    '-indir', '--input_dir', required=True, help='input directory containing the IAV sam file for each cell')
parser.add_argument(
    '-o', '--output', required=True, help='the file name for the output')

args = parser.parse_args()


DIGIT = re.compile('\d+')

patterns = {}

for fn in os.listdir(args.input_dir):
    with open(fn) as f:
        for line in f:
            if not line.startswith('@'):
                pattern = DIGIT.sub('^', line.split()[5])
                if pattern in patterns:
                    patterns[pattern] += 1
                else:
                    patterns[pattern] = 1

with open(args.output, 'a') as f:
    f.write('\n'.join(['%s\t%i' % (k, v) for k, v in patterns.iteritems()]))
