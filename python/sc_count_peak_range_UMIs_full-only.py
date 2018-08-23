import os
import re
import argparse
import pysam
import pandas as pd


parser = argparse.ArgumentParser('Count the number of unique UMIs including both from gapped and normal reads that partially or completely fall into the peak range. Notice that in pysam, the positions are 0-based.')

parser.add_argument(
    'files', nargs='+', help='input sam file for each cell')
parser.add_argument(
    '-bd', '--file_dir', required=True, help='the directory containing corresponding bam and bai files for individual cells')
parser.add_argument(
    '-r', '--peak_range', required=True, type=int, help='the one-side range of peak region to cover in determining the reads mapped to the segment (e.g. 50, which means +-50 from the peak)')
parser.add_argument(
    '-t', '--info_table', required=True, help='the csv file containing the information for the peak (columns: segment, cell1_peak3, cell2_peak3, .... )')
parser.add_argument(
    '-od', '--output_dir', required=True, help='the directory for all the output files')

args = parser.parse_args()

UPPERCASE = re.compile(r'([A-Z])')

PATTERNS = {
    'MNM',
    'SMNM',
    'MNMS',
    'MNMDM',
    'SMNMS',
    'MDMNM',
    'MNMIM',
    'MIMNM',
    'MNMNM',
    'SMDMNM',
    'MNMDMS',
    'SMNMDM',
    'SMNMIM',
    'SMIMNM',
    'MDMNMS',
    'MIMNMS',
    'MNMIMS',
    'MDMIMNM',
    'SMDMNMS',
    'SMNMIMS',
    'SMIMNMS',
    'MDMDMNM',
    'MNMIMIM',
    'SMNMDMS',
    'MNMDMDM',
    'MNMDMIM',
    'MNMNMS',
    'MIMDMNM',
    'MIMIMNM',
    'MNMIMDM',
    'MDMNMIMS',
    'SMIMNMDM',
    'MIMNMIM',
    'MIMNMNM',
    'MDMNMIM',
    'MIMNMDM',
    'SMNMNM',
    'MDMNMDM',
    'MNMIMNM',
    'MNMDMNM',
    'MNMNMIM',
    'MNMNMDM'
}

MAPPINGS = list('MIDNSHP=XB')


def is_header(line):
    return line.startswith('@')


def extract_align_info(split_line):
    sposition1 = int(split_line[3])
    cigar = {'cha': [], 'len': []}
    last = None
    for i, s in enumerate(re.split(UPPERCASE, split_line[5])):
        if s:
            if i % 2 == 0:
                last = int(s)
            else:
                cigar['cha'].append(s)
                cigar['len'].append(last)
    return sposition1, cigar


def extract_umi(split_line):
    for s in split_line:
        if s.startswith('UB:Z:'):
            return s[5:]


def get_gapped_umis(file_path):
    umis = set()
    with open(file_path) as f:
        for line in f:
            if not is_header(line):
                split_line = line.split('\t')
                _, cigar = extract_align_info(split_line)
                pattern = ''.join(cigar['cha'])
                if pattern in PATTERNS:
                    umi = extract_umi(split_line)
                    if umi is not None:
                        umis.add(umi)
    return umis


def seg_sum(file_path, cid, info_table, peak_range, gapped_umis):
    bam_file = pysam.AlignmentFile(file_path, 'rb')
    seg_3prime_sum = {}
    for segment in info_table.segment:
        start = int(info_table.ix[info_table.segment == segment, cid].iloc[0]) - peak_range - 1
        stop = int(info_table.ix[info_table.segment == segment, cid].iloc[0]) + peak_range - 1
        umis = set()
        if start > 0:
            for record in bam_file.fetch(segment, start, stop, until_eof=True):
                pattern = ''.join([MAPPINGS[cigartuple[0]] for cigartuple in record.cigartuples])
                if pattern not in PATTERNS:
                    try:
                        umi = record.get_tag('UB')
                        if umi not in gapped_umis:
                            umis.add(umi)
                    except KeyError:
                        pass
            seg_3prime_sum.update({segment: len(umis)})
        else:
            seg_3prime_sum.update({segment: 0})
    return seg_3prime_sum


def main():
    info_table = pd.read_csv(args.info_table)
    read_sum_master = {}
    for file_ in args.files:
        print 'Processing file: ' + file_
        cid = file_.split('.')[0]
        file_path = os.path.join(args.file_dir, cid + '.sorted.bam')
        gapped_umis = get_gapped_umis(file_)
        sums_counts = seg_sum(file_path, cid, info_table, args.peak_range, gapped_umis)
        read_sum_master.update({cid: sums_counts})
    read_sum_master_df = pd.DataFrame(read_sum_master)
    read_sum_master_df.to_csv(os.path.join(args.output_dir, 'reads_mapped-to-3p-peak-each-segment_per-sample_full-length_UMIcounts.csv'))


if __name__ == '__main__':
    main()
