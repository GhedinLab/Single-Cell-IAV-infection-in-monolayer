import re
import argparse
import collections


parser = argparse.ArgumentParser('Extract all the gapped alignment (corresponding to DIs with one internal deletion) with unique UMIs, export them in sam files and count their numbers. Notice that alignments with two Ns will be exluced.')

parser.add_argument(
    '--gtf_file', required=True, help='gtf file')

parser.add_argument(
    '--input_sam_file', required=True, help='input sam file')

parser.add_argument(
    '--output_sam_file', required=True, help='output sam file')

parser.add_argument(
    '--counts_file', required=True, help='count file (e.g., "counts.txt")')

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
    'MIMDMNM',
    'MIMIMNM',
    'MNMIMDM',
    'MDMNMIMS',
    'SMIMNMDM',
    'MIMNMIM',
    'MDMNMIM',
    'MIMNMDM',
    'MDMNMDM'
}


def is_header(line):
    return line.startswith('@')


def is_in_range(pattern, cigar, split_line, gtf):
    segment = split_line[2]
    info = gtf[segment]
    firstbase = int(split_line[3])

    if pattern == 'MNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] - 1

    if pattern == 'SMNM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] - 1

    if pattern == 'MNMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] - 1

    if pattern == 'MNMDM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] - 1

    if pattern == 'SMNMS':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] - 1

    if pattern == 'MDMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] - 1

    if pattern == 'MNMIM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] - 1

    if pattern == 'MIMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] - 1

    if pattern == 'SMDMNM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] - 1

    if pattern == 'MNMDMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] - 1

    if pattern == 'SMNMDM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] - 1

    if pattern == 'SMNMIM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][5] - 1

    if pattern == 'SMIMNM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] - 1

    if pattern == 'MDMNMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] - 1

    if pattern == 'MIMNMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] - 1

    if pattern == 'MNMIMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] - 1

    if pattern == 'MDMIMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'SMDMNMS':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] - 1

    if pattern == 'SMNMIMS':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][5] - 1

    if pattern == 'SMIMNMS':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] - 1

    if pattern == 'MDMDMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MNMIMIM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] + cigar['len'][6] - 1

    if pattern == 'SMNMDMS':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] - 1

    if pattern == 'MNMDMDM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MNMDMIM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][6] - 1

    if pattern == 'MIMDMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MIMIMNM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][2] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MNMIMDM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MDMNMIMS':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][6] - 1

    if pattern == 'SMIMNMDM':
        lastbase = firstbase + cigar['len'][1] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] + cigar['len'][7] - 1

    if pattern == 'MIMNMIM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][6] - 1

    if pattern == 'MDMNMIM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][6] - 1

    if pattern == 'MIMNMDM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    if pattern == 'MDMNMDM':
        lastbase = firstbase + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5] + cigar['len'][6] - 1

    return (info['min'] <= firstbase) and (lastbase <= info['max'])


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


def main():
    gtf = {}
    with open(args.gtf_file) as f:
        for line in f:
            split_line = line.split('\t')
            gtf[split_line[0]] = {
                'min': int(split_line[3]),
                'max': int(split_line[4])
            }

    headers = []
    with open(args.input_sam_file) as f:
        segs = collections.defaultdict(list)
        for line in f:
            if not is_header(line):
                split_line = line.split('\t')
                segs[split_line[2]].append(split_line)
            else:
                headers.append(line)

    segs = collections.OrderedDict(sorted(segs.items(), key=lambda t: t[0]))
    filtered_segs = collections.OrderedDict()
    for seg, split_lines in segs.iteritems():
        umis = set()
        filtered_split_lines = []
        for split_line in split_lines:
            sposition1, cigar = extract_align_info(split_line)
            pattern = ''.join(cigar['cha'])
            if pattern in PATTERNS:
                umi = extract_umi(split_line)
                if umi is not None:
                    if is_in_range(pattern, cigar, split_line, gtf):
                        if umi not in umis:
                            umis.add(umi)
                            filtered_split_lines.append(split_line)
        filtered_segs[seg] = filtered_split_lines

    with open(args.output_sam_file, 'w') as f:
        f.writelines(headers)
        for seg, filtered_split_lines in filtered_segs.iteritems():
            f.writelines(map(lambda x: '\t'.join(x), filtered_split_lines))

    with open(args.counts_file, 'w') as f:
        for seg, filtered_split_lines in filtered_segs.iteritems():
            f.write('{}\t{}\n'.format(seg, len(filtered_split_lines)))


if __name__ == '__main__':
    main()
