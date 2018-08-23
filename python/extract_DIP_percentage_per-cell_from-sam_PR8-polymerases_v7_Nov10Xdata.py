import os
import re
import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser('Extract the junction coordinates of defective interferring particles for each polymerase segment in each cell from sam file. Note: the given boundaries are corresponding to the last base at the 3prime end of the first proportion of the gapped reads and the 5prime end of the last proportion of the gapped reads.')

parser.add_argument(
    'files', nargs='+', help='input sam file for each cell')
parser.add_argument(
    '-m', '--min_length', required=True, type=int, help='the minimum length of parts of reads mapped to 5 and 3 end separately')
parser.add_argument(
    '-sl', '--skip_length', required=True, type=int, help='the minimum length of the skip region in cigar')
parser.add_argument(
    '-p', '--percentile', required=True, type=int, help='the percentile of boundary range to cover in determining DIP (e.g. 95 or 100)')
parser.add_argument(
    '-od', '--output_dir', required=True, help='the directory for all the output files')

args = parser.parse_args()

uppercase = re.compile(r'([A-Z])')


def is_header(line):
    return line.startswith('@')


# def iav_filter(split_line, info_table):
#    return split_line[2] in info_table.segment.values


def in_range(value, range_5, range_3):
    return range_5 <= value <= range_3


# def intersect(value, length, peak, peak_range):
#    return (in_range(value, peak - peak_range, peak + peak_range) or
#            in_range(value + length, peak - peak_range, peak + peak_range))


def extract_align_info(split_line):
    sposition1 = int(split_line[3])
    cigar = {'cha': [], 'len': []}
    last = None
    for i, s in enumerate(re.split(uppercase, split_line[5])):
        if s:
            if i % 2 == 0:
                last = int(s)
            else:
                cigar['cha'].append(s)
                cigar['len'].append(last)
    return sposition1, cigar


def record_dip_boundary(file_, min_length, skip_length):
    boundary_all = {'AF389115.1': [[], []], 'AF389116.1': [[], []], 'AF389117.1': [[], []]}
    with open(file_) as f:
        for line in f:
            if not is_header(line):
                split_line = line.split('\t')
#                if iav_filter(split_line, info_table):
                sposition1, cigar = extract_align_info(split_line)
                segment = split_line[2]
                if segment in boundary_all.keys():
                    # For each segment ('key'), the 'value' is a nested list including two list, boundary5 and boundary3
                    pattern = ''.join(cigar['cha'])
                    if pattern == 'MNM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if cigar['len'][0] >= min_length and cigar['len'][-1] >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][1] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][0] + cigar['len'][1])
                    if pattern == 'SMNM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if cigar['len'][1] >= min_length and cigar['len'][-1] >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][2] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][1] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][1] + cigar['len'][2])
                    if pattern == 'MNMS':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if cigar['len'][0] >= min_length and cigar['len'][2] >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][1] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][0] + cigar['len'][1])
                    if pattern == 'MNMDM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if cigar['len'][0] >= min_length and (cigar['len'][2] + cigar['len'][-1]) >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][1] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][0] + cigar['len'][1])
                    if pattern == 'SMNMS':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if cigar['len'][1] >= min_length and cigar['len'][3] >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][2] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][1] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][1] + cigar['len'][2])
                    if pattern == 'MDMNM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if (cigar['len'][0] + cigar['len'][2]) >= min_length and cigar['len'][-1] >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][3] >= skip_length:
                                boundary_all[segment][0].append(
                                    sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] - 1)
                                boundary_all[segment][1].append(
                                    sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3])
                    if pattern == 'MNMIM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if cigar['len'][0] >= min_length and (cigar['len'][2] + cigar['len'][-1]) >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][1] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][0] + cigar['len'][1])
                    if pattern == 'MIMNM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if (cigar['len'][0] + cigar['len'][2]) >= min_length and cigar['len'][-1] >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][3] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] + cigar['len'][2] - 1)
                                boundary_all[segment][1].append(
                                    sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][3])
                    if pattern == 'SMDMNM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if (cigar['len'][1] + cigar['len'][3]) >= min_length and cigar['len'][-1] >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][4] >= skip_length:
                                boundary_all[segment][0].append(
                                    sposition1 + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] - 1)
                                boundary_all[segment][1].append(
                                    sposition1 + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4])
                    if pattern == 'MNMDMS':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if cigar['len'][0] >= min_length and (cigar['len'][2] + cigar['len'][4]) >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][1] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][0] + cigar['len'][1])
                    if pattern == 'SMNMDM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if cigar['len'][1] >= min_length and (cigar['len'][3] + cigar['len'][-1]) >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][2] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][1] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][1] + cigar['len'][2])
                    if pattern == 'SMNMIM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if cigar['len'][1] >= min_length and (cigar['len'][3] + cigar['len'][-1]) >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][2] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][1] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][1] + cigar['len'][2])
                    if pattern == 'SMIMNM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if (cigar['len'][1] + cigar['len'][3]) >= min_length and cigar['len'][-1] >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][4] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][1] + cigar['len'][3] - 1)
                                boundary_all[segment][1].append(
                                    sposition1 + cigar['len'][1] + cigar['len'][3] + cigar['len'][4])
                    if pattern == 'MDMNMS':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if (cigar['len'][0] + cigar['len'][2]) >= min_length and cigar['len'][4] >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][3] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] + cigar['len'][1] +
                                                                cigar['len'][2] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][0] + cigar['len'][1] +
                                                                cigar['len'][2] + cigar['len'][3])
                    if pattern == 'MIMNMS':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if (cigar['len'][0] + cigar['len'][2]) >= min_length and cigar['len'][4] >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][3] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] + cigar['len'][2] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][3])
                    if pattern == 'MNMIMS':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if cigar['len'][0] >= min_length and (cigar['len'][2] + cigar['len'][4]) >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][1] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][0] + cigar['len'][1])
                    if pattern == 'MDMIMNM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if (cigar['len'][0] + cigar['len'][2] + cigar['len'][4]) >= min_length and cigar['len'][-1] >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][5] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] + cigar['len'][5])
                    if pattern == 'SMDMNMS':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if (cigar['len'][1] + cigar['len'][3]) >= min_length and cigar['len'][5] >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][4] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4])
                    if pattern == 'SMNMIMS':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if cigar['len'][1] >= min_length and (cigar['len'][3] + cigar['len'][5]) >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][2] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][1] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][1] + cigar['len'][2])
                    if pattern == 'SMIMNMS':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if (cigar['len'][1] + cigar['len'][3]) >= min_length and cigar['len'][5] >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][4] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][1] + cigar['len'][3] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][1] + cigar['len'][3] + cigar['len'][4])
                    if pattern == 'MDMDMNM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if (cigar['len'][0] + cigar['len'][2] + cigar['len'][4]) >= min_length and cigar['len'][-1] >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][5] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5])
                    if pattern == 'MNMIMIM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if cigar['len'][0] >= min_length and (cigar['len'][2] + cigar['len'][4] + cigar['len'][-1]) >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][1] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][0] + cigar['len'][1])
                    if pattern == 'SMNMDMS':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if cigar['len'][1] >= min_length and (cigar['len'][3] + cigar['len'][5]) >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][2] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][1] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][1] + cigar['len'][2])
                    if pattern == 'MNMDMDM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if cigar['len'][0] >= min_length and (cigar['len'][2] + cigar['len'][4] + cigar['len'][-1]) >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][1] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][0] + cigar['len'][1])
                    if pattern == 'MNMDMIM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if cigar['len'][0] >= min_length and (cigar['len'][2] + cigar['len'][4] + cigar['len'][-1]) >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][1] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][0] + cigar['len'][1])
                    if pattern == 'MIMDMNM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if (cigar['len'][0] + cigar['len'][2] + cigar['len'][4]) >= min_length and cigar['len'][-1] >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][5] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5])
                    if pattern == 'MIMIMNM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if (cigar['len'][0] + cigar['len'][2] + cigar['len'][4]) >= min_length and cigar['len'][-1] >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][5] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][4] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][4] + cigar['len'][5])
                    if pattern == 'MNMIMDM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if cigar['len'][0] >= min_length and (cigar['len'][2] + cigar['len'][4] + cigar['len'][-1]) >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][1] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][0] + cigar['len'][1])
                    if pattern == 'MDMNMIMS':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if (cigar['len'][0] + cigar['len'][2]) >= min_length and (cigar['len'][4] + cigar['len'][6]) >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][3] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3])
                    if pattern == 'SMIMNMDM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if (cigar['len'][1] + cigar['len'][3]) >= min_length and (cigar['len'][5] + cigar['len'][-1]) >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][4] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][1] + cigar['len'][3] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][1] + cigar['len'][3] + cigar['len'][4])
                    if pattern == 'MIMNMIM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if (cigar['len'][0] + cigar['len'][2]) >= min_length and (cigar['len'][4] + cigar['len'][-1]) >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][3] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] + cigar['len'][2] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][3])
                    if pattern == 'MDMNMIM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if (cigar['len'][0] + cigar['len'][2]) >= min_length and (cigar['len'][4] + cigar['len'][-1]) >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][3] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3])
                    if pattern == 'MIMNMDM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if (cigar['len'][0] + cigar['len'][2]) >= min_length and (cigar['len'][4] + cigar['len'][-1]) >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][3] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] + cigar['len'][2] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][3])
                    if pattern == 'MDMNMDM':
                        # filter reads based on the length of two parts of soft-clipped reads
                        if (cigar['len'][0] + cigar['len'][2]) >= min_length and (cigar['len'][4] + cigar['len'][-1]) >= min_length:
                            # filter soft-clipped reads based on the length of the skipped region
                            if cigar['len'][3] >= skip_length:
                                boundary_all[segment][0].append(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] - 1)
                                boundary_all[segment][1].append(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3])
    return {k: v for k, v in boundary_all.iteritems() if len(v[0]) > 0}


def generate_info_table(boundary_all, percentile):
    percentile1 = float(100 - percentile) / 2
    percentile2 = float(100 - percentile) / 2 + percentile

    info_table = pd.DataFrame({'segment': boundary_all.keys()})
#    info_table = info_table[info_table.segment.isin(boundary_all.keys())].copy()

    info_table['boundary5_range5'] = ''
    info_table['boundary5_range3'] = ''
    info_table['boundary3_range5'] = ''
    info_table['boundary3_range3'] = ''
    for key in boundary_all:
        boundary5_all = boundary_all[key][0]
        boundary3_all = boundary_all[key][1]
        boundary5_all = np.array(boundary5_all)
        boundary3_all = np.array(boundary3_all)
        boundary5_range_5 = np.percentile(boundary5_all, percentile1)
        boundary5_range_3 = np.percentile(boundary5_all, percentile2)
        boundary3_range_5 = np.percentile(boundary3_all, percentile1)
        boundary3_range_3 = np.percentile(boundary3_all, percentile2)
        info_table.ix[info_table.segment == key, 'boundary5_range5'] = boundary5_range_5
        info_table.ix[info_table.segment == key, 'boundary5_range3'] = boundary5_range_3
        info_table.ix[info_table.segment == key, 'boundary3_range5'] = boundary3_range_5
        info_table.ix[info_table.segment == key, 'boundary3_range3'] = boundary3_range_3
    return info_table


def dip_filter_record(split_line, sposition1, cigar, min_length, skip_length, info_table):
    pattern = ''.join(cigar['cha'])
    if pattern == 'MNM':
        # filter reads based on the length of two parts of soft-clipped reads
        if cigar['len'][0] >= min_length and cigar['len'][-1] >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][1] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][1],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][1]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'SMNM':
        # filter reads based on the length of two parts of soft-clipped reads
        if cigar['len'][1] >= min_length and cigar['len'][-1] >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][2] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][1] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][1] + cigar['len'][2],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][1] - 1
                        boundary3 = sposition1 + cigar['len'][1] + cigar['len'][2]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MNMS':
        # filter reads based on the length of two parts of soft-clipped reads
        if cigar['len'][0] >= min_length and cigar['len'][2] >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][1] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][1],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][1]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MNMDM':
        # filter reads based on the length of two parts of soft-clipped reads
        if cigar['len'][0] >= min_length and (cigar['len'][2] + cigar['len'][-1]) >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][1] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][1],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][1]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'SMNMS':
        # filter reads based on the length of two parts of soft-clipped reads
        if cigar['len'][1] >= min_length and cigar['len'][3] >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][2] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][1] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][1] + cigar['len'][2],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][1] - 1
                        boundary3 = sposition1 + cigar['len'][1] + cigar['len'][2]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MDMNM':
        # filter reads based on the length of two parts of soft-clipped reads
        if (cigar['len'][0] + cigar['len'][2]) >= min_length and cigar['len'][-1] >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][3] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MNMIM':
        # filter reads based on the length of two parts of soft-clipped reads
        if cigar['len'][0] >= min_length and (cigar['len'][2] + cigar['len'][-1]) >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][1] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][1],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][1]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MIMNM':
        # filter reads based on the length of two parts of soft-clipped reads
        if (cigar['len'][0] + cigar['len'][2]) >= min_length and cigar['len'][-1] >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][3] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] + cigar['len'][2] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][3],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] + cigar['len'][2] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][3]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'SMDMNM':
        # filter reads based on the length of two parts of soft-clipped reads
        if (cigar['len'][1] + cigar['len'][3]) >= min_length and cigar['len'][-1] >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][4] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] - 1
                        boundary3 = sposition1 + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MNMDMS':
        # filter reads based on the length of two parts of soft-clipped reads
        if cigar['len'][0] >= min_length and (cigar['len'][2] + cigar['len'][4]) >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][1] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][1],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][1]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'SMNMDM':
        # filter reads based on the length of two parts of soft-clipped reads
        if cigar['len'][1] >= min_length and (cigar['len'][3] + cigar['len'][-1]) >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][2] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][1] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][1] + cigar['len'][2],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][1] - 1
                        boundary3 = sposition1 + cigar['len'][1] + cigar['len'][2]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'SMNMIM':
        # filter reads based on the length of two parts of soft-clipped reads
        if cigar['len'][1] >= min_length and (cigar['len'][3] + cigar['len'][-1]) >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][2] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][1] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][1] + cigar['len'][2],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][1] - 1
                        boundary3 = sposition1 + cigar['len'][1] + cigar['len'][2]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'SMIMNM':
        # filter reads based on the length of two parts of soft-clipped reads
        if (cigar['len'][1] + cigar['len'][3]) >= min_length and cigar['len'][-1] >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][4] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][1] + cigar['len'][3] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][1] + cigar['len'][3] + cigar['len'][4],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][1] + cigar['len'][3] - 1
                        boundary3 = sposition1 + cigar['len'][1] + cigar['len'][3] + cigar['len'][4]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MDMNMS':
        # filter reads based on the length of two parts of soft-clipped reads
        if (cigar['len'][0] + cigar['len'][2]) >= min_length and cigar['len'][4] >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][3] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MIMNMS':
        # filter reads based on the length of two parts of soft-clipped reads
        if (cigar['len'][0] + cigar['len'][2]) >= min_length and cigar['len'][4] >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][3] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] + cigar['len'][2] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][3],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] + cigar['len'][2] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][3]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MNMIMS':
        # filter reads based on the length of two parts of soft-clipped reads
        if cigar['len'][0] >= min_length and (cigar['len'][2] + cigar['len'][4]) >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][1] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][1],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][1]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MDMIMNM':
        # filter reads based on the length of two parts of soft-clipped reads
        if (cigar['len'][0] + cigar['len'][2] + cigar['len'][4]) >= min_length and cigar['len'][-1] >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][5] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] + cigar['len'][5],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][4] + cigar['len'][5]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'SMDMNMS':
        # filter reads based on the length of two parts of soft-clipped reads
        if (cigar['len'][1] + cigar['len'][3]) >= min_length and cigar['len'][5] >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][4] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] - 1
                        boundary3 = sposition1 + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'SMNMIMS':
        # filter reads based on the length of two parts of soft-clipped reads
        if cigar['len'][1] >= min_length and (cigar['len'][3] + cigar['len'][5]) >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][2] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][1] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][1] + cigar['len'][2],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][1] - 1
                        boundary3 = sposition1 + cigar['len'][1] + cigar['len'][2]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'SMIMNMS':
        # filter reads based on the length of two parts of soft-clipped reads
        if (cigar['len'][1] + cigar['len'][3]) >= min_length and cigar['len'][5] >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][4] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][1] + cigar['len'][3] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][1] + cigar['len'][3] + cigar['len'][4],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][1] + cigar['len'][3] - 1
                        boundary3 = sposition1 + cigar['len'][1] + cigar['len'][3] + cigar['len'][4]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MDMDMNM':
        # filter reads based on the length of two parts of soft-clipped reads
        if (cigar['len'][0] + cigar['len'][2] + cigar['len'][4]) >= min_length and cigar['len'][-1] >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][5] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MNMIMIM':
        # filter reads based on the length of two parts of soft-clipped reads
        if cigar['len'][0] >= min_length and (cigar['len'][2] + cigar['len'][4] + cigar['len'][-1]) >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][1] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][1],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][1]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'SMNMDMS':
        # filter reads based on the length of two parts of soft-clipped reads
        if cigar['len'][1] >= min_length and (cigar['len'][3] + cigar['len'][5]) >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][2] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][1] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][1] + cigar['len'][2],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][1] - 1
                        boundary3 = sposition1 + cigar['len'][1] + cigar['len'][2]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MNMDMDM':
        # filter reads based on the length of two parts of soft-clipped reads
        if cigar['len'][0] >= min_length and (cigar['len'][2] + cigar['len'][4] + cigar['len'][-1]) >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][1] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][1],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][1]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MNMDMIM':
        # filter reads based on the length of two parts of soft-clipped reads
        if cigar['len'][0] >= min_length and (cigar['len'][2] + cigar['len'][4] + cigar['len'][-1]) >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][1] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][1],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][1]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MIMDMNM':
        # filter reads based on the length of two parts of soft-clipped reads
        if (cigar['len'][0] + cigar['len'][2] + cigar['len'][4]) >= min_length and cigar['len'][-1] >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][5] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][3] + cigar['len'][4] + cigar['len'][5]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MIMIMNM':
        # filter reads based on the length of two parts of soft-clipped reads
        if (cigar['len'][0] + cigar['len'][2] + cigar['len'][4]) >= min_length and cigar['len'][-1] >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][5] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][4] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][4] + cigar['len'][5],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][4] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][4] + cigar['len'][5]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MNMIMDM':
        # filter reads based on the length of two parts of soft-clipped reads
        if cigar['len'][0] >= min_length and (cigar['len'][2] + cigar['len'][4] + cigar['len'][-1]) >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][1] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][1],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][1]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MDMNMIMS':
        # filter reads based on the length of two parts of soft-clipped reads
        if (cigar['len'][0] + cigar['len'][2]) >= min_length and (cigar['len'][4] + cigar['len'][6]) >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][3] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'SMIMNMDM':
        # filter reads based on the length of two parts of soft-clipped reads
        if (cigar['len'][1] + cigar['len'][3]) >= min_length and (cigar['len'][5] + cigar['len'][-1]) >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][4] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][1] + cigar['len'][3] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][1] + cigar['len'][3] + cigar['len'][4],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][1] + cigar['len'][3] - 1
                        boundary3 = sposition1 + cigar['len'][1] + cigar['len'][3] + cigar['len'][4]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MIMNMIM':
        # filter reads based on the length of two parts of soft-clipped reads
        if (cigar['len'][0] + cigar['len'][2]) >= min_length and (cigar['len'][4] + cigar['len'][-1]) >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][3] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] + cigar['len'][2] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][3],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] + cigar['len'][2] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][3]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MDMNMIM':
        # filter reads based on the length of two parts of soft-clipped reads
        if (cigar['len'][0] + cigar['len'][2]) >= min_length and (cigar['len'][4] + cigar['len'][-1]) >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][3] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MIMNMDM':
        # filter reads based on the length of two parts of soft-clipped reads
        if (cigar['len'][0] + cigar['len'][2]) >= min_length and (cigar['len'][4] + cigar['len'][-1]) >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][3] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] + cigar['len'][2] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][3],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] + cigar['len'][2] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][2] + cigar['len'][3]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary
    if pattern == 'MDMNMDM':
        # filter reads based on the length of two parts of soft-clipped reads
        if (cigar['len'][0] + cigar['len'][2]) >= min_length and (cigar['len'][4] + cigar['len'][-1]) >= min_length:
            # filter soft-clipped reads based on the length of the skipped region
            if cigar['len'][3] >= skip_length:
                segment = split_line[2]
                info = info_table[info_table['segment'] == segment].iloc[0].to_dict()
                if in_range(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] - 1,
                            int(info['boundary5_range5']),
                            int(info['boundary5_range3'])):
                    if in_range(sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3],
                                int(info['boundary3_range5']),
                                int(info['boundary3_range3'])):
                        boundary5 = sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] - 1
                        boundary3 = sposition1 + cigar['len'][0] + cigar['len'][1] + cigar['len'][2] + cigar['len'][3]
                        dip_filtered_boundary = (segment, boundary5, boundary3)
                        return dip_filtered_boundary


# def seg_sum(file_dir, cid, info_table, peak_range):
#    bam_file = pysam.AlignmentFile(os.path.join(file_dir, cid + '.sorted.bam'), 'rb')
#    seg_3prime_sum = {}
#    for segment in info_table.segment:
#        start = int(info_table.ix[info_table.segment == segment, cid].iloc[0]) - peak_range
#        stop = int(info_table.ix[info_table.segment == segment, cid].iloc[0]) + peak_range
#        if start > 0:
#            individual_seg_sum = list(bam_file.fetch(segment, start, stop, until_eof=True))
#            seg_3prime_sum.update({segment: len(individual_seg_sum)})
#        else:
#            seg_3prime_sum.update({segment: 0})
#    return seg_3prime_sum


def main():
#    info_table = pd.read_csv(args.info_table)
#    read_sum_master = {}
    segment_boundary = {'AF389115.1': {}, 'AF389116.1': {}, 'AF389117.1': {}}
    for file_ in args.files:
        print 'Processing file: ' + file_
        cid = os.path.basename(file_).split('.')[0]
#        sums_counts = seg_sum(args.file_dir, cid, info_table, args.peak_range)
#        read_sum_master.update({cid: sums_counts})
        boundary_all = record_dip_boundary(file_, args.min_length, args.skip_length)
        info_table_new = generate_info_table(boundary_all, args.percentile)
        for segment in segment_boundary:
            segment_boundary[segment].update({cid: []})
        if not info_table_new.empty:
            with open(file_) as f:
                for line in f:
                    if not is_header(line):
                        split_line = line.split('\t')
#                        if iav_filter(split_line, info_table):
                        if split_line[2] in segment_boundary.keys(): # Only do for the polymerase segments
                            sposition1, cigar = extract_align_info(split_line)
                            dip_filtered_boundary = dip_filter_record(
                                split_line, sposition1, cigar, args.min_length, args.skip_length, info_table_new)
                            if dip_filtered_boundary is not None:
                                segment_boundary[dip_filtered_boundary[0]][cid].append([dip_filtered_boundary[1],
                                                                                        dip_filtered_boundary[2]])
        else:
            for segment in info_table_new.segment.values:
                segment_boundary[segment].update({cid: [[np.nan, np.nan]]})
#    read_sum_master_df = pd.DataFrame(read_sum_master)
#    read_sum_master_df.to_csv(os.path.join(args.output_dir, 'reads_mapped-to-3p-peak-each-segment_per-sample.csv'))
    boundary5_all = {}
    boundary3_all = {}
    for seg in segment_boundary:
        for cid in segment_boundary[seg]:
            boundary5_all[cid] = pd.Series([iterm[0] for iterm in segment_boundary[seg][cid]])
            boundary3_all[cid] = pd.Series([iterm[1] for iterm in segment_boundary[seg][cid]])
        boundary5_all_df = pd.DataFrame(boundary5_all).transpose().fillna('na')
        boundary3_all_df = pd.DataFrame(boundary3_all).transpose().fillna('na')
        boundary5_all_df.to_csv(os.path.join(args.output_dir, '_'.join((seg, 'boundary_5prime-end.csv'))))
        boundary3_all_df.to_csv(os.path.join(args.output_dir, '_'.join((seg, 'boundary_3prime-end.csv'))))


if __name__ == '__main__':
    main()
