import os
import argparse
import re

parser = argparse.ArgumentParser('Separate sam alignment for individual cells from cell ranger output.')

parser.add_argument(
    '-i', '--input', required=True, help='sam input file')
parser.add_argument(
    '-l', '--cb_list', required=True, help='cell barcode list')
parser.add_argument(
    '-od', '--output_directory', required=True, help='the directory for all the sam output files')

args = parser.parse_args()


cb = re.compile(r'(?<=CB:Z:).+?(?=\t)')


def read_cblist(cblf):
    with open(cblf) as f:
        cbl = [l.rstrip() for l in f.readlines()]
    return cbl


def process_file(input_, cbl, output_dir):
    if os.path.exists(output_dir):
        raise ValueError('Output directory already exists.')
    else:
        os.makedirs(output_dir)
    headers = []
    ids = set()
    with open(input_) as f:
        for l in f:
            if l.startswith('@HD'):
                headers.append(l)
            elif l.startswith('@SQ'):
                headers.append(l)
            elif not l.startswith('@'):
                id = cb.findall(l)
                if len(id) == 1:
                    id = id[0]
                    if id in cbl:
                        if id not in ids:
                            ids.add(id)
                            with open(os.path.join(output_dir, '%s.sam' % id), 'w') as f:
                                f.writelines(headers)
                                f.write(l)
                        else:
                            with open(os.path.join(output_dir, '%s.sam' % id), 'a') as f:
                                f.write(l)


if __name__ == '__main__':
    cbl = read_cblist(args.cb_list)
    process_file(args.input, cbl, args.output_directory)
