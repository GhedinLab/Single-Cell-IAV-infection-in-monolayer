import os
import argparse


# argument parser
parser = argparse.ArgumentParser('Fill IAV sam files from STAR output to exclude multi-mapper (MAPQ != 255).')

parser.add_argument(
    'files', nargs='+', help='input sam file for each cell')
parser.add_argument(
    '-od', '--output_dir', required=True, help='the directory for all the output files')

args = parser.parse_args()


def process_file(file_):
    headers = []
    data = []
    with open(file_) as f:
        cid = os.path.basename(file_).split('.')[0]
        for l in f:
            if l.startswith('@'):
                headers.append(l)
            else:
                if l.split('\t')[4] == '255':
                    data.append(l)
    return headers, data, cid


def main():
    for file_ in args.files:
        headers, data, cid = process_file(file_)
        with open(os.path.join(args.output_dir, '%s.sam' % cid), 'w') as f:
            f.writelines(headers + data)


if __name__ == '__main__':
    main()
