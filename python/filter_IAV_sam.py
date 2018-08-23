import os
import argparse


# argument parser
parser = argparse.ArgumentParser('Extract reads mapped to IAV reference genomes.')

parser.add_argument(
    'files', nargs='+', help='input sam file for each cell')
parser.add_argument(
    '-g', '--reference_genome', required=True, help='IAV reference genomes in fasta format')
parser.add_argument(
    '-od', '--output_dir', required=True, help='the directory for all the output files')

args = parser.parse_args()


def read_fasta_ids(path):
    ids = set()
    with open(path) as f:
        for l in f:
            if l.startswith('>'):
                seqid = l.split(' ')[0][1:]
                ids.add(seqid)
    return ids


def process_file(file_, fasta_ids):
    headers = []
    data = []
    with open(file_) as f:
        cid = os.path.basename(file_).split('-')[0]
        for l in f:
            if l.startswith('@HD'):
                headers.append(l)
            elif l.startswith('@SQ'):
                if l.split('\t')[1][3:] in fasta_ids:
                    headers.append(l)
            elif l.startswith('@PG'):
                headers.append(l)
            elif not l.startswith('@'):
                if l.split('\t')[2] in fasta_ids:
                    data.append(l)
    return headers, data, cid


def main():
    ids = read_fasta_ids(args.reference_genome)
    for file_ in args.files:
        headers, data, cid = process_file(file_, ids)
        with open(os.path.join(args.output_dir, '%s.sam' % cid), 'w') as f:
            f.writelines(headers + data)


if __name__ == '__main__':
    main()
