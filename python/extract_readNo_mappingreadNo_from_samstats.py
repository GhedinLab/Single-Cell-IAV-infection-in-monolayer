import os
import argparse
import re
import pandas as pd

parser = argparse.ArgumentParser('Extract # of reads and # of mapped reads for each cell from samtools stats output.')

parser.add_argument(
    'files', nargs='+', help='input sam file for each cell')
parser.add_argument(
    '-o', '--output', required=True, help='the output file')

args = parser.parse_args()


raw = re.compile(r'(?<=SN\tsequences:\t).+?(?=\n)')
mapped = re.compile(r'(?<=SN\treads mapped:\t).+?(?=\n)')


def retrieve_info(inputs):
    metrics = {}
    for i in inputs:
        cid = os.path.basename(i).split('.')[0][:-2]
        with open(i) as f:
            for l in f:
                if l.startswith('SN\tsequences:'):
                    reads = raw.findall(l)[0]
                elif l.startswith('SN\treads mapped:'):
                    mreads = mapped.findall(l)[0]
        if cid not in metrics:
            metrics.update({cid: [reads, mreads]})
        else:
            cid = '-'.join([cid, '2'])
            metrics.update({cid: [reads, mreads]})
    return metrics


def main():
    metrics = retrieve_info(args.files)
    mdf = pd.DataFrame.from_dict(metrics, orient='index')
    mdf.columns = ['reads', 'mapped_reads']
    mdf.to_csv(args.output)


if __name__ == '__main__':
    main()
