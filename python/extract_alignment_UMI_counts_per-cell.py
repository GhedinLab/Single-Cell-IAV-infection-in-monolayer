import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser('Extract # of the UMI counts of gapped alignments for each cell.')

parser.add_argument(
    'inputs', nargs='+', help='input count txt file for each cell')
parser.add_argument(
    '-o', '--output', required=True, help='the output file')

args = parser.parse_args()


def main():
    for i, input_ in enumerate(args.inputs):
        cid = os.path.basename(input_).split('_')[0]
        if i == 0:
            df = pd.read_csv(input_, sep='\t', names=['seg', cid])
        else:
            new = pd.read_csv(input_, sep='\t', names=['seg', cid])
            df = pd.merge(df, new, on='seg', how='outer')
    df.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()