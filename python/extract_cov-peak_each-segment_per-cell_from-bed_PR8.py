import os
import argparse
import pandas as pd

# argument parser
parser = argparse.ArgumentParser('Extract the peak position in the bed-genomecov output for 3 polymerase segments in PR8.')

parser.add_argument(
    'files', nargs='+', help='input bed genomecov file for each cell')
parser.add_argument(
    '-o', '--output', required=True, help='the output csv file name')

args = parser.parse_args()


def main():
    peaks = {}
    for file_ in args.files:
        cid = os.path.basename(file_).split('.')[0]
        cov = {'AF389115.1': 0, 'AF389116.1': 0, 'AF389117.1': 0}
        pb2_max_cov = 0
        pb1_max_cov = 0
        pa_max_cov = 0
        with open(file_) as f:
            for l in f:
                segment = l.split('\t')[0]
                if int(l.split('\t')[1]) > 1000:
                    if segment == 'AF389115.1':
                        if int(l.split('\t')[2].split('\n')[0]) > pb2_max_cov:
                            pb2_max_cov = int(l.split('\t')[2].split('\n')[0])
                            cov['AF389115.1'] = l.split('\t')[1]
                    elif segment == 'AF389116.1':
                        if int(l.split('\t')[2].split('\n')[0]) > pb1_max_cov:
                            pb1_max_cov = int(l.split('\t')[2].split('\n')[0])
                            cov['AF389116.1'] = l.split('\t')[1]
                    elif segment == 'AF389117.1':
                        if int(l.split('\t')[2].split('\n')[0]) > pa_max_cov:
                            pa_max_cov = int(l.split('\t')[2].split('\n')[0])
                            cov['AF389117.1'] = l.split('\t')[1]
        peaks[cid] = cov
    peaks_df = pd.DataFrame(peaks)
    peaks_df.to_csv(args.output)


if __name__ == '__main__':
    main()
