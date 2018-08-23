import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description='Calculate TPM & FPKM')

parser.add_argument('--counts', required=True, help='raw counts')
parser.add_argument('--lengths', required=True, help='gene lengths')
parser.add_argument('--meanfsize', required=True, help='the mean of fragment size for each library')
parser.add_argument('--output1_TPM', required=True, help='output TPM file')
parser.add_argument('--output2_FPKM', required=True, help='output FPKM file')
parser.add_argument('--skip', type=int, default=8, help='rows of viral segments to skip when calculating library size')

args = parser.parse_args()

eps = 1e-300


def get_mu(col):
    data = pd.read_csv(args.meanfsize)
    data.set_index('Cell', inplace=True)
    col2mu = data['Mean_insert_size'].to_dict()
    return col2mu[col]


def get_tl(l, col):
    return np.maximum(l - get_mu(col) + 1, 0)


def get_tpm_and_fpkm(x, l, skip):
    x[l < eps] = 0
    x = x / x[skip:].sum()

    # calculate fpkm
    fpkm = x / l * 1e9
    fpkm[l < eps] = 0

    # calculate tpm
    tpm = fpkm / fpkm[skip:].sum() * 1e6
    return tpm, fpkm


if __name__ == '__main__':
    counts = pd.read_csv(args.counts, index_col=0)
    lengths = pd.read_csv(args.lengths, index_col=0)

    # make sure counts table does not have a col named 'Length'
    assert 'Length' not in counts.columns

    lengths.set_index('Geneid', inplace=True)

    master = pd.merge(
        counts, lengths, 'left',
        left_index=True, right_index=True)

    # make sure 'Length' col is all not null
    assert master.Length.notnull().all()

    tpms = {}
    fpkms = {}
    for col in master:
        if col != 'Length':
            tpm, fpkm = get_tpm_and_fpkm(
                master[col], get_tl(master['Length'], col), args.skip)
            tpms[col] = tpm
            fpkms[col] = fpkm

    tpms = pd.DataFrame(tpms)
    fpkms = pd.DataFrame(fpkms)
    tpms.to_csv(args.output1_TPM)
    fpkms.to_csv(args.output2_FPKM)
