import argparse

import pandas as pd


parser = argparse.ArgumentParser(description='Group records.')

parser.add_argument(
    '--input_file', required=True, help='input file')

parser.add_argument(
    '--output_file', required=True, help='output file')

args = parser.parse_args()


class Group:
    def __init__(self, row):
        self.rows = [row]

    def get_max(self, key, rows=None):
        if rows is None:
            rows = self.rows
        return max([row[key] for row in rows])

    def get_min(self, key, rows=None):
        if rows is None:
            rows = self.rows
        return min([row[key] for row in rows])

    def add(self, new_row):
        new_rows = self.rows + [new_row]

        if self.get_max('5p', new_rows) - self.get_min('5p', new_rows) > 10:
            return False

        if self.get_max('3p', new_rows) - self.get_min('3p', new_rows) > 10:
            return False

        self.rows = new_rows
        return True

    @property
    def df(self):
        return pd.DataFrame(self.rows)


if __name__ == '__main__':
    data = pd.read_csv(args.input_file, index_col=0)
    data['5p'] = data['boundary_comb'].apply(lambda x: int(x.split('-')[0]))
    data['3p'] = data['boundary_comb'].apply(lambda x: int(x.split('-')[1]))

    data.sort_values(['5p', '3p'], inplace=True)
    data.reset_index(drop=True, inplace=True)

    groups = []
    for i, row in data.iterrows():
        grouped = False
        for group in groups:
            if group.add(row):
                grouped = True
                break
        if not grouped:
            groups.append(Group(row))

    dfs = []
    for i, group in enumerate(groups):
        df = group.df
        df['group'] = i + 1
        dfs.append(df)

    master = pd.concat(dfs, ignore_index=True)
    master.to_csv(args.output_file, index=False)
