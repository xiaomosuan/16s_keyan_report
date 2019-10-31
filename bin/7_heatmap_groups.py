import pandas as pd
import numpy as np
import re
import os
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import read_config as rc


tax_dict = {
    'k': 'Kingdom',
    'p': 'Phylum',
    'c': 'Class',
    'o': 'Order',
    'f': 'Family',
    'g': 'Genus',
    's': 'Species',
}


def args_params():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input_file', help='input_file')
    parser.add_argument('-m', dest='map_file', help='map_file')
    args = parser.parse_args()
    return args.input_file, args.map_file


def get_map():
    maps = rc.get_items('map')
    map_dict = {}
    for i in range(len(maps)):
        map_dict[maps[i][0]] = maps[i][1]

    compares = rc.get_items('compare')
    compare_groups = []
    for i in range(len(compares)):
        aa = compares[i][1].replace(' ', '').split(',')
        compare_groups.append([aa[0], aa[1]])

    return map_dict, compare_groups


def get_df(output_path, pre='ln'):
    df = pd.read_csv(os.path.join(output_path, 'combined_otu_table_m3_std.txt'), sep='\t', index_col='OTU ID').drop('taxonomy', axis=1)
    if pre == 'ln':
        df = np.log(df+1)
    elif pre == 'log10':
        df = np.log10(df+1)

    return df


def get_top50(df, map_dict, output_path):
    df['total'] = df.apply(lambda x: x.sum(), axis=1)
    df.sort_values('total', inplace=True, ascending=False)
    samples = df.columns.tolist()
    samples.pop()
    for tax in ['p', 'c', 'o', 'f', 'g', 's']:
        top_list, tax_list, others_list, num = [], [], [], 0
        for i in range(df.shape[0]):
            if re.match('^{}__'.format(tax), df.index[i]):
                num += 1
                single_list = [x for x in df.iloc[i, :-1]]
                if num <= 50:
                    top_list.append(single_list)
                    tax_list.append(df.index[i])
                else:
                    others_list.append(single_list)

        top_array = np.asarray(top_list)
        draw_cluster_map(top_array, samples, tax_list, tax, map_dict, output_path)


def draw_cluster_map(top_array, samples, tax_list, tax, map_dict, output_path):
    top_array = top_array.T
    df = pd.DataFrame(top_array, columns=tax_list, index=samples)
    df['group'] = df.index.map(map_dict)
    groups = df.pop('group')
    lut = dict(zip(groups.unique(), sns.color_palette("Set1")))
    col_colors = groups.map(lut)
    df = df.T
    # braycurtis euclidean
    sns.clustermap(df, linewidths=0.5, col_colors=col_colors, metric='braycurtis')
    # plt.legend(handles=[], labels=[], loc='best')
    plt.savefig(os.path.join(output_path, 'heatmap_top50_{}.pdf'.format(tax_dict[tax])), bbox_inches='tight')
    plt.close()

    # top_array = top_array.T
    # df = pd.DataFrame(top_array, columns=tax_list, index=samples)
    # df = np.log10(df + 1)
    # df['group'] = df.index.map(map_dict)
    # groups = df.pop('group')
    # lut = dict(zip(groups.unique(), sns.color_palette("Set1")))
    # row_colors = groups.map(lut)
    # sns.clustermap(df, linewidths=0.5, row_colors=row_colors, metric='braycurtis')
    # plt.savefig('heatmap_top50_{}.pdf'.format(tax_dict[tax]), bbox_inches='tight')
    # plt.close()


def main():
    paths = rc.get_items('path')
    map_dict, _ = get_map()
    df = get_df(paths[1][1])
    get_top50(df, map_dict, paths[1][1])


if __name__ == '__main__':
    main()
