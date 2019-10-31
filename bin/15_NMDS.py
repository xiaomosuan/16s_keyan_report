import pandas as pd
import numpy as np
import os
import read_config as rc
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
"""
非度量多维尺度分析(Non-metric multidimensional scaling, NMDS)
"""


def get_config_info():
    paths = rc.get_items('path')

    maps = rc.get_items('map')
    map_dict = {}
    for i in range(len(maps)):
        map_dict[maps[i][0]] = maps[i][1]

    compares = rc.get_items('compare')
    compare_groups = []
    for i in range(len(compares)):
        aa = compares[i][1].replace(' ', '').split(',')
        compare_groups.append(aa)

    return paths[0][1], paths[1][1], map_dict, compare_groups


def get_df(output_path, pre='ln', file='combined_alpha_m1.txt'):
    df = pd.read_csv(os.path.join(output_path, file), sep='\t', index_col='OTU ID').drop('taxonomy', axis=1).T
    if pre == 'ln':
        df = np.log(df+0.00001)
    elif pre == 'log10':
        df = np.log10(df+0.00001)

    return df


def sk_mds(df, n):
    """
    euclidean
    """
    mds = MDS(n_components=n, metric=False)
    reduced_x = mds.fit_transform(df)
    return reduced_x


def draw_scatter_plot(reduced_x, map_dict, compare_group, d, output_path):
    y = d.index.tolist()
    c_list = ['red', 'blue', 'green', 'darkorange', 'grey', 'hotpink', 'skyblue', 'violet', 'yellow', 'black']
    marker_list = ['o', '^', 's', 'd', '*', 'v', 'x', '+', 'p', '8', 'h', '1']
    for i in range(len(compare_group)):
        exec('x_{}, y_{} = [], []'.format(compare_group[i], compare_group[i]))

    for i in range(len(reduced_x)):
        exec('x_{}.append(reduced_x[i][0])'.format(map_dict[y[i]]))
        exec('y_{}.append(reduced_x[i][1])'.format(map_dict[y[i]]))

    plt.figure(figsize=[8, 8])
    for i in range(len(compare_group)):
        exec("plt.scatter(x_{}, y_{}, c='{}', marker='{}', label='{}', s=100)".format(
            compare_group[i], compare_group[i], c_list[i], marker_list[i], compare_group[i]))

    plt.xlabel('NMDS1', fontdict={'family': 'Times New Roman', 'size': 16})
    plt.ylabel('NMDS2', fontdict={'family': 'Times New Roman', 'size': 16})
    plt.xticks(fontproperties='Times New Roman', size=14)
    plt.yticks(fontproperties='Times New Roman', size=14)
    plt.title('NMDS', pad=10, fontdict={'family': 'Times New Roman', 'size': 22})
    plt.legend(prop={'family': 'Times New Roman', 'size': 16})
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'NMDS_{}.pdf'.format('_'.join([x for x in compare_group]))))


def main():
    _, output_path, map_dict, compare_groups = get_config_info()
    df = get_df(output_path, 'ln', 'combined_otu_table_m2_std.txt')
    for i in range(len(compare_groups)):
        df1 = df[df.index.map(map_dict) == compare_groups[i][0]]
        for j in range(1, len(compare_groups[i])):
            df2 = df[df.index.map(map_dict) == compare_groups[i][j]]
            df1 = df1.append(df2)

        reduced_x = sk_mds(df1, 2)
        draw_scatter_plot(reduced_x, map_dict, compare_groups[i], df1, output_path)


if __name__ == '__main__':
    main()
