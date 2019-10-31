import pandas as pd
import numpy as np
import os
import read_config as rc
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
"""
PCoA(principle coordination analysis)主坐标分析
otu表，map表
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


def compute_dist(df, dist_type='braycurtis'):
    dists = []
    for i in range(df.shape[0]):
        d = []
        for j in range(df.shape[0]):
            if i == j:
                dist = 0
            else:
                dist = pdist(df.iloc[[i, j], :], dist_type)[0]

            d.append(dist)
        dists.append(d)
    return dists


def mds(dists, n):
    d = np.asarray(dists)
    d_square = d ** 2
    total_mean = np.mean(d_square)
    column_mean = np.mean(d_square, axis=0)
    row_mean = np.mean(d_square, axis=1)
    B = np.zeros(d_square.shape)
    for i in range(B.shape[0]):
        for j in range(B.shape[1]):
            B[i][j] = -0.5 * (d_square[i][j] - row_mean[i] - column_mean[j] + total_mean)

    # 求特征值及特征向量
    value, vector = np.linalg.eigh(B)
    # 对特征值进行排序，得到排序索引
    value_sorted_indices = np.argsort(value)
    # 提取d个最大特征向量, -d-1前加:才能向左切
    top_vector = vector[:, value_sorted_indices[:-(n + 1):-1]]
    top_value = value[value_sorted_indices[:-(n + 1):-1]]
    top_ratio = [x / np.sum(value) for x in top_value]
    reduced_x = np.dot(top_vector, np.sqrt(np.diag(top_value)))
    return reduced_x, top_ratio


def draw_scatter_plot(reduced_x, top_ratio, map_dict, compare_group, d, output_path):
    y = d.index.tolist()
    c_list = ['red', 'blue', 'green', 'darkorange', 'grey', 'hotpink', 'skyblue', 'violet', 'yellow', 'black']
    marker_list = ['o', '^', 's', 'd', '*', 'v', 'x', '+', 'p', '8', 'h', '1']
    for i in range(len(compare_group)):
        exec('x_{}, y_{} = [], []'.format(compare_group[i], compare_group[i]))

    for i in range(len(reduced_x)):
        exec('x_{}.append(reduced_x[i][0])'.format(map_dict[y[i]]))
        exec('y_{}.append(reduced_x[i][1])'.format(map_dict[y[i]]))

    plt.figure(figsize=(8, 8))
    for i in range(len(compare_group)):
        exec("plt.scatter(x_{}, y_{}, c='{}', marker='{}', label='{}', s=100)".format(
            compare_group[i], compare_group[i], c_list[i], marker_list[i], compare_group[i]))

    plt.xlabel('PC1 ( {:.2f}% )'.format(top_ratio[0] * 100), fontdict={'family': 'Times New Roman', 'size': 16})
    plt.ylabel('PC2 ( {:.2f}% )'.format(top_ratio[1] * 100), fontdict={'family': 'Times New Roman', 'size': 16})
    plt.xticks(fontproperties='Times New Roman', size=14)
    plt.yticks(fontproperties='Times New Roman', size=14)
    plt.title('PCoA', pad=10, fontdict={'family': 'Times New Roman', 'size': 22})
    plt.legend(prop={'family': 'Times New Roman', 'size': 16})
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'PCoA_{}.pdf'.format('_'.join([x for x in compare_group]))))
    plt.close()


def main():
    _, output_path, map_dict, compare_groups = get_config_info()
    df = get_df(output_path, 'ln', 'combined_otu_table_m2_std.txt')

    for i in range(len(compare_groups)):
        df1 = df[df.index.map(map_dict) == compare_groups[i][0]]
        for j in range(1, len(compare_groups[i])):
            df2 = df[df.index.map(map_dict) == compare_groups[i][j]]
            df1 = df1.append(df2)
        # df1.to_csv('1.txt', index=True, header=True, sep='\t')
        dists = compute_dist(df1, dist_type='braycurtis')
        reduced_x, top_ratio = mds(dists, 2)
        draw_scatter_plot(reduced_x, top_ratio, map_dict, compare_groups[i], df1, output_path)


if __name__ == '__main__':
    main()
