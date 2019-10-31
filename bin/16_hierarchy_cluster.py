import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
import read_config as rc
import os


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


def get_group_col(df, map_dict, compare_group):
    c_list = ['red', 'blue', 'green', 'darkorange', 'yellow', 'hotpink', 'grey', 'skyblue', 'violet', 'cadetblue', 'black']
    group_col = {}
    for i in range(len(compare_group)):
        group_col[compare_group[i]] = c_list[i]
    sample_col = {}
    for i in range(len(df.index.tolist())):
        sample_col[df.index.tolist()[i]] = group_col[map_dict[df.index.tolist()[i]]]

    return sample_col


def upgma(df, sample_col, compare_group, output_path):
    dis = pdist(df, metric='braycurtis')
    z = hierarchy.linkage(dis, method='average')
    hierarchy.dendrogram(z, labels=df.index)
    plt.xticks(fontproperties='Times New Roman', size=14, rotation=75)
    plt.title('UPGMA', pad=10, fontdict={'family': 'Times New Roman', 'size': 22})
    ax = plt.gca()
    xlbls = ax.get_xmajorticklabels()
    for lbl in xlbls:
        lbl.set_color(sample_col[lbl.get_text()])
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'UPGMA_{}.pdf'.format('_'.join([x for x in compare_group]))))
    plt.close()


def wpgma(df, sample_col, compare_group, output_path):
    dis = pdist(df, metric='braycurtis')
    z = hierarchy.linkage(dis, method='weighted')
    hierarchy.dendrogram(z, labels=df.index)
    plt.xticks(fontproperties='Times New Roman', size=14, rotation=75)
    plt.title('WPGMA', pad=10, fontdict={'family': 'Times New Roman', 'size': 22})
    ax = plt.gca()
    xlbls = ax.get_xmajorticklabels()
    for lbl in xlbls:
        lbl.set_color(sample_col[lbl.get_text()])
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'WPGMA_{}.pdf'.format('_'.join([x for x in compare_group]))))
    plt.close()


def main():
    _, output_path, map_dict, compare_groups = get_config_info()
    df = get_df(output_path, 'ln', 'combined_otu_table_m2_std.txt')
    for i in range(len(compare_groups)):
        df1 = df[df.index.map(map_dict) == compare_groups[i][0]]
        for j in range(1, len(compare_groups[i])):
            df2 = df[df.index.map(map_dict) == compare_groups[i][j]]
            df1 = df1.append(df2)

        sample_col = get_group_col(df1, map_dict, compare_groups[i])
        upgma(df1, sample_col, compare_groups[i], output_path)
        wpgma(df1, sample_col, compare_groups[i], output_path)


if __name__ == '__main__':
    main()
