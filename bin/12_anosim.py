import pandas as pd
import numpy as np
import os
from scipy.spatial.distance import pdist
import read_config as rc


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


def get_df(pre='ln', file='combined_alpha_m1.txt'):
    paths = rc.get_items('path')
    df = pd.read_csv(os.path.join(paths[1][1], file), sep='\t', index_col='OTU ID').drop('taxonomy', axis=1).T
    if pre == 'ln':
        df = np.log(df+0.00001)
    elif pre == 'log10':
        df = np.log10(df+0.00001)

    return df


def compute_dist(df, map_dict, dist_type):
    dists = []
    values_list = []
    for i in range(df.shape[0]-1):
        d = []
        for j in range(i+1, df.shape[0]):
            dist = pdist(df.iloc[[i, j], :], dist_type)[0]
            d.append(dist)
            values_list.append(dist)
        dists.append(d)

    values_list = sorted(values_list)
    for i in range(len(dists)):
        for j in range(len(dists[i])):
            dists[i][j] = values_list.index(dists[i][j])

    rb, nb, rw, nw = 0, 0, 0, 0
    for i in range(len(dists)):
        for j in range(len(dists[i])):
            if map_dict[df.index.tolist()[len(df.index.tolist())-len(dists[i])-1]] == map_dict[df.index.tolist()[j+len(df.index.tolist())-len(dists[i])]]:
                rw += dists[i][j]
                nw += 1
            else:
                rb += dists[i][j]
                nb += 1

    r = 4 * (rb/nb - rw/nw) / (len(df.index.tolist()) * (len(df.index.tolist())-1))
    p = permutation_test(r, dists, df.index.tolist(), map_dict)
    return r, p


def permutation_test(r, dists, index_list, map_dict):
    m = 0
    for k in range(999):
        index_list = np.random.permutation(index_list)
        rb, nb, rw, nw = 0, 0, 0, 0
        for i in range(len(dists)):
            for j in range(len(dists[i])):
                if map_dict[index_list[len(index_list) - len(dists[i]) - 1]] == map_dict[index_list[j + len(index_list) - len(dists[i])]]:
                    rw += dists[i][j]
                    nw += 1
                else:
                    rb += dists[i][j]
                    nb += 1

        r_test = 4 * (rb / nb - rw / nw) / (len(index_list) * (len(index_list) - 1))
        if r_test > r:
            m += 1

    p = m / 999
    if p < 0.0001:
        p = 0.0001
    return p


def anosim(df, map_dict, compare_groups, output_path):
    with open(os.path.join(output_path, 'anosim_result.txt'), 'w') as f:
        f.write('ps: R值域[-1, 1]，R=0表示没有差异，R>0表示组间差异大于组内差异，R<0表示组内差异大于组间差异\np<0.05表示具有显著性差异\n\n')
        f.write('Group0\tGroup1\tR值\tp-vaule\n')
        for i in range(len(compare_groups)):
            if len(compare_groups[i]) == 2:
                d0 = df[df.index.map(map_dict) == compare_groups[i][0]]
                d1 = df[df.index.map(map_dict) == compare_groups[i][1]]
                d = d0.append(d1, ignore_index=False)
                r, p = compute_dist(d, map_dict, dist_type='braycurtis')
                f.write('{}\t{}\t{:.2f}\t{}\n'.format(compare_groups[i][0], compare_groups[i][1], r, p))


def main():
    paths = rc.get_items('path')
    map_dict, compare_groups = get_map()
    df = get_df('ln', 'combined_otu_table_m2_std.txt')
    anosim(df, map_dict, compare_groups, paths[1][1])


if __name__ == '__main__':
    main()
