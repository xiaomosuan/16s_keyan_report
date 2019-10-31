from scipy import stats
import pandas as pd
import os
import numpy as np
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


def get_df(map_dict, pre='ln', file='combined_alpha_m1.txt'):
    paths = rc.get_items('path')
    df = pd.read_csv(os.path.join(paths[1][1], file), sep='\t', index_col='Sample')
    if pre == 'ln':
        df = np.log(df+0.00001)
    elif pre == 'log10':
        df = np.log10(df+0.00001)

    df = df.reset_index()
    df['Type'] = df['Sample'].map(map_dict)
    df = df.set_index('Type')
    return df


def alpha_stats_kruskal(df, compare_groups, output_path):
    """
    Kruskal-Wallis H-test
    """
    f = open(os.path.join(output_path, 'alpha_kruskal_test.txt'), 'w')
    print('method: Kruskal-Wallis H-test', file=f)
    print('ps: p-value小于0.05表示Group0和Group1在Alpha(xx)中具有显著性差异\n', file=f)
    print('Alpha\tGroup0\tGroup1\tH-statistic\tp-value', file=f)
    for i in ['chao1', 'observed_species', 'shannon', 'simpson']:
        for j in range(len(compare_groups)):
            group_list = []
            for g in range(len(compare_groups[j])):
                group_list.append([x for x in df.loc[compare_groups[j][g], i]])

            w, p = stats.kruskal(group_list[0], group_list[1])
            print('{}\t{}\t{}\t{:.2f}\t{:.4f}'.format(i, compare_groups[j][0], compare_groups[j][1], w, p), file=f)

    f.close()


def alpha_stats_ranksums(df, compare_groups, output_path):
    """
    Wilcoxon rank-sum statistic for two samples
    """
    f = open(os.path.join(output_path, 'alpha_ranksums_test.txt'), 'w')
    print('method: Wilcoxon rank-sum test', file=f)
    print('ps: p-value小于0.05表示Group0和Group1在Alpha(xx)中具有显著性差异\n', file=f)
    print('Alpha\tGroup0\tGroup1\tstatistic\tp-value', file=f)
    for i in ['chao1', 'observed_species', 'shannon', 'simpson']:
        for j in range(len(compare_groups)):
            group_list = []
            for g in range(len(compare_groups[j])):
                group_list.append([x for x in df.loc[compare_groups[j][g], i]])

            w, p = stats.ranksums(group_list[0], group_list[1])
            print('{}\t{}\t{}\t{:.2f}\t{:.4f}'.format(i, compare_groups[j][0], compare_groups[j][1], w, p), file=f)

    f.close()


def alpha_stats_f_oneway(df, compare_groups, output_path):
    """
    one-way ANOVA tests
    """
    f = open(os.path.join(output_path, 'alpha_f_oneway_test.txt'), 'w')
    print('method: one-way ANOVA tests', file=f)
    print('ps: p-value小于0.05表示Group0和Group1在Alpha(xx)中具有显著性差异\n', file=f)
    print('Alpha\tGroup0\tGroup1\tF-value\tp-value', file=f)
    for i in ['chao1', 'observed_species', 'shannon', 'simpson']:
        for j in range(len(compare_groups)):
            group_list = []
            for g in range(len(compare_groups[j])):
                group_list.append([x for x in df.loc[compare_groups[j][g], i]])

            w, p = stats.f_oneway(group_list[0], group_list[1])
            print('{}\t{}\t{}\t{:.2f}\t{:.4f}'.format(i, compare_groups[j][0], compare_groups[j][1], w, p), file=f)

    f.close()


def alpha_stats_mannwhitneyu(df, compare_groups, output_path):
    """
    Mann-Whitney rank test on samples x and y
    """
    f = open(os.path.join(output_path, 'alpha_mannwhitneyu_test.txt'), 'w')
    print('method: Mann-Whitney rank test', file=f)
    print('ps: p-value小于0.05表示Group0和Group1在Alpha(xx)中具有显著性差异\n', file=f)
    print('Alpha\tGroup0\tGroup1\tU-statistic\tp-value', file=f)
    for i in ['chao1', 'observed_species', 'shannon', 'simpson']:
        for j in range(len(compare_groups)):
            group_list = []
            for g in range(len(compare_groups[j])):
                group_list.append([x for x in df.loc[compare_groups[j][g], i]])

            w, p = stats.mannwhitneyu(group_list[0], group_list[1])
            print('{}\t{}\t{}\t{:.2f}\t{:.4f}'.format(i, compare_groups[j][0], compare_groups[j][1], w, p), file=f)

    f.close()


def alpha_stats_levene(df, compare_groups, output_path):
    """
    Levene test for equal variances
    """
    f = open(os.path.join(output_path, 'alpha_levene_test.txt'), 'w')
    print('method: Levene test', file=f)
    print('ps: p-value小于0.05表示Group0和Group1在Alpha(xx)中方差是不相等的\n', file=f)
    print('Alpha\tGroup0\tGroup1\tstatistic\tp-value', file=f)
    for i in ['chao1', 'observed_species', 'shannon', 'simpson']:
        for j in range(len(compare_groups)):
            group_list = []
            for g in range(len(compare_groups[j])):
                group_list.append([x for x in df.loc[compare_groups[j][g], i]])

            w, p = stats.levene(group_list[0], group_list[1])
            print('{}\t{}\t{}\t{:.2f}\t{:.4f}'.format(i, compare_groups[j][0], compare_groups[j][1], w, p), file=f)

    f.close()


def alpha_stats_ttest_ind(df, compare_groups, output_path):
    """
    T-test for the means of two independent samples of scores.
    """
    f = open(os.path.join(output_path, 'alpha_ttest_ind_test.txt'), 'w')
    print('method: ttest_ind ', file=f)
    print('ps: p-value小于0.05表示Group0和Group1在Alpha(xx)中方差是不相等的\n', file=f)
    print('Alpha\tGroup0\tGroup1\tt-statistic\tp-value', file=f)
    for i in ['chao1', 'observed_species', 'shannon', 'simpson']:
        for j in range(len(compare_groups)):
            group_list = []
            for g in range(len(compare_groups[j])):
                group_list.append([x for x in df.loc[compare_groups[j][g], i]])

            w, p = stats.ttest_ind(group_list[0], group_list[1], equal_var=False)
            print('{}\t{}\t{}\t{:.2f}\t{:.4f}'.format(i, compare_groups[j][0], compare_groups[j][1], w, p), file=f)

    f.close()


def main():
    paths = rc.get_items('path')
    map_dict, compare_groups = get_map()
    df = get_df(map_dict, 'ln')
    alpha_stats_kruskal(df, compare_groups, paths[1][1])
    alpha_stats_ranksums(df, compare_groups, paths[1][1])
    alpha_stats_f_oneway(df, compare_groups, paths[1][1])
    alpha_stats_mannwhitneyu(df, compare_groups, paths[1][1])
    alpha_stats_levene(df, compare_groups, paths[1][1])
    alpha_stats_ttest_ind(df, compare_groups, paths[1][1])


if __name__ == '__main__':
    main()
