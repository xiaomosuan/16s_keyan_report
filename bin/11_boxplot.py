import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import read_config as rc


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


def get_df(map_dict, pre='ln', file='combined_alpha_m1.txt'):
    paths = rc.get_items('path')
    df = pd.read_csv(os.path.join(paths[1][1], file), sep='\t', index_col='Sample')
    if pre == 'ln':
        df = np.log(df+0.00001)
    elif pre == 'log10':
        df = np.log10(df+0.00001)

    df = df.reset_index()
    df['Type'] = df['Sample'].map(map_dict)
    return df


def draw_boxplot_all(df, compare_group, output_path):
    fig, axes = plt.subplots(2, 2)
    axe1, axe2, axe3, axe4 = axes.ravel()
    plt.suptitle('Alpha diversity')
    if len(compare_group) > 2:
        plt.subplots_adjust(wspace=0.5, hspace=1)
    else:
        plt.subplots_adjust(wspace=0.5, hspace=0.2)

    sns.boxplot(x='Type', y='shannon', data=df, width=0.4, order=compare_group, ax=axe1)
    sns.boxplot(x='Type', y='chao1', data=df, width=0.4, order=compare_group, ax=axe2)
    sns.boxplot(x='Type', y='observed_species', data=df, width=0.4, order=compare_group, ax=axe3)
    sns.boxplot(x='Type', y='simpson', data=df, width=0.4, order=compare_group, ax=axe4)
    axe1.set_xlabel('')
    axe2.set_xlabel('')
    axe3.set_xlabel('')
    axe4.set_xlabel('')
    if len(compare_group) > 2:
        axe1.set_xticklabels(labels=compare_group, rotation=75)
        axe2.set_xticklabels(labels=compare_group, rotation=75)
        axe3.set_xticklabels(labels=compare_group, rotation=75)
        axe4.set_xticklabels(labels=compare_group, rotation=75)

    plt.savefig(os.path.join(output_path, 'boxplot_all_{}.pdf'.format('_'.join([x for x in compare_group]))), bbox_inches='tight')
    plt.close()


def draw_boxplot(df, compare_group, output_path):
    for i in ['chao1', 'observed_species', 'shannon', 'simpson']:
        plt.figure(figsize=(6, 6))
        sns.boxplot(x='Type', y=i, data=df, width=0.4, order=compare_group)
        plt.xlabel('')
        plt.ylabel('{}'.format(i), fontdict={'family': 'Times New Roman', 'size': 16})

        plt.yticks(fontproperties='Times New Roman', size=14)
        plt.title('{}'.format(i), pad=10, fontdict={'family': 'Times New Roman', 'size': 20})
        if len(compare_group) > 2:
            plt.xticks(fontproperties='Times New Roman', size=14, rotation=75)
        else:
            plt.xticks(fontproperties='Times New Roman', size=14)

        plt.tight_layout()
        plt.savefig(os.path.join(output_path, 'boxplot_{}_{}.pdf'.format(i, '_'.join([x for x in compare_group]))))
        plt.close()


def main():
    input_path, output_path, map_dict, compare_groups = get_config_info()
    df = get_df(map_dict, 'ln')
    for i in range(len(compare_groups)):
        df1 = df[df['Type'] == compare_groups[i][0]]
        for j in range(1, len(compare_groups[i])):
            df2 = df[df['Type'] == compare_groups[i][j]]
            df1 = df1.append(df2)

        draw_boxplot_all(df, compare_groups[i], output_path)
        draw_boxplot(df, compare_groups[i], output_path)


if __name__ == '__main__':
    main()
