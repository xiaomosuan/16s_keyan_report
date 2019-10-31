import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cross_decomposition import PLSRegression #偏最小二乘法
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


def get_df(output_path, pre='ln', file='combined_otu_table_m2_std.txt'):
    df = pd.read_csv(os.path.join(output_path, file), sep='\t', index_col='OTU ID').drop('taxonomy', axis=1).T
    if pre == 'ln':
        df = np.log(df + 0.00001) ##为了处理0？
    elif pre == 'log10':
        df = np.log10(df + 0.00001)

    return df


def pls_da(df, n):
    X = df.iloc[:, 0:-1]
    y = df.iloc[:, -1]
    # y_class = pd.get_dummies(y)
    plsda = PLSRegression(n_components=n)
    reduced_x = plsda.fit_transform(X, y)
    return reduced_x[0]


def draw_scatter_plot(reduced_x, map_dict, compare_group, df, output_path):
    y = df.index.tolist()
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

    plt.xlabel('COMPONENT 1', fontdict={'family': 'Times New Roman', 'size': 16})
    plt.ylabel('COMPONENT 2', fontdict={'family': 'Times New Roman', 'size': 16})
    plt.xticks(fontproperties='Times New Roman', size=14)
    plt.yticks(fontproperties='Times New Roman', size=14)
    plt.title('PLSDA', pad=10, fontdict={'family': 'Times New Roman', 'size': 22})
    plt.legend(prop={'family': 'Times New Roman', 'size': 16})
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'PLSDA_{}.pdf'.format('_'.join([x for x in compare_group]))))
    plt.close()


def main():
    _, output_path, map_dict, compare_groups = get_config_info()
    df = get_df(output_path, 'ln', '/Users/congliu/combined_otu_table_m2_std.txt')
    df['Type'] = df.index.map(map_dict)
    print(compare_groups) #[['Normal', 'CKD5']]

    for i in range(len(compare_groups)): #i==0
        df1 = df[df.index.map(map_dict) == compare_groups[i][0]]
        for j in range(1, len(compare_groups[i])):
            df2 = df[df.index.map(map_dict) == compare_groups[i][j]]
            df1 = df1.append(df2) #df1===df?

        compare_dict = {}
        for k in range(len(compare_groups[i])):
            compare_dict[compare_groups[i][k]] = k

        df1['group'] = df1['Type'].map(compare_dict)
        df1 = df1.drop('Type', axis=1)
        reduced_x = pls_da(df1, 2)
        #print(len(reduced_x))
        draw_scatter_plot(reduced_x, map_dict, compare_groups[i], df1, output_path)


if __name__ == '__main__':
    main()
