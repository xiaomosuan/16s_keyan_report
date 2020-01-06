import pandas as pd
import numpy as np
import re
import os
import matplotlib.pyplot as plt
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
print ("wo hao diao,hhhhhhhh")


def get_top10(df, output_path):
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
                if num <= 10:
                    top_list.append(single_list)
                    tax_list.append(df.index[i])
                else:
                    others_list.append(single_list)

        if num > 10:
            tax_list.append('other {}'.format(tax_dict[tax]))
            others_array = np.asarray(others_list)
            others = np.sum(others_array, axis=0).tolist()
            top_list.append(others)

        top_array = np.asarray(top_list)
        col_sum = np.sum(top_array, axis=0)
        b_array = np.zeros(top_array.shape)
        for i in range(b_array.shape[0]):
            for j in range(b_array.shape[1]):
                b_array[i][j] = top_array[i][j] / col_sum[j]

        draw_stack_bar(b_array, samples, tax_list, tax, output_path)


def draw_stack_bar(b_array, samples, tax_list, tax, output_path):
    c_list = ['red', 'blue', 'green', 'darkorange', 'yellow', 'hotpink', 'grey', 'skyblue', 'violet', 'cadetblue', 'black']
    plt.figure(figsize=[12, 8])
    plt.bar(samples, b_array[0], align='center', color=c_list[0], label=tax_list[0])
    for i in range(1, b_array.shape[0]):
        if i == 1:
            bottom = b_array[0]
        else:
            bottom = bottom + b_array[i-1]

        plt.bar(samples, b_array[i], bottom=bottom, align='center', color=c_list[i], label=tax_list[i])

    plt.xlabel('')
    plt.ylabel('Relative Abundance', fontdict={'family': 'Times New Roman', 'size': 16})
    plt.xticks(fontproperties='Times New Roman', size=12, rotation=75)
    plt.yticks(fontproperties='Times New Roman', size=14)
    plt.title('Top10 {}'.format(tax_dict[tax]), pad=10, fontdict={'family': 'Times New Roman', 'size': 22})
    plt.legend(loc=[1.01, 0.3], prop={'family': 'Times New Roman', 'size': 16})
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'stack_top10_{}.pdf'.format(tax_dict[tax])))


def main():
    paths = rc.get_items('path')
    df = pd.read_csv(os.path.join(paths[1][1], 'combined_otu_table_m3_std.txt'), sep='\t', index_col='OTU ID').drop('taxonomy', axis=1)
    get_top10(df, paths[1][1])


if __name__ == '__main__':
    main()
