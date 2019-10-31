import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import read_config as rc


def rank_abundance(input_path, output_path):
    plt.figure(figsize=(8, 7))
    num = 0
    for roots, dirs, files in os.walk(input_path):
        if len(re.findall('qiime_silva99v2$', roots)) > 0:
            for file in files:
                if file == 'otu_table_m1.txt':
                    num += 1
                    sample_id = roots.split('/')[-3]
                    df = pd.read_csv(os.path.join(roots, file), sep='\t').drop(['OTU ID', 'taxonomy'], axis=1)
                    # df[sample_id] = df[sample_id] / 1000000
                    df_sum = df[sample_id].sum()
                    df[sample_id] = df[sample_id].apply(lambda x: x / df_sum)
                    df.sort_values(sample_id, inplace=True, ascending=False)
                    df = df.reset_index().drop('index', axis=1)
                    plt.semilogy(df.index, df[sample_id])

    plt.xlim(0, 4000)
    plt.xticks(fontproperties='Times New Roman', size=14, rotation=0)
    plt.yticks(fontproperties='Times New Roman', size=14)
    plt.xlabel('OTU rank', fontdict={'family': 'Times New Roman', 'size': 16})
    plt.ylabel('Relative Abundance', fontdict={'family': 'Times New Roman', 'size': 16})
    plt.title('Rank-abundance distribution curve', pad=10, fontdict={'family': 'Times New Roman', 'size': 22})
    if num <= 20:
        plt.legend(loc=[1.01, 0.3], prop={'family': 'Times New Roman', 'size': 12})
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'rank_abundance.pdf'))
    plt.close()


def main():
    paths = rc.get_items('path')
    rank_abundance(paths[0][1], paths[1][1])


if __name__ == '__main__':
    main()
