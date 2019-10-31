import os
import re
import random
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import read_config as rc


def rarefaction_curve(input_path):
    plot_data = []
    for roots, dirs, files in os.walk(input_path):
        if len(re.findall('qiime_silva99v2$', roots)) > 0:
            for file in files:
                if file == 'otu_table_m2.txt':
                    reads_list = []
                    with open(os.path.join(roots, file), 'r') as f:
                        for line in f:
                            if len(re.findall('^(?:#|)OTU ID', line)) > 0:
                                sample_id = re.findall('^(?:#|)OTU ID\t(.*)\ttaxonomy', line)[0]
                                continue

                            aa = line.strip().split('\t')
                            if int(aa[1]) == 1:
                                continue

                            for i in range(int(aa[1])):
                                reads_list.append(aa[0])

                    random.seed(1)
                    random.shuffle(reads_list)
                    n_otu, n_read, otu = 0, 0, {}
                    x, y = [], []
                    for i in range(len(reads_list)):
                        n_read += 1
                        if reads_list[i] not in otu:
                            n_otu += 1
                            otu.setdefault(reads_list[i], 0)
                            otu[reads_list[i]] += 1
                            x.append(n_read)
                            y.append(n_otu)

                    x_new = np.arange(x[0], x[-1], len(x)/10)
                    fun = interpolate.interp1d(x, y, kind='linear')
                    y_new = fun(x_new)
                    plot_data.append([x_new.tolist(), y_new.tolist(), sample_id])
                    # plot_data.append([x, y, sample_id])

    sort_plot_data = [value for index, value in sorted(enumerate(plot_data), key=lambda plot_data:plot_data[1][2])]
    return sort_plot_data


def draw_plot(plot_data, output_path):
    plt.figure(figsize=(10, 8))
    for i in range(len(plot_data)):
        plt.plot(plot_data[i][0], plot_data[i][1], label=plot_data[i][2])

    # plt.xlim(0, 30000)
    plt.xticks(fontproperties='Times New Roman', size=14)
    plt.yticks(fontproperties='Times New Roman', size=14)
    plt.title('Rarefaction curve', pad=10, fontdict={'family': 'Times New Roman', 'size': 22})
    if len(plot_data) <= 20:
        plt.legend(loc=[1.01, 0.1], prop={'family': 'Times New Roman', 'size': 14})

    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'rarefaction_curve.pdf'))
    plt.close()


def main():
    paths = rc.get_items('path')
    plot_data = rarefaction_curve(paths[0][1])
    draw_plot(plot_data, paths[1][1])


if __name__ == '__main__':
    main()
