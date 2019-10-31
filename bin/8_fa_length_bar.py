import os
import re
import numpy as np
import matplotlib.pyplot as plt
import read_config as rc


def compute_length_ratio(input_path, output_path):
    reads_length = np.zeros(8)
    for roots, dirs, files in os.walk(input_path):
        for file in files:
            if file == 'final.fa':
                with open(os.path.join(roots, file), 'r') as f1:
                    for line in f1:
                        if re.match('^>', line):
                            continue
                        line = line.strip()
                        if 200 <= len(line) < 250:
                            reads_length[0] += 1
                        elif 250 <= len(line) < 300:
                            reads_length[1] += 1
                        elif 300 <= len(line) < 350:
                            reads_length[2] += 1
                        elif 350 <= len(line) < 400:
                            reads_length[3] += 1
                        elif 400 <= len(line) < 450:
                            reads_length[4] += 1
                        elif 450 <= len(line) < 500:
                            reads_length[5] += 1
                        elif 500 <= len(line) < 550:
                            reads_length[6] += 1
                        elif 550 <= len(line) < 600:
                            reads_length[7] += 1

    all_reads = np.sum(reads_length)
    reads_ratio = np.zeros(reads_length.shape)
    for i in range(len(reads_length)):
        reads_ratio[i] = reads_length[i] / all_reads * 100

    draw_bar(reads_ratio, output_path)


def draw_bar(reads_length, output_path):
    k = [x for x in range(250, 650, 50)]
    plt.bar(k, reads_length, align='center', color='blue', width=50.1)
    plt.xlabel('Sequence length(bp)', fontdict={'family': 'Times New Roman', 'size': 16})
    plt.ylabel('Percent of Sequences(%)', fontdict={'family': 'Times New Roman', 'size': 16})
    plt.title('Effective sequence length distribution', pad=10, fontdict={'family': 'Times New Roman', 'size': 22})
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'reads_length_bar.pdf'))


def main():
    paths = rc.get_items('path')
    compute_length_ratio(paths[0][1], paths[1][1])


if __name__ == '__main__':
    main()
