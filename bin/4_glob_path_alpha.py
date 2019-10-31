import os
import re
import read_config as rc
"""
合并目录下样品的alpha表
input_path/sample/03_Alpha/alpha_qiime_silva99v2_m1.txt
"""


def glob_all_otu(input_path, output_path):
    f_out = open(os.path.join(output_path, 'combined_alpha_m1.txt'), 'w')
    print('Sample\tchao1\tobserved_species\tshannon\tsimpson\tgoods_coverage', file=f_out)
    for roots, paths, files in os.walk(input_path):
        if len(re.findall('03_Alpha$', roots)) > 0:
            for file in files:
                if file == 'alpha_qiime_silva99v2_m1.txt':
                    with open(os.path.join(roots, file), 'r') as f:
                        for line in f:
                            if len(re.findall('chao1', line)) > 0:
                                continue

                            line = line.strip()
                            print(line, file=f_out)

    f_out.close()


def main():
    paths = rc.get_items('path')
    glob_all_otu(paths[0][1], paths[1][1])


if __name__ == '__main__':
    main()
