import os
import re
import read_config as rc
"""
合并目录下样品的OTU表
input_path/sample/01_OTU/qiime_silva99v2/otu_table_m2_std.txt
"""


def glob_all_otu(input_path):
    all_dict, sample_list = {}, []
    for roots, paths, files in os.walk(input_path):
        if len(re.findall('qiime_silva99v2$', roots)) > 0:
            for file in files:
                if file == 'otu_table_m2_std.txt':
                    with open(os.path.join(roots, file), 'r') as f:
                        for line in f:
                            if len(re.findall('^OTU ID', line)) > 0:
                                sample_id = re.findall('^OTU ID\t(.*)\ttaxonomy', line)[0]
                                sample_list.append(sample_id)
                                continue

                            aa = line.strip().split('\t')
                            par1 = all_dict
                            par1 = par1.setdefault(aa[-1], {})
                            par1[sample_id] = aa[1]

    sort_sample = sorted(sample_list)
    sort_dict = sorted(all_dict.items(), key=lambda x: x[0])
    return sort_sample, sort_dict


def combined_otu(sort_sample, sort_dict, output_path):
    f = open(os.path.join(output_path, 'combined_otu_table_m2_std.txt'), 'w')
    print('OTU ID\t{}\ttaxonomy'.format('\t'.join(sort_sample)), file=f)
    num = 0
    for k, v in sort_dict:
        num += 1
        print('OTU{}'.format(num), end='\t', file=f)
        for i in range(len(sort_sample)):
            if sort_sample[i] in v:
                print(v[sort_sample[i]], end='\t', file=f)
            else:
                print(0, end='\t', file=f)

        print(k, file=f)

    f.close()


def main():
    paths = rc.get_items('path')
    sort_sample, sort_dict = glob_all_otu(paths[0][1])
    combined_otu(sort_sample, sort_dict, paths[1][1])


if __name__ == '__main__':
    main()
