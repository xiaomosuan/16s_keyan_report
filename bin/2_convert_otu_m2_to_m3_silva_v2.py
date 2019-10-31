import re
import os
import numpy as np
import read_config as rc


def m2_to_m3(output_path):
    f = open(os.path.join(output_path, 'combined_otu_table_m3_std.txt'), 'w')
    f_ln = open(os.path.join(output_path, 'combined_otu_table_m3_std_ln.txt'), 'w')
    m3 = {}
    sample_list = []
    with open(os.path.join(output_path, 'combined_otu_table_m2_std.txt'), 'r') as file:
        for line in file:
            line = line.strip().split('\t')
            if len(re.findall('^OTU ID', line[0])) > 0:
                for i in range(1, len(line)-1):
                    sample_list.append(line[i])

                print('OTU ID\t{}\ttaxonomy'.format('\t'.join(sample_list)), file=f)
                print('OTU ID\t{}\ttaxonomy'.format('\t'.join(sample_list)), file=f_ln)
                continue

            if len(re.findall('^Unassigned$', line[-1])) > 0:
                continue

            taxs = line[-1].split('; ')
            for i in range(len(taxs)):
                if len(re.findall('__$', taxs[i])) > 0:
                    continue

                par1 = m3
                par1 = par1.setdefault(taxs[i], {})
                for j in range(1, len(line) - 1):
                    par2 = par1
                    par2.setdefault(sample_list[j-1], 0)
                    par2[sample_list[j-1]] += float(line[j])

    for k1, v1 in m3.items():
        print(k1, end='\t', file=f)
        print(k1, end='\t', file=f_ln)
        for i in range(len(sample_list)):
            v = m3[k1][sample_list[i]]
            v_ln = np.log(m3[k1][sample_list[i]] + 1)
            print(v, end='\t', file=f)
            print(v_ln, end='\t', file=f_ln)

        print(k1, file=f)
        print(k1, file=f_ln)

    f.close()
    f_ln.close()


def main():
    paths = rc.get_items('path')
    m2_to_m3(paths[1][1])


if __name__ == '__main__':
    main()
