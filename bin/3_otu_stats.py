import re
import os
import read_config as rc


name_dict = {
    'k': '界(Kingdom)',
    'p': '门(Phylum)',
    'c': '纲(Class)',
    'o': '目(Order)',
    'f': '科(Family)',
    'g': '属(Genus)',
    's': '种(Species)',
}


def get_tax_num(output_file):
    tax_type = ['k', 'p', 'c', 'o', 'f', 'g', 's']
    sample_list, tax_dict = [], {}
    with open(os.path.join(output_file, 'combined_otu_table_m3_std.txt'), 'r') as f:
        for line in f:
            aa = line.strip().split('\t')
            if len(re.findall('^OTU ID', aa[0])) > 0:
                for i in range(1, len(aa) - 1):
                    sample_list.append(aa[i])
                continue

            for tax in tax_type:
                if re.match('^{}__'.format(tax), aa[0]):
                    par = tax_dict
                    par = par.setdefault('{}'.format(tax), {})
                    for i in range(1, len(aa) - 1):
                        par.setdefault(sample_list[i - 1], 0)
                        if float(aa[i]) != 0:
                            par[sample_list[i - 1]] += 1

    f_out = open(os.path.join(output_file, 'taxonomy_stats.txt'), 'w')
    print('Taxonomy\t{}'.format('\t'.join(sample_list)), file=f_out)
    for tax in tax_type:
        print(name_dict[tax], end='\t', file=f_out)

        for i in range(len(sample_list)):
            if sample_list[i] in tax_dict[tax]:
                print(tax_dict[tax][sample_list[i]], end='', file=f_out)
            else:
                print(0, end='', file=f_out)

            if i < len(sample_list) - 1:
                print(end='\t', file=f_out)
            else:
                print(end='\n', file=f_out)

    f_out.close()


def main():
    paths = rc.get_items('path')
    get_tax_num(paths[1][1])


if __name__ == '__main__':
    main()
