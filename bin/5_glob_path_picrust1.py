import os
import re
import read_config as rc
"""
合并目录下样品的OTU表
input_path/sample/02_Picrust1/kegg_predicted_metagenomes.L3.txt
input_path/sample/02_Picrust1/cog_predicted_metagenomes.L3.txt
"""


def glob_all_otu(input_path):
    all_dict_kegg, tax_kegg, all_dict_cog, tax_cog, sample_list = {}, {}, {}, {}, []
    for roots, paths, files in os.walk(input_path):
        if len(re.findall('02_Picrust1$', roots)) > 0:
            for file in files:
                if file == 'kegg_predicted_metagenomes.L3.txt':
                    with open(os.path.join(roots, file), 'r') as f:
                        for line in f:
                            if re.match('# Constructed from biom file', line):
                                continue

                            if len(re.findall('^(?:#|)OTU ID', line)) > 0:
                                sample_id = re.findall('^(?:#|)OTU ID\t(.*)\tKEGG_Pathways', line)[0]
                                sample_list.append(sample_id)
                                continue

                            aa = line.strip().split('\t')
                            par1 = all_dict_kegg
                            par1 = par1.setdefault(aa[0], {})
                            tax_kegg[aa[0]] = aa[-1]
                            par1[sample_id] = aa[1]

                if file == 'cog_predicted_metagenomes.L2.txt':
                    with open(os.path.join(roots, file), 'r') as f:
                        for line in f:
                            if re.match('# Constructed from biom file', line):
                                continue

                            if len(re.findall('^(?:#|)OTU ID', line)) > 0:
                                sample_id = re.findall('^(?:#|)OTU ID\t(.*)\tCOG_Category', line)[0]
                                sample_list.append(sample_id)
                                continue

                            aa = line.strip().split('\t')
                            par2 = all_dict_cog
                            par2 = par2.setdefault(aa[0], {})
                            tax_cog[aa[0]] = aa[-1]
                            par2[sample_id] = aa[1]

    sort_sample = sorted(sample_list)
    sort_dict_kegg = sorted(all_dict_kegg.items(), key=lambda x: x[0])
    sort_dict_cog = sorted(all_dict_cog.items(), key=lambda x: x[0])
    return sort_sample, sort_dict_kegg, tax_kegg, sort_dict_cog, tax_cog


def combined_otu(sort_sample, sort_dict_kegg, tax_kegg, sort_dict_cog, tax_cog, output_path):
    f = open(os.path.join(output_path, 'combined_kegg.txt'), 'w')
    print('OTU ID\t{}\ttaxonomy'.format('\t'.join(sort_sample)), file=f)
    for k, v in sort_dict_kegg:
        print(k, end='\t', file=f)
        for i in range(len(sort_sample)):
            if sort_sample[i] in v:
                print(v[sort_sample[i]], end='\t', file=f)
            else:
                print(0, end='\t', file=f)

        print(tax_kegg[k], file=f)

    f.close()

    f = open(os.path.join(output_path, 'combined_cog.txt'), 'w')
    print('OTU ID\t{}\ttaxonomy'.format('\t'.join(sort_sample)), file=f)
    for k, v in sort_dict_cog:
        print(k, end='\t', file=f)
        for i in range(len(sort_sample)):
            if sort_sample[i] in v:
                print(v[sort_sample[i]], end='\t', file=f)
            else:
                print(0, end='\t', file=f)

        print(tax_cog[k], file=f)

    f.close()


def main():
    paths = rc.get_items('path')
    sort_sample, sort_dict_kegg, tax_kegg, sort_dict_cog, tax_cog = glob_all_otu(paths[0][1])
    combined_otu(sort_sample, sort_dict_kegg, tax_kegg, sort_dict_cog, tax_cog, paths[1][1])


if __name__ == '__main__':
    main()
