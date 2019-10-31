#!/usr/bin/env bash

##################
pipeline_path='/home/pipeline/New_keyan_report'
########start#######
python $pipeline_path/bin/1_glob_path_otu.py
python $pipeline_path/bin/2_convert_otu_m2_to_m3_silva_v2.py
python $pipeline_path/bin/3_otu_stats.py
python $pipeline_path/bin/4_glob_path_alpha.py
python $pipeline_path/bin/5_glob_path_picrust1.py
#python $pipeline_path/bin/6_alpha_stats.py
#python $pipeline_path/bin/7_heatmap_groups.py
python $pipeline_path/bin/8_fa_length_bar.py
python $pipeline_path/bin/9_rank_abundance.py
#python $pipeline_path/bin/10_stack_bar.py
#python $pipeline_path/bin/11_boxplot.py
#python $pipeline_path/bin/12_anosim.py
#python $pipeline_path/bin/13_PCA.py
#python $pipeline_path/bin/14_PCoA.py
#python $pipeline_path/bin/15_NMDS.py
#python $pipeline_path/bin/16_hierarchy_cluster.py
python $pipeline_path/bin/17_rarefaction_curve.py

