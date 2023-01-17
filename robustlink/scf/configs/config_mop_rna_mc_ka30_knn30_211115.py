#!/usr/bin/env python3
"""An example configuration file
"""
import sys
import os
import collections

# Assuming the cell order in the metadata tables are the same as those in the gene level matrices
# The output knn matrices follow such order as well 

dir_path = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(dir_path, '../../data')
outdir   = os.path.join(dir_path, '../../results')
ka_smooth = 30 
knn = 30
date = 211115

# # fixed dataset configs
# sys.path.insert(0, DATA_DIR)
# from __init__datasets import *

# # Configs  
name = 'mop_rna_mc_ka{}_knn{}_{}'.format(ka_smooth, knn, date,)
output_pcX_all = outdir + '/pcX_all_{}.npy'.format(name)
output_cells_all = outdir + '/cells_all_{}.npy'.format(name) 
output_imputed_data_format = outdir + '/imputed_data_{}_{{}}.npy'.format(name)
output_clst_and_umap = outdir + '/intg_summary_{}.tsv'.format(name)
output_figures = outdir + '/figures/{}_{{}}.{{}}'.format(name)
output_cluster_centroids = outdir + '/centroids_{}.pkl'.format(name)

save_knn = True # new required arguments (7/27/2020) 
output_knn_within = outdir + "/knn_within_{}_{{}}.npz".format(name)
output_knn_across = outdir + "/knn_across_{}_{{}}_{{}}.npz".format(name)
# end of new required arguments (7/27/2020)

# required for downsamp (8/7/2020)
output_cells = outdir + "/cells_{{}}_{}.npy".format(name)

meta_f = os.path.join(data_dir, '{0}_metadata.tsv')
hvftrs_f = os.path.join(data_dir, '{0}_hvfeatures.{1}')
hvftrs_gene = os.path.join(data_dir, '{0}_hvfeatures.gene')
hvftrs_cell = os.path.join(data_dir, '{0}_hvfeatures.cell')

mods_selected = [
    'mc',
    'rna',
    ]

features_selected = ['mc']
# check features
for features_modality in features_selected:
    assert (features_modality in mods_selected)

# within modality
ps = {'mc': 0.9,
      'atac': 0.1,
      'rna': 0.7,
     }
drop_npcs = {
      'mc': 0,
      'atac': 0,
      'rna': 0,
     }
ka_smooth = ka_smooth # default: 5

# across modality
cross_mod_distance_measure = 'correlation' # cca
knn = knn 
relaxation = 3
n_cca = 30

# PCA
npc = 50

# clustering
k = 30 # default: 30
resolutions = [0.1, 1, 2, 4]
# umap
umap_neighbors = 60
min_dist = 0.5

# meta settings
mods = (
    'mc',
    'atac',
    'rna',
)

Mod_info = collections.namedtuple('Mod_info', [
    'mod', 
    # 'name',
    'mod_category',
    'norm_option',
    'mod_direction', # +1 or -1
    'cell_col',
    # 'cluster_col', 
    # 'annot_col', 
    # 'category_col', # neuron or not
    'global_mean', # in general or mch
    'global_mean_mcg', # only for mcg
    # 'total_reads',
    # 'color',
    # 'species',
])

# settngs
settings_mc = Mod_info(
    mods[0],
    # 'DNA methylation',
    'mc',
    'mc', 
    -1, # negative direction 
    'cell',
    # 'SubCluster',
    # 'SubCluster', #'major_clusters' 'sub_cluster' 
    # 'SubCluster',
    'CH_Rate',
    'CG_Rate',
    # 'FinalReads',
    # 'C5',
    # 'mouse',
)

settings_atac = Mod_info(
    mods[1],
    # 'ATAC Seq',
    'atac',
    'tpm', 
    +1, # direction 
    'cell',
    # 'cluster',
    # 'cluster',
    # 'cluster',
    '',
    '',
    # '',
    # 'C2',
    # 'mouse',
)

settings_rna = Mod_info(
    mods[2],
    # '10X V3 cells',
    'rna',
    'cpm', 
    +1, # direction 
    'cell',
    # 'cluster_id',
    # 'cluster_label',
    # 'class_label',
    '',
    '',
    # '',
    # 'C6',
    # 'mouse',
)

settings = collections.OrderedDict({
    mods[0]: settings_mc,
    mods[1]: settings_atac,
    mods[2]: settings_rna,
})

