#!/usr/bin/env python3
"""SingleCellFusion main rontine"""

from .__init__ import *

import collections
import itertools
import sys
import pickle
import argparse

from . import basic_utils
from . import SCF_utils

def add_args(parser):
    """
    """
    # parser.add_argument("-c", "--config_py", help="Configuration file", required=True)

    required = parser.add_argument_group('required')
    optional = parser.add_argument_group('optional') 
    advanced = parser.add_argument_group('advanced')

    required.add_argument("-i", "--data_dir", help='data dir', required=True)
    required.add_argument("-o", "--outdir", help='outdir', required=True)
    required.add_argument(
        "-id", "--input_datasets", 
        type=str,
        nargs="+",
        required=True,
        help='''(list of str)
             list of .h5ad files (AnnData) in the data_dir 
             '''
    )
    required.add_argument(
        "-im", "--input_modalities", 
        type=str,
        nargs="+",
        required=True,
        help='''(list of str)
             Data modalities chosen from 'rna', 'atac', or 'mc'. This should be 
             listed in the same order as input_datasets. 
             '''
    )
    required.add_argument(
        "-fd", "--feature_datasets", 
        type=str,
        nargs="+",
        required=True,
        help='''(list of str)
             Dataset(s) whose features all other datasets will impute into.
             This should be a subset of --input_datasets.
             Enter multiple datasets as a space-separated list of filenames.
             '''
    )

    optional.add_argument("-on", "--outname", help='output files will contain this name', type=str, default='scf')
    optional.add_argument("--ka_smooth", default=30)
    optional.add_argument("--knn", help=">1", type=int, default=30)
    optional.add_argument("-s", "--subsample_fraction", help="0~1", type=float, default=1)
    optional.add_argument("-sn", "--subsample_times", help=">1", type=int, default=1)
    # constraint kNN across modalities
    optional.add_argument(
        "--relaxation", 
        type=float,
        default=3,
        help='''(float)
             A value between 1 to infinity. 
             This is a parameter that constraints the number of neighbors a cell is allowed to receive.
             Assume dataset 1 has N1 cells, dataset 2 has N2 cells. To find k neighbors in dataset 2 for 
             every cell in dataset 1 means on average each cell in dataset 2 receives (kN1/N2) connections.
             However, not all cells in dataset 2 gets the same number of connections. We therefore set an 
             upper bound for the number of connections a cell in dataset 2 can receive to be:
                (kN1/N2)*relaxation
             where relaxation >= 1. Relaxation=1 enforces a hard limit that every cell receives 
             the same number of nearest neighbors, while relaxation=infinity approaches traditional kNN.
             '''
    )

    advanced.add_argument("--drop_npcs",     type=int, default=0)
    advanced.add_argument(
        "--smoothing_fractions", 
        nargs="+",
        type=float,
        default=[0.7, 0.1, 0.9], # rna, atac, mc
        help='''(list of floats) 
             A list of three values between 0 to 1 that controls the relative contribution
             from the cell itself vs. its neighbors in within-dataset smoothing, 
             specified for 'rna', 'atac', 'mc' data, respectively.
             '''
    )
    advanced.add_argument(
        "--num_pcs", 
        type=int,
        default=50,
        help='''(integer)
             Number of Principal Components to keep for each dataset 
             for smoothing and for clustering/embedding after imputation.
             '''
    )

    return parser

def parse_filename(data_file):
    """turn a xxx/xxx/XXXX.h5ad into XXXX 
    """
    dataset_name = os.path.basename(data_file)
    if dataset_name.endswith('.h5ad'):
        dataset_name = dataset_name[:-len('.h5ad')]
    else:
        raise ValueError("filenames don't have the format xxxx.h5ad")
    return dataset_name

def modality_default_options(mod):
    """
    """
    if mod == 'mc':
        mod_direction = -1
        # norm_option = 'mc'
    elif mod == 'rna':
        mod_direction = 1
        # norm_option = 'cpm'
    elif mod == 'atac':
        mod_direction = 1
        # norm_option = 'tpm'
    else:
        raise ValueError("choose from ['mc', 'rna', 'atac']")
    return mod_direction 

def subsampling(mods_selected, settings, metas, gxc_hvftrs, p=0, n=0):
    """Do many datasets at the same time
    p - fraction of cells from each dataset to be included
    """
    if p < 1:
        metas_sub = collections.OrderedDict()
        gxc_hvftrs_sub = collections.OrderedDict()
        for mod in mods_selected: 
            # subsample meta
            if p != 0: 
                cells_included = metas[mod].index.values[np.random.rand(len(metas[mod]))<p]
            elif n != 0:
                if n > len(metas[mod]):
                    n = len(metas[mod])
                    
                cells_included = metas[mod].sample(n=n).index.values
            metas_sub[mod] = metas[mod].loc[cells_included]

            # subsample gxc_hvftrs
            if settings[mod].mod_category == 'mc':
                gxc_hvftrs_sub[mod] = gxc_hvftrs[mod][cells_included]
                logging.info("{} {} {}".format(mod, metas_sub[mod].shape, gxc_hvftrs_sub[mod].shape))
                continue

            cells_included_idx = basic_utils.get_index_from_array(gxc_hvftrs[mod].cell, cells_included)
            gxc_hvftrs_sub[mod] = GC_matrix(
                                            gxc_hvftrs[mod].gene,
                                            cells_included,
                                            gxc_hvftrs[mod].data.tocsc()[:, cells_included_idx],
                                            )
            logging.info("{} {} {}".format(mod, metas_sub[mod].shape, gxc_hvftrs_sub[mod].data.shape))
        return metas_sub, gxc_hvftrs_sub
    else:
        return metas, gxc_hvftrs


def main(args):

    # normal
    subsample_fraction = args.subsample_fraction
    subsample_times = args.subsample_times

    # to fix
    # dir_path = os.path.dirname(os.path.realpath(__file__))
    data_dir = args.data_dir # os.path.join(dir_path, '../../data')
    outdir   = args.outdir # os.path.join(dir_path, '../../results')
    ka_smooth = args.ka_smooth # 
    knn = args.knn

    # # Configs  
    name = args.outname
    output_pcX_all = outdir + '/pcX_all_{}.npy'.format(name)
    output_cells_all = outdir + '/cells_all_{}.npy'.format(name) 
    output_imputed_data_format = outdir + '/imputed_data_{}_{{}}.npy'.format(name)

    # new required arguments (7/27/2020)
    save_knn = True  
    output_knn_within = outdir + "/knn_within_{}_{{}}.npz".format(name)
    output_knn_across = outdir + "/knn_across_{}_{{}}_{{}}.npz".format(name)
    # end of new required arguments (7/27/2020)

    # required for downsamp (8/7/2020)
    output_cells = outdir + "/cells_{{}}_{}.npy".format(name)

    meta_f      = os.path.join(data_dir, '{0}_metadata.tsv')
    hvftrs_f    = os.path.join(data_dir, '{0}_hvfeatures.{1}')
    hvftrs_gene = os.path.join(data_dir, '{0}_hvfeatures.gene')
    hvftrs_cell = os.path.join(data_dir, '{0}_hvfeatures.cell')

    # mods_selected = [
    #     'mc',
    #     'rna',
    #     ]
    # features_selected = ['mc']
    mods_selected = args.input_datasets
    features_selected = args.feature_datasets

    # check features
    for features_modality in features_selected:
        assert (features_modality in mods_selected)

    # within modality
    ps = {mod: frac for mod, frac in zip(['rna', 'atac', 'mc'], args.smoothing_fractions)}
    drop_npcs = {
        'rna':  args.drop_npcs,
        'atac': args.drop_npcs,
        'mc':   args.drop_npcs,
        }

    # across modality
    cross_mod_distance_measure = 'correlation'
    knn = args.knn 
    relaxation = args.relaxation
    n_cca = -1 # 30
    # PCA
    npc = args.num_pcs

    ## dataset
    # meta settings
    mods = (
        'mc',
        'atac',
        'rna',
    )
    Mod_info = collections.namedtuple('Mod_info', [
        'mod', 
        'mod_category',
        'norm_option',
        'mod_direction', # +1 or -1
        'cell_col',
        'global_mean', # in general or mch
        'global_mean_mcg', # only for mcg
    ])
    # settngs
    settings_mc = Mod_info(
        mods[0],
        'mc',
        'mc', 
        -1, # negative direction 
        'cell',
        'CH_Rate',
        'CG_Rate',
    )
    settings_atac = Mod_info(
        mods[1],
        'atac',
        'tpm', 
        +1, # direction 
        'cell',
        '',
        '',
    )
    settings_rna = Mod_info(
        mods[2],
        'rna',
        'cpm', 
        +1, # direction 
        'cell',
        '',
        '',
    )
    settings = collections.OrderedDict({
        mods[0]: settings_mc,
        mods[1]: settings_atac,
        mods[2]: settings_rna,
    })

    ####!!actually parse things -- instead of import config.py
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # end of configurations

    ### ---- fixed after ----
    # ## Read in data 
    logging.info('* Begin integration')

    metas = collections.OrderedDict()
    for mod in mods_selected:
        metas[mod] = pd.read_csv(meta_f.format(mod), sep="\t").reset_index().set_index(settings[mod].cell_col)
        logging.info("Metadata {} {}".format(mod, metas[mod].shape))

    gxc_hvftrs = collections.OrderedDict()
    for mod in mods_selected:
        if settings[mod].mod_category == 'mc':
            f_mat = hvftrs_f.format(mod, 'tsv')
            gxc_hvftrs[mod] = pd.read_csv(f_mat, sep='\t', header=0, index_col=0) 
            logging.info("Feature matrix {} {}".format(mod, gxc_hvftrs[mod].shape))
            assert np.all(gxc_hvftrs[mod].columns.values == metas[mod].index.values) # make sure cell name is in the sanme order as metas (important if save knn mat)
            continue
            
        f_mat = hvftrs_f.format(mod, 'npz')
        f_gene = hvftrs_gene.format(mod)
        f_cell = hvftrs_cell.format(mod)
        _gxc_tmp = basic_utils.load_gc_matrix(f_gene, f_cell, f_mat)
        _gene = _gxc_tmp.gene
        _cell = _gxc_tmp.cell
        _mat = _gxc_tmp.data

        gxc_hvftrs[mod] = GC_matrix(_gene, _cell, _mat)
        assert np.all(gxc_hvftrs[mod].cell == metas[mod].index.values) # make sure cell name is in the sanme order as metas (important if save knn mat)
        logging.info("Feature matrix {} {}".format(mod, gxc_hvftrs[mod].data.shape))
    logging.info('Done reading data')


    # subsampling 80% of cells from each modality
    for i in range(subsample_times):
        output_tag = ".{}".format(i)

        logging.info("Subsampling {}/{}".format(i+1, subsample_times))
        metas_sub, gxc_hvftrs_sub = subsampling(
                    mods_selected, settings, metas, gxc_hvftrs, p=subsample_fraction,
                    )
        for mod in mods_selected:
            fout = (output_cells+output_tag).format(mod)
            cells_mod = metas_sub[mod].index.values
            np.save(fout, cells_mod)

        # ## run SCF
        pcX_all, cells_all = SCF_utils.core_scf_routine(mods_selected, features_selected, settings, 
                                                        metas_sub, gxc_hvftrs_sub, # sub 
                                                        ps, drop_npcs,
                                                        cross_mod_distance_measure, knn, relaxation, n_cca,
                                                        npc,
                                                        output_pcX_all+output_tag, output_cells_all+output_tag,
                                                        output_imputed_data_format+output_tag,
                                                        ka_smooth=ka_smooth,
                                                        save_knn=save_knn,
                                                        output_knn_within=output_knn_within+output_tag,
                                                        output_knn_across=output_knn_across+output_tag,
                                                        )
        logging.info('Done integration into a common PC space')

    return 

if __name__ == '__main__':
    log = basic_utils.create_logger()
    parser = argparse.ArgumentParser()
    add_args(parser)
    args = parser.parse_args()
    main(args)
