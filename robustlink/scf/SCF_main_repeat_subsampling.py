#!/usr/bin/env python3
"""SingleCellFusion main rontine"""

from .__init__ import *

import collections
import argparse
import anndata

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

    optional.add_argument("-tag", "--nametag", help='output files will contain this tag', type=str, default='scfusion')
    optional.add_argument("--ka_smooth", type=int, default=30)
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

def get_default_mod_direction(mod):
    """
    """
    if mod == 'mc':
        mod_direction = -1
    elif mod == 'rna':
        mod_direction = 1
    elif mod == 'atac':
        mod_direction = 1
    else:
        raise ValueError("choose from ['mc', 'rna', 'atac']")
    return mod_direction 

def get_default_mod_normalization(mod):
    """
    """
    if mod == 'mc':
        mod_norm = 'mc'
    elif mod == 'rna':
        mod_norm = 'cpm'
    elif mod == 'atac':
        mod_norm = 'tpm'
    else:
        raise ValueError("choose from ['mc', 'rna', 'atac']")
    return mod_norm 

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

    subsample_fraction = args.subsample_fraction
    subsample_times = args.subsample_times

    data_dir = args.data_dir 
    outdir   = args.outdir  
    ka_smooth = args.ka_smooth #
    knn = args.knn

    # # Configs  
    name = args.nametag

    mods_selected = [dtst[:-len('.h5ad')] for dtst in args.input_datasets if dtst.endswith('.h5ad')]
    assert len(mods_selected) == len(args.input_datasets)
    features_selected = [dtst[:-len('.h5ad')] for dtst in args.feature_datasets if dtst.endswith('.h5ad')]
    assert len(features_selected) == len(args.feature_datasets)
    # check features
    for features_modality in features_selected:
        assert (features_modality in mods_selected)
    input_modalities = args.input_modalities
    
    # settings
    Mod_info = collections.namedtuple('Mod_info', [
        'mod', 
        'mod_category',
        'mod_direction', # +1 or -1
    ])
    settings = collections.OrderedDict({})
    for mod, mod_cat in zip(mods_selected, input_modalities):
        settings[mod] = Mod_info(
            mod, 
            mod_cat,
            get_default_mod_direction(mod_cat),
            )

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

    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # end of configurations

    # ## Read in data 
    logging.info('* Begin integration')

    metas = collections.OrderedDict()
    gxc_hvftrs = collections.OrderedDict()
    for mod in mods_selected:
        # read in the h5ad file
        adata = anndata.read(os.path.join(data_dir, mod+'.h5ad'))
        # metadata (cell level)
        metas[mod] = adata.var.copy() 
        logging.info("Metadata {} {}".format(mod, metas[mod].shape))

        # the feature matrix
        if settings[mod].mod_category == 'mc':
            gxc_hvftrs[mod] = pd.DataFrame(adata.X, index=adata.obs.index.values, columns=adata.var.index.values)
            logging.info("Feature matrix {} {}".format(mod, gxc_hvftrs[mod].shape))
            continue

        gxc_hvftrs[mod] = basic_utils.gc_matrix_from_anndata(adata)
        logging.info("Feature matrix {} {}".format(mod, gxc_hvftrs[mod].data.shape))
    logging.info('Done reading data')

    # subsampling a given fraction of cells from each modality
    for i in range(subsample_times):
        # # specifying output files (one set per subsampling)
        output_pcX_all    = os.path.join(outdir, f'{name}_s{i}_fused_pc'              + f'.h5ad') #  

        # required for downsamp (8/7/2020)
        output_cells      = os.path.join(outdir, f'{name}_s{i}_cells_{{}}'            + f'.txt') # mod 
        # new required arguments (7/27/2020)
        save_knn = True  
        output_knn_within = os.path.join(outdir, f'{name}_s{i}_knn_within_{{}}'       + f'.npz') # mod 
        output_knn_across = os.path.join(outdir, f'{name}_s{i}_knn_across_{{}}_{{}}'  + f'.npz') # modx, mody 
        # # output files end

        logging.info("Subsampling {}/{}".format(i+1, subsample_times))
        metas_sub, gxc_hvftrs_sub = subsampling(
                    mods_selected, settings, metas, gxc_hvftrs, p=subsample_fraction,
                    )
        for mod in mods_selected:
            fout = output_cells.format(mod)
            cells_mod = metas_sub[mod].index.values
            np.savetxt(fout, cells_mod, fmt='%s')

        # ## run SCF
        pcX_all, cells_all = SCF_utils.core_scf_routine(mods_selected, features_selected, settings, 
                                                        metas_sub, gxc_hvftrs_sub, # sub 
                                                        ps, drop_npcs,
                                                        cross_mod_distance_measure, knn, relaxation, n_cca,
                                                        npc,
                                                        output_pcX_all, 
                                                        ka_smooth=ka_smooth,
                                                        save_knn=save_knn,
                                                        output_knn_within=output_knn_within,
                                                        output_knn_across=output_knn_across,
                                                        )
        logging.info('Done integration into a common PC space')
    return 

if __name__ == '__main__':
    log = basic_utils.create_logger()
    parser = argparse.ArgumentParser()
    add_args(parser)
    args = parser.parse_args()
    main(args)
