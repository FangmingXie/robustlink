#!/usr/bin/env python
# coding: utf-8

import sys
from multiprocessing import Pool,cpu_count
from scipy import sparse
import pickle
import argparse
import logging
import pandas as pd
import numpy as np

sys.path.insert(0, "../")
import utils
import enhancer_gene_utils

logger = utils.create_logger()

def pipe_corr_analysis_atac(
        common_modx_cells, common_mody_cells,
        cell_cell_knn_xaxis, cell_cell_knn_yaxis,
        common_genes, common_enhancer_regions,
        X, Y, 
        modx_clsts, knn_xy, 
        enhancer_gene_to_eval,
        output_corrs,
        corr_type='pearsonr',
        force=False,
        num_metacell_limit=0,
    ):
    """
    """
    # new cells  
    common_modx_cells_updated = np.intersect1d(common_modx_cells, cell_cell_knn_xaxis)
    common_mody_cells_updated = np.intersect1d(common_mody_cells, cell_cell_knn_yaxis)

    # make sure the original matrices have the correct index
    x_idx = utils.get_index_from_array(common_modx_cells, common_modx_cells_updated)
    y_idx = utils.get_index_from_array(common_mody_cells, common_mody_cells_updated)
    X = X.tocsc()[:, x_idx] 
    Y = Y.tocsc()[:, y_idx]

    # make sure knn_xy, knn_xx have the right cell index
    cell_idx_xaxis = utils.get_index_from_array(cell_cell_knn_xaxis, common_modx_cells_updated)
    cell_idx_yaxis = utils.get_index_from_array(cell_cell_knn_yaxis, common_mody_cells_updated)
    knn_xy = knn_xy.tocsr()[cell_idx_xaxis,:].tocsc()[:,cell_idx_yaxis] # x-by-y
    modx_clsts = modx_clsts.reindex(common_modx_cells_updated)

    logging.info("{}_{}_{}_{}".format(knn_xy.shape, modx_clsts.shape, X.shape, Y.shape,))

    for clst_col in modx_clsts.columns: 
        logging.info(clst_col)
        output_corr = output_corrs.format(clst_col)
        if not force and os.path.isfile(output_corr):
            logging.info("skip {}, already exists...".format(output_corr))
            continue # skip the existing file

        # choose one clustering to proceed
        uniq_labels = np.sort(modx_clsts[clst_col].unique()) 
        logging.info("Number of metacells: {}".format(len(uniq_labels)))
        if num_metacell_limit > 0 and len(uniq_labels) > num_metacell_limit:
            logging.info("skip {}, exceeding max num_metacell_limit...".format(len(uniq_labels)))
            continue

        knn_xz = enhancer_gene_utils.turn_cluster_labels_to_knn(modx_clsts[clst_col].values, 
                                            uniq_labels,
                                           )

        # Dec 21,2020
        # gene by metacell (counts)
        gc_rna = X.dot(knn_xz).todense() 
        # normalization (logCPM)
        gc_rna = utils.logcpm(pd.DataFrame(gc_rna)).values

        # Dec 21,2020
        enh_lengths = pd.Series((common_enhancer_regions['end']-common_enhancer_regions['start']).values)
        # enhancer by metacell (counts)
        knn_yz = knn_xy.T.dot(knn_xz)
        ec_atac = Y.dot(knn_yz).todense() 
        # normalization (logTPM)
        ec_atac = utils.logtpm(pd.DataFrame(ec_atac), enh_lengths).values
        logging.info("{} {}".format(gc_rna.shape, ec_atac.shape,))

        # corr analysis
        output_corr = output_corrs.format(clst_col)
        (to_correlate, corrs, corrs_shuffled, corrs_shuffled_cells) = enhancer_gene_utils.compute_enh_gene_corrs(
            gc_rna, ec_atac, 
            common_genes, np.arange(len(ec_atac)),
            enhancer_gene_to_eval['gene'].values, 
            enhancer_gene_to_eval['ens'].values, 
            output_file=output_corr, corr_type=corr_type, chunksize=100000, verbose_level=0,
            )
    return

def wrap_corr_analysis_atac(
        mod_x, mod_y, 
        input_name_tag, i_sub,
        corr_type='pearsonr',
        force=False,
        num_metacell_limit=0,
    ):
    """
    """
    # (i, k, --r)
    output_corrs = './results/{}_{}_{{}}_{}_corrs.pkl'.format(input_name_tag, i_sub, corr_type)

    # input enh-gene tables, gene-by-cell, enhancer-by-cell matrices
    input_enh_gene_table = './data/counts/enhancer_gene_pairs_1mbp.tsv' 
    input_bundle_dirc = './data/counts'
    bundle_fnames = (
        'cell_rna.pkl',
        'cell_atac.pkl',

        'gene_rna.pkl',
        'enh.pkl',

        'mat_rna.pkl',
        'mat_atac.pkl',
    )


    # for knn_xx
    input_knn_dirc = './results'
    input_modx_clsts = [
        'clusterings_{}_{}_sub{}.tsv.gz'.format(mod_x, input_name_tag, i_sub),
    ]

    # for knn_xy
    input_knn_xy = 'knn_across_{}_{}_{}.npz.{}.npz'.format(input_name_tag, mod_x, mod_y, i_sub) 
    input_knn_cells_xaxis = 'cells_{}_{}.npy.{}.npy'.format(mod_x, input_name_tag, i_sub)
    input_knn_cells_yaxis = 'cells_{}_{}.npy.{}.npy'.format(mod_y, input_name_tag, i_sub)

    # # Load data 
    # input_bundle
    with utils.cd(input_bundle_dirc):
        bundle = []
        for fname in bundle_fnames:
            #  save all as pickle file
            with open(fname, "rb") as fh:
                item = pickle.load(fh)
            bundle.append(item)
            logging.info("{}_{}_{}".format(type(item), item.shape, fname))

    (common_modx_cells, common_mody_cells, 
     common_genes, common_enhancer_regions,
     X, Y, 
    ) = bundle

    # input knn networks 
    with utils.cd(input_knn_dirc):
        # for knn_xx 
        # modx_clsts = pd.read_csv(input_modx_clsts, sep='\t',index_col=0)
        modx_clsts = pd.concat([
            pd.read_csv(fname, sep='\t',index_col=0)
            for fname in input_modx_clsts
        ], axis=1)
        # for knn_xy 
        knn_xy = sparse.load_npz(input_knn_xy)  
        cell_cell_knn_xaxis = np.load(input_knn_cells_xaxis, allow_pickle=True)
        cell_cell_knn_yaxis = np.load(input_knn_cells_yaxis, allow_pickle=True)

        logging.info("{} {} {} {}".format(modx_clsts.shape, knn_xy.shape, 
                                          cell_cell_knn_xaxis.shape, 
                                          cell_cell_knn_yaxis.shape,
                                         )
                    )

    # enhancer-gene linkage
    enhancer_gene_to_eval = pd.read_csv(input_enh_gene_table, sep='\t')

    pipe_corr_analysis_atac(
        common_modx_cells, common_mody_cells,
        cell_cell_knn_xaxis, cell_cell_knn_yaxis,
        common_genes, common_enhancer_regions,
        X, Y, 
        modx_clsts, knn_xy, 
        enhancer_gene_to_eval,
        output_corrs,
        corr_type=corr_type,
        force=force,
        num_metacell_limit=num_metacell_limit,
    )
    return 

def create_parser():
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-modx', '--mod_x', required=True, type=str)
    parser.add_argument('-mody', '--mod_y', required=True, type=str)
    parser.add_argument('-tag', '--input_name_tag', required=True, type=str)
    parser.add_argument('-isub', '--i_sub', type=str, help="[0-9]")
    parser.add_argument('-ct', '--corr_type', choices=['pearsonr', 'spearmanr'],  
                                              default='pearsonr',
                                              help="choose from pearsonr or spearmanr")
    parser.add_argument('-f', '--force', action='store_true', help='to overwrite existing file')
    parser.add_argument('-n', '--num_metacell_limit', default=0, type=int, help='max num of metacells')
    return parser


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    parser = create_parser()
    args = parser.parse_args()

    # output setting
    # run this with each combination of (i_sub, knn)
    mod_x = args.mod_x
    mod_y = args.mod_y
    input_name_tag = args.input_name_tag
    i_sub = args.i_sub
    corr_type = args.corr_type
    force = args.force
    num_metacell_limit = args.num_metacell_limit

    logging.info(",".join([mod_x, mod_y, input_name_tag, str(i_sub), corr_type, str(force), str(num_metacell_limit)]))

    wrap_corr_analysis_atac(
        mod_x, mod_y, 
        input_name_tag, i_sub,
        corr_type=corr_type,
        force=force,
        num_metacell_limit=num_metacell_limit,
    )

