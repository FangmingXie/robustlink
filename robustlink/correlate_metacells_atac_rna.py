#!/usr/bin/env python
# coding: utf-8

import os
from scipy import sparse
import argparse
import logging
import pandas as pd
import numpy as np
import anndata

from robustlink import utils
from robustlink import enhancer_gene_utils

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
        output = enhancer_gene_utils.compute_enh_gene_corrs(
            gc_rna, ec_atac, 
            common_genes, np.arange(len(ec_atac)),
            enhancer_gene_to_eval['gene'].values, 
            enhancer_gene_to_eval['enh'].values, 
            output_file="", corr_type=corr_type, chunksize=100000, verbose_level=0,
            )

        # save results
        (to_correlate, corrs, corrs_shuffled, corrs_shuffled_cells) = output
        res_corrs = enhancer_gene_to_eval[to_correlate].copy()
        res_corrs['corr'] = corrs 
        res_corrs['corr_shuff'] = corrs_shuffled
        res_corrs['corr_shuff_cells'] = corrs_shuffled_cells
        res_corrs.to_csv(output_corr, sep='\t', header=True, index=True)

    return 

def wrap_corr_analysis_atac(
        out_dir,
        f_egpairs, f_gene, f_enh,
        scfusion_dir, mod_x, mod_y, input_name_tag, i_sub,
        corr_type='pearsonr',
        force=False,
        num_metacell_limit=0,
    ):
    """
    """
    # output: (i, k, --r)
    output_corrs = os.path.join(out_dir, f'{input_name_tag}_s{i_sub}_corr_{corr_type}_{{}}.tsv') 
            
    # for knn_xx
    input_modx_clsts =       [os.path.join(scfusion_dir, f'{input_name_tag}_s{i_sub}_metacells_{mod_x}.tsv.gz')]
    # for knn_xy
    input_knn_cells_xaxis  =  os.path.join(scfusion_dir, f'{input_name_tag}_s{i_sub}_cells_{mod_x}.txt')
    input_knn_cells_yaxis  =  os.path.join(scfusion_dir, f'{input_name_tag}_s{i_sub}_cells_{mod_y}.txt')
    input_knn_xy =            os.path.join(scfusion_dir, f'{input_name_tag}_s{i_sub}_knn_across_{mod_x}_{mod_y}.npz') 
    ### end of file paths

    # # Load data 
    # load enhancer-gene linkage
    enhancer_gene_to_eval = pd.read_csv(f_egpairs, sep='\t')

    # load count matrices
    adata_gene = anndata.read(f_gene)
    adata_enh  = anndata.read(f_enh)

    common_modx_cells = adata_gene.var.index.values
    common_genes = adata_gene.obs.index.values
    X = adata_gene.X

    common_mody_cells  = adata_enh.var.index.values 
    common_enhancer_regions = adata_enh.obs # check if this works
    Y = adata_enh.X

    # load knn networks 
    # for knn_xx 
    modx_clsts = pd.concat([
        pd.read_csv(fname, sep='\t',index_col=0)
        for fname in input_modx_clsts
    ], axis=1)
    # for knn_xy 
    knn_xy = sparse.load_npz(input_knn_xy)  
    cell_cell_knn_xaxis = np.loadtxt(input_knn_cells_xaxis, dtype=str)
    cell_cell_knn_yaxis = np.loadtxt(input_knn_cells_yaxis, dtype=str)

    logging.info("{} {} {} {}".format(modx_clsts.shape, knn_xy.shape, 
                                        cell_cell_knn_xaxis.shape, 
                                        cell_cell_knn_yaxis.shape,
                                        )
                )

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

def add_args(parser):
    """
    """
    parser.add_argument('--tolink',         required=True, type=str, help='a table of enhancer-gene pairs to be assessed for linkages')
    parser.add_argument('--countdata_gene', required=True, type=str, help='gene by cell count matrix; h5ad format')
    parser.add_argument('--countdata_enh',  required=True, type=str, help='enhancer by cell count matrix; h5ad format')
    parser.add_argument('-o', '--out_dir',  required=True, type=str)

    # to pull scFusion results (kNN graphs)
    parser.add_argument('--scfusion_dir',  required=True)
    parser.add_argument('--fusiondata_rna',  required=True, type=str, help='data for scFusion (h5ad file) -- need the filename to pull scFusion results')
    parser.add_argument('--fusiondata_mc',   required=True, type=str, help='data for scFusion (h5ad file) -- need the filename to pull scFusion results')
    parser.add_argument('-tag',  '--input_name_tag', required=True, type=str)
    parser.add_argument('-isub', '--i_sub', type=str, help="[0-9]")

    parser.add_argument('-ct', '--corr_type', choices=['pearsonr', 'spearmanr'],  
                                              default='pearsonr',
                                              help="choose from pearsonr or spearmanr")
    parser.add_argument('-f', '--force', action='store_true', help='to overwrite existing file')
    parser.add_argument('-n', '--num_metacell_limit', default=0, type=int, help='max num of metacells')
    return parser

def main(args):
    # output setting
    # run this with each combination of (i_sub, knn)
    f_egpairs = args.tolink
    f_gene = args.countdata_gene
    f_enh  = args.countdata_enh

    out_dir = args.out_dir
    scfusion_dir = args.scfusion_dir

    mod_x = os.path.basename(args.fusiondata_rna).replace('.h5ad', '')
    mod_y = os.path.basename(args.fusiondata_mc ).replace('.h5ad', '')

    input_name_tag = args.input_name_tag
    i_sub = args.i_sub
    corr_type = args.corr_type
    force = args.force
    num_metacell_limit = args.num_metacell_limit

    logging.info(",".join([mod_x, mod_y, input_name_tag, str(i_sub), corr_type, str(force), str(num_metacell_limit)]))

    wrap_corr_analysis_atac(
        out_dir,
        f_egpairs, f_gene, f_enh,
        scfusion_dir, mod_x, mod_y, input_name_tag, i_sub,
        corr_type=corr_type,
        force=force,
        num_metacell_limit=num_metacell_limit,
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    add_args(parser)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    main(args)