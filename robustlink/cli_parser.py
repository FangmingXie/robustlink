"""Command line interface is defined here.
"""
DESCRIPTION="""
robustlink is a computational tool that links candidate enhancers to genes using single-cell transcriptome and epigenome datasets. 
"""

EPILOG="""
Contributors: Fangming Xie, Ethan J. Armand, Wayne I. Doyle, Eran Mukamel.
Contact: Eran Mukamel (emukamel@ucsd.edu).
"""

import argparse
import os

def create_parser():
    """
    """
    parser = argparse.ArgumentParser(
        prog="robustlink",
        description=DESCRIPTION,
        epilog=EPILOG,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    required = parser.add_argument_group('required')
    optional = parser.add_argument_group('optional') 
    advanced = parser.add_argument_group('advanced')
    
    ## ARGUMENTS DIRECTLY FED INTO SingleCellFusion CLI 
    # Input/Output Dataset Settings
    required.add_argument(
        "-i", "--input_datasets", 
        metavar="xx.h5ad",
        type=str,
        nargs="+",
        required=True,
        help='''(list of str) 
             Paths to .h5ad files, each containing a cell-by-gene feature matrix, 
             cell IDs and gene IDs. Cell IDs should be unique within each .h5ad file, 
             Gene IDs should be shared or partially shared across files. 
             Multiple inputs should be listed as a space seperated list of filenames.
             '''
    )
    required.add_argument(
        "-im", "--input_modalities", 
        metavar="rna/atac/mc",
        type=str,
        nargs="+",
        required=True,
        help='''(list of str)
             Data modalities chosen from 'rna', 'atac', or 'mc'. This should be 
             listed in the same order as input_datasets. 
             '''
    )
    # may need this in the future
    # parser.add_argument(
    #     "-im", "--input_meta", 
    #     type=str,
    #     required=True,
    #     help="(list of str) Input metadata csv file",
    # )

    required.add_argument(
        "-f", "--feature_datasets", 
        metavar="xx.h5ad",
        type=str,
        nargs="+",
        required=True,
        help='''(list of str)
             Dataset(s) whose features all other datasets will impute into.
             This should be a subset of --input_datasets.
             Enter multiple datasets as a space-separated list of filenames.
             The features of these datasets will
             be the features kept in the output imputed data table.",
             '''
    )
    optional.add_argument(
        "-o", "--output_dir", 
        metavar="DIR",
        type=str,
        default="./results",
        help='''(str)
             Directory to store output files
             '''
    )
    optional.add_argument(
        "-op", "--output_prefix", 
        type=str,
        default="SingleCellFusion",
        help='''(str)
             The output files will contain this prefix.
             '''
    )

    # constraint kNN across modalities
    optional.add_argument(
        "--nearest_neighbors", 
        type=int,
        default=20,
        help='''(integer)
             Number of nearest neighbors used to impute data
             '''
    )
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
    optional.add_argument(
        "--precomputed_pca_file", 
        type=str,
        default='',
        help='''(str)
             Precomputed PCA matrix (tab separated table; text file or gzipped)
             with the first row as the header, and the first column as the cell_id
             Each following rows are cell features, and columns are PCs.

             Providing this file will by-pass SingleCellFusion integration, 
             and do clustering and UMAP just on this matrix instead.
             '''
    )
    optional.add_argument(
        "--use_netUMAP", 
        action='store_true',
        help='''(bool)
             Include this argument to use Net-UMAP from Pegasus (Li et al. 2020)
             Net-UMAP is an approximate but fast algorithm for UMAP.
             It runs traditional UMAP on a subset of cells, 
             then it uses deep neural network to learn embedding for all cells.
             The package pegasus is required.
             '''
    )
    optional.add_argument(
        "--use_tsne", 
        action='store_true',
        help='''(bool)
             Include this argument to use tSNE instead of UMAP
             '''
    )

    # within modality smoothing 
    advanced.add_argument(
        "--num_pcs", 
        type=int,
        default=50,
        help='''(integer)
             Number of Principal Components to keep for each dataset 
             for smoothing and for clustering/embedding after imputation.
             '''
    )
    advanced.add_argument(
        "--smoothing_fractions", 
        nargs="+",
        type=float,
        default=[0.7, 0.1, 0.9],
        help='''(list of floats) 
             A list of three values between 0 to 1 that controls the relative contribution
             from the cell itself vs. its neighbors in within-dataset smoothing, 
             specified for 'rna', 'atac', 'mc' data, respectively.
             '''
    )

    # Arguments for Clustering
    advanced.add_argument(
        "--leiden_n_neighbors", 
        type=int,
        default=30,
        help='''(integer) 
             Number of nearest neighbors to form in the integrated space, 
             the resulting nearest neighbor graph is used for Leiden clustering.
             It is passed into the python package leidenalg.
             '''
    )
    advanced.add_argument(
        "--leiden_resolutions", 
        type=list,
        default=[0.1, 0.2, 0.4, 0.8],
        help='''(list of floats) 
             A list of resolutions to be used for Leiden Clustering.
             It is passed into the python package leidenalg.
             '''
    )

    # Arguments for UMAP
    advanced.add_argument(
        "--umap_n_neighbors", 
        type=int,
		default=60,
        help='''(integer)
             Number of neighbors for UMAP. It is passed into the python package umap.UMAP(n_neighbors).
             '''
    )
    advanced.add_argument(
        "--umap_min_dist", 
        type=float,
        default=0.5,
        help='''(float)
             Minimum distance for UMAP. It is passed into the python package umap.UMAP(min_dist).
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
