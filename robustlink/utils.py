"""
"""
import logging
import numpy as np
import pandas as pd
import os

def logcpm(counts):
    """
    Args:
        - gene-cell matrix (pandas DataFrame)
    """
    cov = counts.sum(axis=0)
    logcpm = np.log10(counts.divide(cov, axis=1)*1000000 + 1)
    return logcpm

def logtpm(counts, gene_lengths):
    """
    Args:
        - gene-cell matrix (pandas DataFrame)
        - gene_lengths: a series indexed by gene_id
    """
    tpm = counts.divide(gene_lengths.loc[counts.index], axis=0)
    cov = tpm.sum(axis=0)
    logtpm = np.log10((tpm.divide(cov, axis=1))*1000000 + 1)
    return logtpm

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def create_logger(name='log'):
    """
    args: logger name
    return: a logger object
    """
    logging.basicConfig(
        format='%(asctime)s %(message)s', 
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.INFO)
    return logging.getLogger(name)

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def get_mcc_lite_v4(df_c, df_mc, base_call_cutoff, sufficient_coverage_fraction=1, fillna=True):
    """
    """
    # a gene is sufficiently covered in % of cells 
    condition = (df_c > base_call_cutoff).sum(axis=1) >= sufficient_coverage_fraction*(df_c.shape[1]) 

    # get mcc matrix with kept bins and nan values for low coverage sites
    df_c_nan = df_c.copy()
    df_c_nan[df_c < base_call_cutoff] = np.nan
    df_mcc = df_mc.loc[condition]/df_c_nan.loc[condition]

    # imputation (missing value -> mean value of all cells)
    if fillna:
        logging.info('Imputing data... (No effect if sufficient_coverage_fraction=1)')
        means = df_mcc.mean(axis=1)
        fill_value = pd.DataFrame({col: means for col in df_mcc.columns})
        df_mcc.fillna(fill_value, inplace=True)
    
    return df_mcc

def get_index_from_array(arr, inqs, na_rep=-1):
    """Get index of array
    """
    arr = np.array(arr)
    arr = pd.Series(arr).reset_index().set_index(0)
    idxs = arr.reindex(inqs)['index'].fillna(na_rep).astype(int).values
    return idxs

def import_single_textcol(fname, header=None, col=0):
    return pd.read_csv(fname, header=header, sep='\t')[col].values

def export_single_textcol(fname, array):
    with open(fname, 'w') as f:
        f.write('\n'.join(array)+'\n')

def combine_legends(axs_flat):
    """Combine legends from subplots (axs.flat)
    """
    handles_all = [] 
    labels_all = []
    for ax in axs_flat:
        handles, labels = ax.get_legend_handles_labels()
        handles_all.append(handles)
        labels_all.append(labels)
    
    return np.hstack(handles_all), np.hstack(labels_all)
