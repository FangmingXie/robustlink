#
from scipy import optimize
import matplotlib.pyplot as plt
import tqdm
import pickle
from scipy import stats
from scipy import sparse
import pandas as pd
import numpy as np
import time

from robustlink import utils

def set_venn_scale(ax, true_area, ref_area):
    s = np.sqrt(ref_area/true_area)
    ax.set_xlim(-s, s)
    ax.set_ylim(-s, s)
    return

def turn_cluster_labels_to_knn(cluster_labels, uniq_labels):
    """
    """

    clst_idx = utils.get_index_from_array(uniq_labels, cluster_labels)
    m, n = len(clst_idx), len(uniq_labels)
    _i = np.arange(m)
    _j = clst_idx
    _data = [1]*m
    knn = sparse.coo_matrix((_data, (_i, _j)), shape=(m, n))
    return knn

def row_dot_product_norm_by_numcol(X_zscore, Y_zscore, x_idx, y_idx,
                                   chunksize=100000, verbose_level=100000):
    """compute (X_zscore[x_idx]*Y_zscore[y_idx]).mean(axis=1)
    correlation values given matched x_idx and y_idx...
    """
    ti = time.time()

    assert len(x_idx) == len(y_idx)
    num_pairs = len(x_idx)
    corrs = []
    for pair_idx in utils.chunks(np.arange(num_pairs), chunksize):
        if verbose_level and pair_idx[0] % verbose_level == 0:
            print(pair_idx[0], time.time()-ti)

        _res = (X_zscore[x_idx[pair_idx]] *
                Y_zscore[y_idx[pair_idx]]).mean(axis=1)
        corrs.append(_res)
    corrs = np.hstack(corrs)
    return corrs

def compute_enh_gene_corrs(gc_rna, ec_mccg,
                           genes, enhancers,
                           pairs_gene, pairs_enh,
                           enhancer_groups=[],
                           output_file='', corr_type='pearsonr', shuff_enhs=False, **kwargs,
                           ):
    """Compute correlations and 2 shuffled correlations
    inputs:
        gc_rna: enh-by-cell RNA matrix
        ec_mccg: enh-by-cell mC matrix

        genes: gene list in the order of the mat gc_rna
        enhancers: enhancer list in the order of the mat ec_mccg
        pairs_gene: pairs to evalute (gene)
        pairs_enh: pairs to evaluate (enh)
        (referencing genes and enhancers,
        if a gene/enh is not in the genes/enhancers list, do not compute the correlations
        )
        enhancer_groups: an array with the same length as enhancers (and in the same order as enhancers)
    outputs:
        to_correlate # pairs computed correlations
        corrs, corrs_shuffled, corrs_shuffled_cells
    """
    print("{} chosen!".format(corr_type))
    if corr_type == 'pearsonr':
        pass
    elif corr_type == 'spearmanr':
        # ranking - across cells
        gc_rna = pd.DataFrame(np.array(gc_rna)).rank(axis=1).values
        ec_mccg = pd.DataFrame(np.array(ec_mccg)).rank(axis=1).values

    gc_rna_zscore = stats.zscore(
        np.array(gc_rna), axis=1, ddof=0, nan_policy='propagate')
    ec_mccg_zscore = stats.zscore(
        np.array(ec_mccg), axis=1, ddof=0, nan_policy='propagate')

    # correlate e-g according to a e-g table
    gene_idx = utils.get_index_from_array(genes, pairs_gene)
    enh_idx = utils.get_index_from_array(
        enhancers, pairs_enh)  # be careful here!
    to_correlate = ~np.logical_or(gene_idx == -1, enh_idx == -1)

    gene_idx = gene_idx[to_correlate]
    enh_idx = enh_idx[to_correlate]

    # corr
    corrs = row_dot_product_norm_by_numcol(
        gc_rna_zscore, ec_mccg_zscore, gene_idx, enh_idx, **kwargs)

    # corr shuffled cells
    corrs_shuffled_cells = row_dot_product_norm_by_numcol(
        gc_rna_zscore[:, np.random.permutation(gc_rna_zscore.shape[1])],
        ec_mccg_zscore,
        gene_idx, enh_idx, **kwargs)

    # corr shuffled genes (break up the pairs)
    gene_idx_uniq = np.unique(gene_idx)
    shuff_genes = {
        gene: gene_shuff for gene, gene_shuff in
        zip(gene_idx_uniq,
            gene_idx_uniq[np.random.permutation(len(gene_idx_uniq))])
    }
    gene_idx_shuff = np.array([shuff_genes[gene] for gene in gene_idx])

    corrs_shuffled = row_dot_product_norm_by_numcol(
        gc_rna_zscore,
        ec_mccg_zscore,
        gene_idx_shuff, enh_idx, **kwargs)

    if not shuff_enhs:
        output = (to_correlate, corrs, corrs_shuffled, corrs_shuffled_cells)
        # save and return
        if output_file:
            with open(output_file, 'wb') as fh:
                pickle.dump(output, fh)
        return output
    else:
        # corr shuffled enhancers
        # shuffle enhancers within each group
        if len(enhancer_groups) > 0:
            # control for the groupings
            enh_idx_uniq = np.unique(enh_idx) # unique enhancers
            enh_idx_groups = enhancer_groups[enh_idx] # enh_idx - enh_idx_groups
            enh_grps_uniq = np.unique(enh_idx_groups) # unique groups
            shuff_enhs = {}
            for enh_grp in enh_grps_uniq:
                enh_idx_grp = enh_idx[enh_idx_groups==enh_grp]
                enh_idx_grp_uniq = np.unique(enh_idx_grp)
                shuff_enhs.update({
                    enh: enh_shuff for enh, enh_shuff in
                    zip(enh_idx_grp_uniq,
                        enh_idx_grp_uniq[np.random.permutation(len(enh_idx_grp_uniq))])
                })
            enh_idx_shuff = np.array([shuff_enhs[enh] for enh in enh_idx])
            corrs_shuffled_enhs_bygroups = row_dot_product_norm_by_numcol(
                gc_rna_zscore,
                ec_mccg_zscore,
                gene_idx, enh_idx_shuff, **kwargs)
        else:
            corrs_shuffled_enhs_bygroups = 0

        # shuffle all enhancers
        enh_idx_uniq = np.unique(enh_idx)
        shuff_enhs = {
            enh: enh_shuff for enh, enh_shuff in
            zip(enh_idx_uniq,
                enh_idx_uniq[np.random.permutation(len(enh_idx_uniq))])
        }
        enh_idx_shuff = np.array([shuff_enhs[enh] for enh in enh_idx])
        corrs_shuffled_enhs = row_dot_product_norm_by_numcol(
            gc_rna_zscore,
            ec_mccg_zscore,
            gene_idx, enh_idx_shuff, **kwargs)

        output = (to_correlate, corrs, corrs_shuffled, corrs_shuffled_cells, 
                corrs_shuffled_enhs,
                corrs_shuffled_enhs_bygroups,
                )
        # save and return
        if output_file:
            with open(output_file, 'wb') as fh:
                pickle.dump(output, fh)
        return output

DEBUG_MODE = False
def get_r_threshold_smart(bins, fdr, fdr_threshold, side='left'):
    """given fdr function (bins, fdr) and fdr_threshold
    get r threshold for a given fdr_treshold
    """
    # get r_threshold
    # remove nan
    flag = 1  # default no solution

    isnan = np.isnan(fdr)
    _y = fdr[~isnan]
    _x = bins[1:][~isnan]

    # find r threshold
    def f(_x_func): return np.interp(_x_func, _x, _y) - fdr_threshold

    # solve
    if side == 'left':
        bins_sub = bins[bins < 0]
    elif side == 'right':
        bins_sub = bins[bins > 0]
    else:
        raise ValueError("left or right")

    r_min = bins_sub[np.argmin(f(bins_sub))]
    if f(r_min)*f(0) <0:
        sol = optimize.root_scalar(f, bracket=(r_min, 0))
        if sol:
            r_threshold = sol.root
            flag = 0

    if flag == 1:
        r_threshold = np.nan
        print("failed to detect r_threshold:")

    # if DEBUG_MODE and flag == 1:
    if DEBUG_MODE:
        fig, ax = plt.subplots()
        ax.plot(bins, f(bins))
        if flag == 0:
            print(r_threshold, f(r_threshold))
            ax.scatter([r_threshold], [f(r_threshold)])
        # ax2 = ax.twinx()
        # ax2.plot(bins[1:], f(bins[1:]) - f(bins[:-1]), color='C1')
        plt.show()

    return r_threshold


def cumfrac_to_pval(cumfracs, pval_type):
    """Null hypothese
    """
    if pval_type == 'left':
        pvals = cumfracs
    elif pval_type == 'right':
        pvals = 1 - cumfracs
    elif pval_type == 'both':
        alpha = np.minimum(cumfracs, 1-cumfracs)
        pvals = 2*alpha
    else:
        raise ValueError
    return pvals


def cumfrac_to_pval_obs_simple(cumfracs, pval_type):
    """Observations
    Assuming null is symmetric given the bin range
    """
    if pval_type == 'left':
        pvals = cumfracs
    elif pval_type == 'right':
        pvals = 1 - cumfracs
    elif pval_type == 'both':
        alpha = np.minimum(cumfracs, 1-cumfracs)
        alpha_star = alpha[::-1]  # (this implies bins_conj == bins[::-1])
        pvals = alpha + alpha_star
    else:
        raise ValueError
    return pvals


def cumfrac_to_pval_obs_general(bins, cumfracs, cumfracs_null, pval_type):
    """Observations
    """
    if pval_type == 'left':
        pvals = cumfracs
    elif pval_type == 'right':
        pvals = 1 - cumfracs
    elif pval_type == 'both':
        alpha_null = np.minimum(cumfracs_null, 1-cumfracs_null)
        alpha = np.minimum(cumfracs, 1-cumfracs)
        # look for bins_conj such that alpha_null(bins_conj) == alpha_null(bins)
        # bins_conj = get_conjugate(bins, alpha_null)
        # alpha_star = # alpha(bins_conj)
        # pvals = alpha + alpha_star
        raise ValueError("Not yet implemented")
    else:
        raise ValueError
    return pvals


def get_significance_stats(
            pairs,
            corrs, corrs_shuffled, corrs_shuffled_cells,
            pval_type_shuffled, pval_type_shuffled_cells,
            bins=np.linspace(-1, 1,101),
            distance_threshold=1e5, fdr_threshold=0.2,
            positive_side=False,
            return_pval=False,
            return_cdf=False,
        ):
    """assuming significant negative correlation
    pairs should be a df with at least 3 columns: gene, enh, dist
    """
    # align all to negative side
    assert pval_type_shuffled in ['left', 'both', 'right']  # one/two side test
    assert pval_type_shuffled_cells == 'both'  # two sides test

    if positive_side:
        # prevent modifying the original one
        corrs = corrs.copy()
        corrs_shuffled = corrs_shuffled.copy()
        corrs_shuffled_cells = corrs_shuffled_cells.copy()

        corrs *= (-1)
        corrs_shuffled *= (-1)
        corrs_shuffled_cells *= (-1)

    # dists
    dists = pairs['dist'].values
    label_cond = dists < distance_threshold
    track = corrs[label_cond]

    # total numbers with the condition
    num_total_pairs = len(pairs[label_cond])
    num_total_genes = len(pairs[label_cond]['gene'].unique())
    num_total_enhs = len(pairs[label_cond]['enh'].unique())

    # hist_shuff
    hist_shuff, _ = np.histogram(corrs_shuffled, bins=bins, density=True)
    cdf_shuff = np.cumsum(hist_shuff)/np.sum(hist_shuff)
    pval_shuff = cumfrac_to_pval(cdf_shuff, pval_type_shuffled)
    # hist_shuff_cells
    hist_shuff_cells, _ = np.histogram(corrs_shuffled_cells, bins=bins, density=True)
    cdf_shuff_cells = np.cumsum(hist_shuff_cells)/np.sum(hist_shuff_cells)
    pval_shuff_cells = cumfrac_to_pval(cdf_shuff_cells, pval_type_shuffled_cells)
    # hist
    hist, _ = np.histogram(track, bins=bins, density=True)
    cdf = np.cumsum(hist)/np.sum(hist)
    # modified on 12/23/2020
    pval_obs_shuff = cumfrac_to_pval_obs_simple(cdf, pval_type_shuffled)
    pval_obs_shuff_cells = cumfrac_to_pval_obs_simple(cdf, pval_type_shuffled_cells)

    # fdr
    fdr_linked = (pval_shuff+1e-7)/(pval_obs_shuff+1e-7)
    fdr_correlated = (pval_shuff_cells+1e-7)/(pval_obs_shuff_cells+1e-7)

    # get all fdrs
    track_fdr_linked = np.interp(track, bins[1:], fdr_linked)
    track_fdr_correlated = np.interp(track, bins[1:], fdr_correlated)

    # get r_threshold
    if pval_type_shuffled == 'both':
        r_threshold_linked = np.nan
        r_threshold_linked_left = get_r_threshold_smart(bins, fdr_linked, fdr_threshold, side='left')
        r_threshold_linked_right = get_r_threshold_smart(bins, fdr_linked, fdr_threshold, side='right')
    else:
        r_threshold_linked = get_r_threshold_smart(bins, fdr_linked, fdr_threshold, side=pval_type_shuffled)
        r_threshold_linked_left = np.nan
        r_threshold_linked_right = np.nan

    r_threshold_correlated_left = get_r_threshold_smart(bins, fdr_correlated, fdr_threshold, side='left')
    r_threshold_correlated_right = get_r_threshold_smart(bins, fdr_correlated, fdr_threshold, side='right')

    # stats
    # all sig info
    linked_table = pairs[label_cond][track_fdr_linked < fdr_threshold]
    num_linked_pairs = len(linked_table)
    num_linked_genes = len(linked_table['gene'].unique())
    num_linked_enhs = len(linked_table['enh'].unique())
    # all sig info
    correlated_table = pairs[label_cond][track_fdr_correlated < fdr_threshold]
    num_correlated_pairs = len(correlated_table)
    num_correlated_genes = len(correlated_table['gene'].unique())
    num_correlated_enhs = len(correlated_table['enh'].unique())

    output = {
        'dist_th': distance_threshold,
                'num_total_pairs': num_total_pairs,
                'num_total_genes': num_total_genes,
                'num_total_enhs': num_total_enhs,

                'r_th_linked': r_threshold_linked,
                'r_th_linked_left': r_threshold_linked_left,
                'r_th_linked_right': r_threshold_linked_right,

                'r_th_correlated_left': r_threshold_correlated_left,
                'r_th_correlated_right': r_threshold_correlated_right,

                'num_linked_pairs': num_linked_pairs,
                'num_linked_genes': num_linked_genes,
                'num_linked_enhs': num_linked_enhs,
                'linked_table': linked_table,

                'num_correlated_pairs': num_correlated_pairs,
                'num_correlated_genes': num_correlated_genes,
                'num_correlated_enhs': num_correlated_enhs,
                'correlated_table': correlated_table,
        }

    if return_pval:
        output['bins'] = bins
        output['linked_pval'] = pval_shuff
        output['correlated_pval'] = pval_shuff_cells
    if return_cdf:
        output['bins'] = bins
        output['linked_cdf'] = cdf_shuff
        output['correlated_cdf'] = cdf_shuff_cells
    return output

def get_corr_stats(iterator_both, enhancer_gene_to_eval, col_orders,
    bins=np.linspace(-1, 1,101),
    distance_threshold=1e5,
                   fdr_threshold=0.2,
                   r_min=-1,
                   r_max=0,
    ):
    """Wrap-up the routine of calculating corr stats, and compare between mC, ATAC, and both
    used in jupyter notebook visualizations
    """
    pval_type_shuffled, pval_type_shuffled_cells = 'left', 'both'
    positive_sides = [False, True]
    res = []
    for idx, row in tqdm.tqdm(iterator_both.iterrows()):
        fname_mc, fname_atac = (row['fname_mc'], row['fname_atac'],)

        res_2cases = []
        for i, (label, fname, positive_side) in enumerate(zip(
            ['mc', 'atac'],
                            [fname_mc, fname_atac],
                            [False, True],
            )):
            try:
                with open(fname, 'rb') as fh:
                    to_correlate, corrs, corrs_shuffled, corrs_shuffled_cells = pickle.load(fh)
            except:
                continue

            pairs = enhancer_gene_to_eval[to_correlate].copy()
            label_cond = pairs['dist'].values < distance_threshold

            res_1case = get_significance_stats(
                pairs,
                corrs, corrs_shuffled, corrs_shuffled_cells,
                pval_type_shuffled, pval_type_shuffled_cells,
                bins=bins,
                distance_threshold=distance_threshold, fdr_threshold=fdr_threshold,
                positive_side=positive_side,
                return_cdf=False,
                return_pval=False,
            )

            # record
            res_1case['id_total_pairs'] = pairs[label_cond].index.values

            if isinstance(res_1case['linked_table'], pd.DataFrame):
                res_1case['id_linked_pairs'] = res_1case['linked_table'].index.values
            else:
                res_1case['id_linked_pairs'] = np.array([])

            if isinstance(res_1case['correlated_table'], pd.DataFrame):
                res_1case['id_correlated_pairs'] = res_1case['correlated_table'].index.values
            else:
                res_1case['id_correlated_pairs'] = np.array([])

            res_1case = pd.Series(res_1case)[col_orders]
            res_1case.index = res_1case.index + '_{}'.format(label)
            res_2cases.append(res_1case)

        # coommon
        res_2cases = pd.concat(res_2cases)

        # linked both
        common_pair_ids = np.intersect1d(res_2cases['id_linked_pairs_mc'],
                                         res_2cases['id_linked_pairs_atac'],
                                         )
        common_gene_ids = enhancer_gene_to_eval.loc[common_pair_ids, 'gene'].unique()
        common_enh_ids = enhancer_gene_to_eval.loc[common_pair_ids, 'enh'].unique()
        res_2cases['num_linked_pairs_both'] = len(common_pair_ids)
        res_2cases['num_linked_genes_both'] = len(common_gene_ids)
        res_2cases['num_linked_enhs_both'] = len(common_enh_ids)

        # correlated both
        common_pair_ids = np.intersect1d(res_2cases['id_correlated_pairs_mc'],
                                         res_2cases['id_correlated_pairs_atac'],
                                         )
        common_gene_ids = enhancer_gene_to_eval.loc[common_pair_ids, 'gene'].unique()
        common_enh_ids = enhancer_gene_to_eval.loc[common_pair_ids, 'enh'].unique()
        res_2cases['num_correlated_pairs_both'] = len(common_pair_ids)
        res_2cases['num_correlated_genes_both'] = len(common_gene_ids)
        res_2cases['num_correlated_enhs_both'] = len(common_enh_ids)

        # total both
        common_pair_ids = np.intersect1d(res_2cases['id_total_pairs_mc'],
                                         res_2cases['id_total_pairs_atac'],
                                         )
        common_gene_ids = enhancer_gene_to_eval.loc[common_pair_ids, 'gene'].unique()
        common_enh_ids = enhancer_gene_to_eval.loc[common_pair_ids, 'enh'].unique()
        res_2cases['num_total_pairs_both'] = len(common_pair_ids)
        res_2cases['num_total_genes_both'] = len(common_gene_ids)
        res_2cases['num_total_enhs_both'] = len(common_enh_ids)

        res.append(pd.Series(res_2cases))
    res = pd.DataFrame(res, index=iterator_both.index)
    print(res.shape)
    return res
