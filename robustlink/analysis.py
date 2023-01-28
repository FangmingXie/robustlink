"""For analysis of correlation results
"""
from robustlink._init_plots import *
from robustlink import enhancer_gene_utils

KB = 1000
def p25(x):
    return np.nanpercentile(x, 25)

def p75(x):
    return np.nanpercentile(x, 75)

def estimate_frac_tps(pvalues, bin_n=100, frac_bin_null=0.20):
    """Use the median of the last 5% (frac_bin_null) of bins to estimate null level
    """
    y = np.sort(pvalues)
    bin_edges = np.linspace(0, 1, bin_n)
    bin_width = 1.0/bin_n
    
    bin_counts, bin_edges = np.histogram(y, bin_edges)
    y_norm = bin_counts/(bin_width*bin_counts.sum())
    null_level = np.median(y_norm[-int(frac_bin_null*bin_n):])
    frac_tp = np.clip(1 - null_level, 0, 1)
    
    return frac_tp, null_level, bin_counts, bin_edges, y_norm
    
def plot_pval_dist(ax, frac_tp, null_level, bin_counts, bin_edges, y_norm, fillcolor='C0'):
    """
    """
    _x = bin_edges[:-1]
    ax.plot(_x, y_norm, color='k')
    ax.plot([0, 1], [null_level]*2, linestyle='--', color='k')
    ax.fill_between(_x, null_level, y_norm, 
                    where=y_norm>null_level, alpha=1, color=fillcolor)
    ax.fill_between(_x, 0, np.minimum(y_norm, null_level), alpha=0.5, color='lightgray')
    ax.text(0, null_level+0.2*(1-null_level), 
            "{:.1f}%".format(frac_tp*100), 
            color='white',
            fontsize=18)
    sns.despine(ax=ax)
    return ax

def pipe_plot_pval_dist(ax, pvalues, bin_n=100, frac_bin_null=0.20, fillcolor='C0'):
    """
    """
    # fit distribution 
    frac_tp, null_level, bin_counts, bin_edges, y_norm = estimate_frac_tps(pvalues, bin_n=bin_n, frac_bin_null=frac_bin_null)
    # plot it out 
    plot_pval_dist(ax, 
                   frac_tp, null_level, bin_counts, bin_edges, y_norm, 
                   fillcolor=fillcolor,
                  )
    return 

def quantile_norm(array):
    return pd.Series(array).rank(pct=True, method='average').values

def estimate_frac_tps_vs_dists(dists_kb, res_corrs, res_stats, link_type, FLIP_CORR_SIGN, 
                               bin_n=51, frac_bin_null=0.1,
                              ):
    """Estimate the fraction of true positives with many different distance segments
    """
    frac_tps = []
    num_tps = []
    total_nums = []
    for idx in np.arange(len(dists_kb)):
        dist_kb = dists_kb[idx]
        if idx == 0:
            cond = (res_corrs['dist'] < dist_kb*KB)
        else:
            dist_kb_prev = dists_kb[idx-1]
            cond = ((res_corrs['dist'] < dist_kb*KB) & 
                    (res_corrs['dist'] >= dist_kb_prev*KB))

        if FLIP_CORR_SIGN: 
            corr_sign = -1
        else:
            corr_sign = 1

        pval = np.interp(corr_sign*res_corrs[cond]['corr'], 
                                res_stats['bins'][1:], res_stats[link_type])
        frac_tp, null_level, bin_counts, bin_edges, y_norm = estimate_frac_tps(pval, bin_n=bin_n, frac_bin_null=frac_bin_null)
        frac_tps.append(frac_tp)
        total_nums.append(len(pval))
        num_tps.append(len(pval)*frac_tp)
    return dists_kb, frac_tps, num_tps, total_nums

def estimate_num_sigs_vs_dists(fdr_th, dists_kb, res_corrs, res_stats, link_type, FLIP_CORR_SIGN):
    """Estimate the number of significant pairs with many different distance segments
    """
    frac_tps = []
    num_tps = []
    total_nums = []
    
    sig_gene_set = set()
    cum_sig_genes = []
    
    for idx in np.arange(len(dists_kb)):
        dist_kb = dists_kb[idx]
        if idx == 0:
            cond = (res_corrs['dist'] < dist_kb*KB)
        else:
            dist_kb_prev = dists_kb[idx-1]
            cond = ((res_corrs['dist'] < dist_kb*KB) & 
                    (res_corrs['dist'] >= dist_kb_prev*KB))

        if FLIP_CORR_SIGN: 
            corr_sign = -1
        else:
            corr_sign = 1

        pval = np.interp(corr_sign*res_corrs[cond]['corr'], 
                                res_stats['bins'][1:], res_stats[link_type])
        fdr = pval/quantile_norm(pval)
        
        total_num = len(pval)
        total_nums.append(total_num)
        
        sig_cond = fdr<fdr_th
        sig_genes = res_corrs[cond][sig_cond]['gene']
        sig_gene_set = sig_gene_set.union(set(list(sig_genes.values)))
        cum_sig_gene = len(sig_gene_set)
        cum_sig_genes.append(cum_sig_gene)
        
        num_tp = sig_cond.sum()
        num_tps.append(num_tp)
        
        frac_tp = num_tp/total_num
        frac_tps.append(frac_tp)
        
    return dists_kb, frac_tps, num_tps, total_nums, cum_sig_genes

class CorrRes():
    def __init__(self, res_corrs, pcorr, label="enh-gene", color="C0", lightcolor="lightblue"):
        """
        """
        self.res_corrs = res_corrs
        assert isinstance(pcorr, bool)
        self.pcorr = pcorr # True or False
        self.label = label
        self.color = color
        self.lightcolor = lightcolor
        self.res_stats = None
        
        if pcorr: # positive correlation has a sign of -1 by convention
            self.sign = -1 
        else:
            self.sign = 1
        
    def test_significance(self, fdr=0.2, dist_th=1e5, nbins=501, 
                          output_linked=None, output_correlated=None): # 100kb
        pcorr = self.pcorr
        res_corrs = self.res_corrs
        pval_type_shuffled, pval_type_shuffled_cells = 'left', 'both'
        res_stats = enhancer_gene_utils.get_significance_stats(
                                    res_corrs[['gene', 'enh', 'dist']],
                                    res_corrs['corr'],
                                    res_corrs['corr_shuff'],
                                    res_corrs['corr_shuff_cells'],
                                    pval_type_shuffled, pval_type_shuffled_cells,
                                    bins=np.linspace(-1,1,nbins),
                                    distance_threshold=dist_th,
                                    fdr_threshold=fdr,
                                    positive_side=pcorr,
                                    return_pval=True,
                                    return_cdf=False,
                                )
        if output_linked is not None:
            res_stats['linked_table'].to_csv(output_linked, sep="\t", header=True, index=False)
            print(f"saved to {output_linked}")
        if output_correlated is not None:
            res_stats['correlated_table'].to_csv(output_correlated, sep="\t", header=True, index=False)
            print(f"saved to {output_correlated}")
            
        self.res_stats = res_stats
    
    def estimate_dist_dependence(self, fdr=0.2):
        """Estimate the fraction of true links among all pairs (as a function of distance)
        """
        dists_kb = np.arange(20, 500+1, 20)
        dists_kb_plot = np.hstack([[2], dists_kb])
        
        res_corrs = self.res_corrs
        res_stats = self.res_stats
        if self.pcorr: # atac - positive correlation -- need flipping by convention
            _flip = True
        else:
            _flip = False
        
        for link_type in ['linked_pval', 'correlated_pval']:
            _, frac_tps, num_tps, total_nums = estimate_frac_tps_vs_dists(
                    dists_kb, res_corrs, res_stats, link_type, FLIP_CORR_SIGN=_flip, 
                    bin_n=51, frac_bin_null=0.1,)
            _, sig_frac_tps, sig_num_tps, sig_total_nums, cum_sig_genes = estimate_num_sigs_vs_dists(
                    fdr, 
                    dists_kb, res_corrs, res_stats, link_type, FLIP_CORR_SIGN=_flip)

            # organize results
            mats = np.vstack([
                        dists_kb, 
                        frac_tps, num_tps, total_nums,
                        sig_frac_tps, sig_num_tps, sig_total_nums, cum_sig_genes,
                        ]).T 
            cols = ['dist', 
                    'frac_tp', 'num_pos', 'num_total', 
                    'sig_frac_tp', 'sig_num_pos', 'sig_num_total', 'cum_sig_genes',
                   ] 
            res = pd.DataFrame(mats, columns=cols) 
            if link_type == 'linked_pval': self.distdep_linked = res
            if link_type == 'correlated_pval': self.distdep_correlated = res
            
        return 
        
    def plot_corr_vs_dist(self, ax=None):
        """
        """
        res_corrs = self.res_corrs
        label = self.label
        color = self.color
        if ax is None:
            fig, ax = plt.subplots(figsize=(6,4))
        
        n = 100
        _dists = np.linspace(2*KB, 1000*KB, n)
        res_corrs_median = (res_corrs[['corr']].groupby(pd.cut(res_corrs['dist'], _dists)) 
                                               .agg([np.median, np.mean, p25, p75,])['corr']
                            )

        _x = _dists[1:]
        ax.hlines(0, 0, np.max(_x), linestyle='--', color='k')

        _yme = res_corrs_median['median'].values
        _ylo = res_corrs_median['p25'].values
        _yhi = res_corrs_median['p75'].values
        ax.plot(_x, _yme, linewidth=3, label=label, color=color)
        ax.fill_between(_x, _ylo, _yhi, alpha=0.2, color=color)
        sns.despine(ax=ax)

        ax.set_title('All enhancer-gene pairs\n2kb - 1Mb')
        ax.set_ylabel('Spearman correlation\n(median +/- interquartile)')
        ax.set_xlabel('Enhancer - gene (TSS) distance')
        ax.legend(bbox_to_anchor=(1,0), loc='lower right')
        ax.xaxis.set_major_formatter(mtick.EngFormatter())
        ax.set_ylim([-0.5, 0.5])
        ax.grid(False)
    
    def plot_pval(self, ax=None, dist_th=5e5, bin_n=51, frac_bin_null=0.05, pval_type='linked'):
        """
        """
        label = self.label
        color = self.color
        sign  = self.sign
        res_corrs = self.res_corrs
        res_stats = self.res_stats
        
        if ax is None:
            fig, ax = plt.subplots(figsize=(8,6)) 
        
        if pval_type == 'linked':
            xlabel = 'p (shuff regions as null)'
            stat_col = 'linked_pval'
        elif pval_type == 'correlated':
            xlabel = 'p (shuff metacells as null)'
            stat_col = 'correlated_pval'
        else:
            raise ValueError("choose from linked or correlated")
            
        pvals = np.interp(sign*res_corrs[res_corrs['dist']<dist_th]['corr'], 
                          res_stats['bins'][1:], res_stats[stat_col])
        pipe_plot_pval_dist(ax, pvals, bin_n=bin_n, frac_bin_null=frac_bin_null, fillcolor=color)
        
        ax.set_ylim([0, 2])
        ax.annotate(label, xy=(1,1), xycoords='axes fraction', ha='right', va='top')
        ax.set_xlabel(xlabel)
        ax.set_ylabel('Density')
    
    def plot_dist_dep(self, col, link_type, ax=None):
        """
        col: choose from 
            - frac_tp:     fraction of true positives
            - sig_frac_tp: significant fraction of positives
            - num_pos:     number of positives
            - sig_num_pos: significant number of positives
            - cum_sig_genes: number of significant genes
            
        link_type: choose from
            - linked
            - correlated
        """
        col_labels = {
            'frac_tp': 'Estimated fraction of true links', 
            'sig_frac_tp': 'Fraction of individually\n significant pairs', 
            'num_pos': 'Estimated cumulative\n number of true links',
            'sig_num_pos': 'Cumulative number of\n significant pairs',
            'cum_sig_genes': 'Cumulative number of\n significant genes',
            }
        
        if ax is None:
            fig, ax = plt.subplots(figsize=(5, 4))
            
        df = getattr(self, f"distdep_{link_type}")
        if   link_type == 'linked':     marker = 'o'
        elif link_type == 'correlated': marker = '^'
        
        _x = df['dist'] 
        _y = df[col]
        if col in ['num_pos', 'sig_num_pos']:
            _y = np.cumsum(_y)
        if col in ['num_pos', 'sig_num_pos', 'cum_sig_genes']:
            ax.set_yscale('log')
            ax.yaxis.set_major_formatter(mtick.EngFormatter()) 
        if col in ['frac_tp', 'sig_frac_tp']:
            ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, _: '{:.0%}'.format(x))) 
            
        ylabel = col_labels[col]
        ax.plot(_x, _y, f'-{marker}', color=self.color, label=f"{self.label} {link_type}", markersize=7)
        ax.set_xlabel('Enhancer - TSS distance (kb)')
        ax.set_ylabel(ylabel)
        sns.despine(ax=ax)
    
    def plot_corr_distribution(self, ax=None, 
                               bins=np.linspace(-1,1,101),
                              ):
        """
        """
        res_corrs = self.res_corrs
        res_stats = self.res_stats
        sign = self.sign
        if ax is None:
            fig, ax = plt.subplots(figsize=(7,4))
        fontsize = ax.xaxis.label.get_fontsize()
                             
        with sns.axes_style('ticks', {'axes.grid': False}):
            ax.set_title(self.label, pad=20)
            labels = ['shuffled metacells',
                      'shuffled genes',
                      'pairs <500kb',
                      'pairs <100kb',
                      ]
            colors = ['gray', 'black', self.lightcolor, self.color]
            corr_tracks = [
                   res_corrs['corr_shuff_cells'].values,
                   res_corrs['corr_shuff'].values,
                   res_corrs.loc[res_corrs['dist']<=500*KB, 'corr'].values,
                   res_corrs.loc[res_corrs['dist']<=100*KB, 'corr'].values,
                  ]
            vertical_lines = [
                sign*res_stats['r_th_linked'],
                sign*res_stats['r_th_correlated_left'],
                sign*res_stats['r_th_correlated_right'],
            ]

            histy_max = 0
            for j, (_x, label, color) in enumerate(zip(
                    corr_tracks, labels, colors, )):
                # go over columns
                label_comp = '{} ({})'.format(label, len(_x))
                g = ax.hist(_x, bins=bins, histtype='step', label=label, color=color, density=True)
                histy, histx, _ = g
                if j == 0:
                    histy0_max = np.max(histy)
                    histx0_max = bins[np.argmax(histy)]
                elif j > 0:
                    histy_max = max(histy_max, np.max(histy)) 

            ax.set_ylim([0, 1.3*histy_max])
            if histy0_max > 1.3*histy_max:
                # text
                text_config = {
                    'xy': (histx0_max, 1.3*histy_max), 
                    'ha': 'center', 'va': 'bottom', 
                    'xytext': (0, 0),
                    'textcoords': 'offset points',
                    'fontsize': 0.7*fontsize,
                }
                ax.annotate("{:.2f}".format(histy0_max), **text_config)

            # labels
            ax.set_xlim([-1, 1])
            ax.set_xlabel('Spearman correlation')
            ax.set_ylabel('Density')
            sns.despine(ax=ax)

            # line ticks
            ax.grid(which='major', axis='x', linestyle='--')
            ax.xaxis.set_major_formatter(mtick.StrMethodFormatter('{x:.2f}'))
            ax.set_xticks(np.sort(np.hstack([[-1, 0, 1], vertical_lines])))

            # horizontal lines
            lineys = [1.1*histy_max, 1.2*histy_max, 1.2*histy_max]
            linecolors = ['k', 'gray', 'gray']
            texts = ['linked', 'correlated', 'correlated']
            vas = ['top', 'bottom', 'bottom']
            offsets = [(0.2*fontsize, -0.2*fontsize), 
                       (0.2*fontsize, +0.2*fontsize), 
                       (0.2*fontsize, +0.2*fontsize), 
                      ]
            for xcoord, linecolor, liney, text, va, offset in zip(
                vertical_lines, linecolors, lineys, texts, vas, offsets):
                if xcoord < 0:
                    _x = -1
                    xmin, xmax = -1, xcoord
                    ha = 'left'
                else:
                    _x = 1
                    xmin, xmax = xcoord, 1 
                    ha = 'right'
                # line
                ax.hlines(liney, xmin=xmin, xmax=xmax, color=linecolor, linestyle='-')
                # text
                text_config = {
                    'xy': (_x, liney), 
                    'ha': ha, 'va': va, 
                    'xytext': offset,
                    'textcoords': 'offset points',
                    'fontsize': fontsize,
                }
                ax.annotate(text, **text_config)
                
            handles, labels = ax.get_legend_handles_labels()
            handles = [mpl.lines.Line2D([], [], c=h.get_edgecolor()) for h in handles]
            ax.legend(handles, labels, bbox_to_anchor=(1,1), loc='upper left')

    def plot_corr_bimodal(self, othr, dist_th=1e5):
        """
        """
        assert isinstance(othr, CorrRes)

        res1_corrs = self.res_corrs
        res2_corrs = othr.res_corrs
        res1_stats = self.res_stats
        res2_stats = othr.res_stats
        
        assert np.all(res1_corrs.index.values == res2_corrs.index.values)
        
        _x = res1_corrs[res1_corrs['dist']<dist_th]['corr'].values
        _y = res2_corrs[res2_corrs['dist']<dist_th]['corr'].values
        
        r1 = self.sign*res1_stats['r_th_linked']
        r2 = othr.sign*res2_stats['r_th_linked']

        sets = [
            set(res1_stats['linked_table'].index.values),
            set(res2_stats['linked_table'].index.values),
            ]
        num_sig_both = len(sets[0] & sets[1])
        num_sig1 = len(sets[0]) - num_sig_both
        num_sig2 = len(sets[1]) - num_sig_both

        # plot
        fig, ax = plt.subplots(figsize=(6, 6))
        g = ax.hexbin(_x, _y,  
                    gridsize=(100,100),
                    extent=(-1,1,-1,1),
                    cmap='rocket_r', 
                    bins='log', # log10(i+1)
                    rasterized=True,
                    )
        ax.axvline(r1, color='gray', linestyle='--', zorder=2)
        ax.axhline(r2, color='gray', linestyle='--', zorder=2)
        ax.set_aspect('equal')
        ax.set_xlim([-1,1])
        ax.set_ylim([-1,1])  
        ax.set_xticks([-1, r1, 1])
        ax.set_yticks([-1, r2, 1])
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter("%.2f"))
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.2f"))
        ax.set_xlabel(f'{self.label} correlation\n(Spearman)')
        ax.set_ylabel(f'{othr.label} correlation\n(Spearman)')
        ax.set_title(f'All enhancer-gene pairs (2-{dist_th/1000:.1f} kb)', pad=0)

        # annotate
        ax.annotate("{} linked\n(mCG-RNA only)".format(num_sig1), 
                    xy=(0,0), xycoords='axes fraction',
                    xytext=(0, -2*ax.xaxis.label.get_fontsize()), textcoords='offset points',
                    ha='right', va='top',
                    fontsize=ax.title.get_fontsize(),
                )
        ax.annotate("{} linked\n(ATAC-RNA only)".format(num_sig2), 
                    xy=(1,1), xycoords='axes fraction',
                    xytext=(0, +ax.title.get_fontsize()), textcoords='offset points',
                    ha='left', va='bottom',
                    fontsize=ax.title.get_fontsize(),
                )
        ax.annotate("{} linked\n(both)".format(num_sig_both), 
                    xy=(0,1), xycoords='axes fraction',
                    xytext=(0, +ax.title.get_fontsize()), textcoords='offset points',
                    ha='right', va='bottom',
                    fontsize=ax.title.get_fontsize(),
                )

        # lines
        ax.plot([-1.2,-0.9], [-1.2, -0.9], color='k', clip_on=False)
        ax.plot([ 1.1, 0.9], [ 1.1,  0.9], color='k', clip_on=False)
        ax.plot([-1.1,-0.9], [ 1.1,  0.9], color='k', clip_on=False)

        # cbar
        cbar = fig.colorbar(g, ax=ax, 
                            fraction=0.05, aspect=10,
                            label='Num. pairs per pixel\n(1/100 length)')
        cbar.ax.yaxis.set_major_formatter(mtick.EngFormatter())
        plt.show()