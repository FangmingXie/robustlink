"""Import commonly used libraries"""

import numpy as np
import pandas as pd
import collections

# matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['pdf.fonttype'] = 42 # editable text in matplotlib
mpl.rcParams['svg.fonttype'] = 'none'
# mpl.rcParams['font.family'] = 'sans-serif'
# mpl.rcParams['font.sans-serif'] = ['Arial']

import matplotlib.ticker as mtick
# PercentFormat = mtick.FuncFormatter(lambda y, _: '{:.1%}'.format(y))
# ScalarFormat = mtick.ScalarFormatter()

# seaborn
import seaborn as sns
sns.set_style('ticks', rc={'axes.grid':True})
sns.set_context('talk')

# # set matplotlib formats
# from IPython.display import set_matplotlib_formats
# set_matplotlib_formats('retina')

# empty rectangle (for legend)
EMPTY_RECTANGLE = mpl.patches.Rectangle((0,0), 1, 1, fill=False, edgecolor='none',
		                                 visible=False)

DEFAULT_COLORBAR_KWS = {'fraction': 0.05, 
						'shrink': 0.4, 
						'aspect': 5, 
						}
