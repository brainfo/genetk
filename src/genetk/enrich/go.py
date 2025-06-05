## plot
## from pre-made enrichr tables

import pandas as pd
import numpy as np
import seaborn as sns
import __future__
import matplotlib.pyplot as plt
# import plutils as plu

# def plu_savefig(ax, save, plot_size=None, legend_out=True, right_pad_frac=0.3, mplrc=None):
#     plu.reset_mpl_style(mplrc)
#     plu.save_single_ax(ax, save, plot_size=plot_size, legend_out=True, right_pad_frac=0.3)

def barh(df, term, direction, overlap, save, return_ax=False, mplrc=None, plot_size=(0.443125, 0.418125)):
    if isinstance(df, str):
        df = pd.read_csv(df,sep='\t',decimal='.')
    else:
        df = df.copy()
    df['Adjusted P-value'].astype(float)

    if term == 'MP':
        term_name = [' '.join(segment.split(' ')[:-1]) for segment in df.Term]
    else:
        term_name = [segment.split('(')[0] for segment in df.Term]

    ratio = df.Overlap.apply(lambda x: eval(compile(x, '<string>', 'eval', __future__.division.compiler_flag)))
    n_overlap = [int(overlap.split('/')[0]) for overlap in df.Overlap]

    bonferroni = -np.log(pd.to_numeric(df['Adjusted P-value'], errors='coerce'))

    if direction == 'up':
        cmap='YlOrRd'
    else:
        cmap='YlGnBu'
    if overlap=='exact':
        d = {'Term':term_name, 'Overlap':n_overlap, '-logp':bonferroni}
    else:
        d = {'Term':term_name, 'Overlap':ratio, '-logp':bonferroni}

    df_toplot = pd.DataFrame(data=d)
    df_toplot = df_toplot.sort_values(by='-logp', ascending=False)
    norm = plt.Normalize(bonferroni.min(), bonferroni.max())
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    ax = sns.barplot(x='Overlap', y='Term', hue='-logp', data=df_toplot, palette=cmap, dodge=False, legend=False)
    sub_ax = plt.axes([0.96, 0.55, 0.02, 0.3]) # add a small custom axis
    ax.figure.colorbar(sm, cax=sub_ax)
    if return_ax:
        return ax
    # if save:
    #     plu_savefig(ax, save, plot_size=plot_size, legend_out=True, right_pad_frac=0.3, mplrc=mplrc)