## plot
## from pre-made enrichr tables

import pandas as pd
import numpy as np
import seaborn as sns
import __future__
import matplotlib.pyplot as plt
import plutils as plu

def plu_savefig(ax, save, plot_size=None, legend_out=True, right_pad_frac=0.3, mplrc=None):
    plu.reset_mpl_style(mplrc)
    plu.save_single_ax(ax, save, plot_size=plot_size, legend_out=True, right_pad_frac=0.3)

def make_plot_data(df, term, overlap):
    """
    Prepare data for plotting. Returns DataFrame with columns: Term, x, -logp
    x is ratio or n_overlap depending on overlap param.
    """
    if isinstance(df, str):
        df = pd.read_csv(df, sep='\t', decimal='.')
    else:
        df = df.copy()
    df['Adjusted P-value'] = df['Adjusted P-value'].astype(float)

    if term == 'MP':
        term_name = [' '.join(segment.split(' ')[:-1]) for segment in df.Term]
    else:
        term_name = [segment.split('(')[0] for segment in df.Term]

    ratio = df.Overlap.apply(lambda x: eval(compile(x, '<string>', 'eval', __future__.division.compiler_flag)))
    n_overlap = [int(overlap.split('/')[0]) for overlap in df.Overlap]
    bonferroni = -np.log(pd.to_numeric(df['Adjusted P-value'], errors='coerce'))

    if overlap == 'exact':
        x = n_overlap
    else:
        x = ratio

    return pd.DataFrame({'Term': term_name, 'x': x, '-logp': bonferroni})


def barh(df, term, direction, overlap, save, return_ax=False, mplrc=None, plot_size=(0.443125, 0.418125)):
    df_toplot = make_plot_data(df, term, overlap)
    if direction == 'up':
        cmap = 'YlOrRd'
    else:
        cmap = 'YlGnBu'
    norm = plt.Normalize(df_toplot['-logp'].min(), df_toplot['-logp'].max())
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    ax = sns.barplot(x='x', y='Term', hue='-logp', data=df_toplot, palette=cmap, dodge=False, legend=False)
    sub_ax = plt.axes([0.96, 0.55, 0.02, 0.3])  # add a small custom axis
    ax.figure.colorbar(sm, cax=sub_ax)
    if return_ax:
        return ax
    if save:
        plu_savefig(ax, save, plot_size=plot_size, legend_out=True, right_pad_frac=0.3, mplrc=mplrc)


def scatter_okot(df, term, direction, overlap, save=None, return_ax=False, mplrc=None, plot_size=(0.443125, 0.418125)):
    """
    Scatter (dot) plot version of barh, using same data logic.
    """
    df_toplot = make_plot_data(df, term, overlap)
    if direction == 'up':
        cmap = 'YlOrRd'
    else:
        cmap = 'YlGnBu'
    norm = plt.Normalize(df_toplot['-logp'].min(), df_toplot['-logp'].max())
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    ax = sns.scatterplot(x='x', y='Term', hue='-logp', data=df_toplot, palette=cmap, legend=False, s=80)
    sub_ax = plt.axes([0.96, 0.55, 0.02, 0.3])  # add a small custom axis
    ax.figure.colorbar(sm, cax=sub_ax)
    if return_ax:
        return ax
    if save:
        plu_savefig(ax, save, plot_size=plot_size, legend_out=True, right_pad_frac=0.3, mplrc=mplrc)