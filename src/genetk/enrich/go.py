## plot
## from pre-made enrichr tables

import pandas as pd
import numpy as np
import seaborn as sns
import __future__
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = ['sans-serif']
matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams.update({'font.size': 18})

# def plot_GO(fname, name, term):
    
#     df = pd.read_csv(fname,sep='\t')
#     df['Adjusted P-value'].astype(float)
#     df = df.loc[df['Adjusted P-value']<0.05]
#     if term == 'MP':
#         term_name = [' '.join(segment.split(' ')[:-1]) for segment in df.Term]
#     else:
#         term_name = [segment.split('(')[0] for segment in df.Term]
#     term_name.reverse()
#     count = df.Overlap.tolist()
#     count.reverse()
#     bonferroni = -np.log(df['Adjusted P-value'])
#     bonferroni = bonferroni.tolist()
#     bonferroni.reverse()
#     fig = make_subplots(
#         rows=1, cols=2,
#         shared_xaxes=True,
#         vertical_spacing=0.03,
#         specs=[[{"type": "bar"},{"type": "scatter"}]]
#     )

#     fig.add_trace(go.Bar(
#         y=term_name,
#         x=count,
#         name='Count',
#         orientation='h'
#     ))
    
#     fig.add_trace(
#         go.Scatter(
#             x=bonferroni,
#             y=term_name,
#             name="-log(padj)",
#             line=dict(dash='dot')
#         ),
#         row=1, col=2
#     )


#     fig.update_layout(template='simple_white',title_text=name,xaxis = dict(showgrid=False, ticks='inside',mirror=False,showline=True),
#                 yaxis = dict(showgrid=False, ticks='inside',mirror=False,showline=True),yaxis2=dict(side = 'right',showticklabels=False),
#                 font=dict(size=18),
#                   width=2000,
#                   height=800,
#                  title_x=0.5,
#                  showlegend=True)
#     fig.write_image(f'figures/{name}.pdf', scale=1, engine="kaleido", font_family = "Arial")

def plot_GO_one(fname, name, term):
    
    df = pd.read_csv(fname,sep='\t',decimal='.')
    df['Adjusted P-value'].astype(float)
    df = df.loc[df['Adjusted P-value']<0.05]
    if term == 'MP':
        term_name = [' '.join(segment.split(' ')[:-1]) for segment in df.Term]
    else:
        term_name = [segment.split('(')[0] for segment in df.Term]
    #term_name.reverse()
    ratio = df.Overlap.apply(lambda x: eval(compile(x, '<string>', 'eval', __future__.division.compiler_flag)))
    #count = pd.to_numeric(df.Overlap, errors='coerce').tolist()
    #count.reverse()
    bonferroni = -np.log(pd.to_numeric(df['Adjusted P-value'], errors='coerce'))
    #bonferroni.reverse()
    
    d = {'Term':term_name, 'Overlap':ratio, '-logp':bonferroni}
    df_toplot = pd.DataFrame(data=d)
    norm = plt.Normalize(bonferroni.min(), bonferroni.max())
    sm = plt.cm.ScalarMappable(cmap="YlOrRd", norm=norm)
    sm.set_array([])
    fig=plt.figure(figsize=(4, 3))
    ax = sns.barplot(x='Overlap', y='Term', hue='-logp', data=df_toplot, palette='YlOrRd', dodge=False)
    #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.get_legend().remove()
    fig.subplots_adjust(right=0.95) # create space on the right hand side
    sub_ax = plt.axes([0.96, 0.55, 0.02, 0.3]) # add a small custom axis
    ax.figure.colorbar(sm, cax=sub_ax)
    plt.savefig(f'figures/{name}.pdf', bbox_inches='tight', dpi=900, transparent=True)