from os.path import join
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patheffects import Stroke
import matplotlib as mpl
import colorsys as cs
from sklearn.preprocessing import (MinMaxScaler,
                                   StandardScaler,
                                  )

def ContinousBlocks(df):
    """
    Identify continous blocks of data or gaps,
    returning a `pandas.Series` object with blocks of
    1 or 0 values indicating data or gaps,
    respectively.
    
    Parameters
    ----------
    df: pandas.DataFrame
       Dataframe with a column named 'n_day' in which
       the day number of the encounter is coded

    Returns
    -------
    cb: pandas.Series
    """
    
    cb = (pd
          .Series([1
                   if np.isnan(d)
                   else 0
                   for d
                   in df.day_n
                  ]
                 )
         )
    return cb
	
def MaxDays(df,
            final_buffer=42,
           ):
    """
    Identifies the intended length of the data 
    collection and adding a buffer of only 42
    days if the disengagement (last gap), occurred
    before said length. Returns a `pandas.DataFrame`
    object with the padded encounter dates.
    
    Parameters
    ----------
    df: pandas.DataFrame
       Dataframe with a column named 'encounter_date'
       in which the date of the encounter is coded
    final_buffer: int, default=42
        Number of days in which the data is extended
        should disengagement occurrs.

    Returns
    -------
    max_days: pandas.DataFrame
    """

    study_length = {'clintouch': 84,
                    'careloop': 84,
                    'empower': 365,
                   }
    OffSet = pd.DateOffset
    dataset = df.dataset.values[0]
    min_date, last_date = (df
                           .encounter_date
                           .sort_values()
                           .iloc[[0,-1]]
                          )
    theoretical_max = (min_date
                       +OffSet(days=study_length
                               [dataset]
                              )
                      )
    encounters_length = ((last_date-min_date)
                         .days
                        )
    max_date = (min(last_date+OffSet(final_buffer),
                    theoretical_max,
                   )
                if (encounters_length
                    <= study_length[dataset]
                   )
                else
                last_date
               )
    date_range = pd.date_range(min_date,
                               max_date,
                              )
    max_days = (pd
                .DataFrame({'encounter_date':date_range,
                           },
                          )
               )

    return max_days
	
def MissedEncounters(continous_blocks,
                     incremental_step=0.5,
                     min_scale_lim=1,
                     max_scale_lim=7,
                    ):
    """
    Find the continous blocks and calculates the
    `missed_encounter` synthetic construct, returning
    an array with these construct values.
    
    Parameters
    ----------
    continous_blocks: pandas.Series
       Array with blockes of either 1 or 0 values
       mapping the gaps and the data within the
       number of dates.
    incremental_step: float, default=0.5
        Increment per missing encounter
    min_scale_lim: int, default=1
        Default and minimum value of the array
    max_scale_lim: int, default=7
        Maximum value within the array

    Returns
    -------
    me: pandas.Series
    """    
    # expanding window calculations:
    # similar to rolling but the window
    # increases till including the whole array
    me = (continous_blocks
          .expanding()
          # if value of 1 leaves incremental
          # values otherwise it leaves the value
          # in array
          .apply(lambda s:
                 s.replace(to_replace=1,
                           method='bfill',
                          )
                 .sum()
                 *incremental_step
                )
           )

    # limits to 7 to match
    # Likert's scale
    me = me.map(lambda x: min_scale_lim
                if x<incremental_step
                else min(x,
                         max_scale_lim-1
                        )+1,
               )
    return me
	
def Disengagement(missed_encounters):
    """
    Calculation of the values of the `disengagement`
    synthetic construct, returning an array with 
    binary values of either 1 or 7.
    
    Parameters
    ----------
    missed_encounters: pandas.Series
        Array with values of `missed_encounters`
        within the range of 1 to 7.

    Returns
    -------
    d: pandas.Series
    """        
    d = (missed_encounters.map(lambda x:
                                 1 if x<7
                                 else 7
                                )
                          )

    return d
    
def GraphBuilding(nodes, 
                  edges,
                  directional=False,
                 ):
    """
    Production of networkx graph objects based on the 
    CSV files resulting of the R-mlVAR process.

    Parameters
    ----------
    nodes: dict
       Nested dictionary containing the nodes' labels
       as keys in the first level and the attributes
       of the node in a second level dictionary. 
    edges: dict
       Nested dictionary containing the edges' labels 
       as keys in the first level, and the attributes 
       of the node in a second level dictionary. 
    directional: bool
       To define whether the network is either 
       directional or adirectional.

    Returns
    -------
    g: networkx graph object
    """
    if not directional:
        g = nx.Graph()
    else:
        g = nx.DiGraph()
    g.add_nodes_from(node for node in nodes.items())
    g.add_edges_from([(c0, 
                       c1,
                       d,
                      )
                      for (c0, c1), d 
                      in edges.items()
                     ],
                    )
    return g
    
def StrengthCentrality(G, node, attribute):
    """
    Calculation of node strength centrality, which is
    defined as the sum of the absolute values of 
    weights (or strengths) of all of the nodes edges
    connected to the given node.

    Parameters
    ----------
    G: graph
       A networkx graph
    node: object
       A node in the networkx graph
    attribute: object
       The edge attribute used to quantify node strength

    Returns
    -------
    output: float
    """
    output = 0.0
    for edge in G.edges(node):
        output += np.abs(G.get_edge_data(*edge)[attribute])
    return output
    
def mlVARNetworkPlots(networks,
                      out_folder=r'.',
                     ):
    """
    Plotting of the three networks produced by the 
    R-mlVAR process.

    Parameters
    ----------
    networks: dict
       A networkx graph
    node: object
       A node in the networkx graph
    attribute: object
       The edge attribute used to quantify node strength

    Returns
    -------
    PNG and SVG plot files.
    """
    
    mosaic=[[0,1,2]]
    fig, axes = plt.subplot_mosaic(mosaic=mosaic,
                                   figsize=[30,11],
                                  )

    for i, (g_type, g) in enumerate(networks):
        for k, v in {'_':'\n',
                    }.items():
            g = ModifyNodesNames(g,
                                 k,
                                 v,
                                )
        GraphPlot(g,
                  g_type,
                  g_type.title(),
                  axes,
                  fig,
                  i,
                 )

    fig_legend = []
    for n, c in g.nodes.data('color'):
        fig_legend.append(mpl
                          .lines
                          .Line2D([], [],
                                  color=c,
                                  marker='o',
                                  lw=0.0,
                                  markersize=35,
                                  label=n.replace('_', ' '),
                                 )
                         )
    fig.tight_layout(h_pad=0.0,
                     w_pad=-1.0,
                    )

    filename = '3Networks_AllData_Graphs'
    for ext in ['png', 'svg'][:]: 
        out_path = join(out_folder,
                        f'{filename}.{ext}',
                       )
        plt.savefig(out_path,
                    bbox_inches='tight',
                    dpi=300,
                   )
                   
def ModifyNodesNames(g,
                     k,
                     v,
                    ):
    """
    Changes the label of a node within the networkx graph
    object, among nodes and edges.
    
    Parameters
    ----------
    g: networkx graph object
        A networkx graph
    k: str
       Label of the node to modify
    v: str
       Modified label of the node

    Returns
    -------
    h: networkx graph object with modified node label.
    """
    
    if 'DiGraph' not in str(g):
        h = nx.Graph()
    else:
        h = nx.DiGraph()

    h.add_nodes_from([(n
                       .replace(k, v),
                       d)
                      for n, d
                      in sorted(g
                                .nodes(data=True)
                               )
                     ]
                    )
    h.add_edges_from([(n0
                       .replace(k, v),
                       n1
                       .replace(k, v),
                       d)
                      for n0, n1, d
                      in g.edges(data=True)
                     ]
                    )
    return h
    
def GraphPlot(g,
              g_type,
              title,
              axes,
              fig,
              i,
              top_pct=90,
             ):
    '''
    Plotting of networks based on a Networkx Graph object.

    Parameters
    ----------
    g: graph
       A networkx graph.
    g_type: str
       Type of graph, either directed or indirected.
    title: str
       Title of the plot
    axes: list of matplotlib axes
        In case of mutiple networks in a plot.
    i: int
        index for the location of the plot in a
        list of axes.
    top_pct: int
            Integer between 0 and 100 limiting the
            edges to plot to a specified top
            percentile.

    Returns
    -------
    `matplotlib.axes.Axes`: plot of a graph.
    '''
    MMS_ewidth = (MinMaxScaler(feature_range=(0.0,
                                             20.0,
                                            )
                             )
                 .fit_transform
                )
    MMS_ealpha = (MinMaxScaler(feature_range=(0.1,
                                             0.8,
                                            )
                             )
                 .fit_transform
                )
    MMS_nsize = (MinMaxScaler(feature_range=(8_000,
                                             24_000,
                                            )
                             )
                 .fit_transform
                )
    print(g)
    if 'DiGraph' not in str(g):
        arrowstyle='-'
        connectionstyle=('arc3'
                         ',rad=-0.125'
                        )
    else:
        arrowstyle='-|>'
        connectionstyle=('arc3'
                         ',rad=0.125'
                        )
    ax = axes[i]
    ax.set_aspect(1)
    ax.axis("off")
    ax.margins(0.1)
    # ax.set_gid(f'id_{i}')
    pos = nx.circular_layout(g)

    dict_edges = {(n0,n1):s
                 for n0, n1, s
                 in (g
                     .edges
                     .data('similarity')
                    )
                }
    edge_similarity = np.fromiter(dict_edges.values(),
                                    dtype=float
                                   )
    edges_pctil = EdgeWeightsPctil(g,
                                   top_pct=top_pct,
                                  )
    loops = {(n0,n1):v
             for (n0,n1), v
             in dict_edges.items()
             if n0==n1
            }
    nonloops = list(zip(*{(n0,n1):v
                     for (n0,n1), v
                     in dict_edges.items()
                     if (n0,n1) not in loops.keys()
                    }.items()
                  ))
    loops = list(zip(*loops.items()))
    dict_width = {}
    dict_alpha = {}
    for items in [nonloops, loops]:
        if len(items) != 0:
            keys, values = items
            values = np.fromiter(values,
                                 dtype=float
                                )
            edge_width = (MMS_ewidth(values
                                     .reshape(-1,1),
                                    )
                          .flatten()
                         )
            [dict_width.update({k:v})
             for k,v
             in zip(keys, edge_width)
            ]
            edge_alpha = (MMS_ealpha(values
                                     .reshape(-1,1),
                                    )
                          .flatten()
                         )
            [dict_alpha.update({k:v})
             for k,v
             in zip(keys, edge_alpha)
            ]

    edge_width = list({k:dict_width[k]
                  for k
                  in dict_edges.keys()
                 }.values())
    edge_alpha = list({k:(dict_alpha[k]
                          if k in edges_pctil
                          or k[0]==k[1]
                          else 0.0
                         )
                       for k
                       in dict_edges.keys()
                      }.values())

    size_values = np.array([v
                            for _, v
                            in (g
                                .nodes
                                .data('strength')
                               )
                           ]
                          )
    node_size = (MMS_nsize(size_values
                           .reshape(-1,1)
                         )
                 .flatten()
                )
    node_color = [spectrumHSL(i/len(list(g)))
                  for i, n
                  in enumerate(list(g))
                 ]

    edge_color = [[0.0]*3 if s > 0
                  else [0.5,0.0,0.0]
                  for s
                  in edge_similarity
                 ]
    nx.draw_networkx(g,
                     pos=pos,
                     node_size=node_size,
                     nodelist=list(g),
                     with_labels=True,
                     node_color=[1.0]*4,
                     linewidths=8.0,
                     edgecolors=node_color,
                     width=edge_width,
                     edge_color=[c+[a]
                                 for c,a
                                 in zip(edge_color,
                                        edge_alpha,
                                       )
                                ],
                     arrows=True,
                     arrowstyle=arrowstyle,
                     connectionstyle=connectionstyle,
                     ax=ax,
                    )
    [patch
     .set_path_effects([Stroke(joinstyle='miter',
                               capstyle='round',
                              )
                       ]
                      )
     for patch in ax.patches
    ]

    node_dict = {k:v
             for k,v
             in zip(list(g),
                    node_size
                   )
            }

    [t.set_fontsize((node_dict[t.get_text()]
                       ** 0.5
                       )
                       * 1.8
                      / len(t.get_text()
                            .split('\n')
                            [-1]
                           )
                     )
     for t in ax.texts
     if t.get_text() in list(g)
    ]
    transAxFig = (mpl
                  .transforms
                  .blended_transform_factory(ax.transAxes,
                                             fig.transFigure,
                                            )
                 )
    ax.set_title(title,
                 fontsize=40,
                 x=0.5,
                 y=0.975,
                 ha='center',
                 va='top',
                 transform=transAxFig
                )
    if g_type == 'temporal':
        loops = {n0:i
                 for i, (n0, n1)
                 in enumerate(g.edges())
                 if n0==n1
                }
        node_size = {n:s
                     for n,s
                     in zip(list(g),
                            node_size,
                           )
                    }
        for node, index in loops.items():
            loop = ax.patches[index]
            alpha = loop.get_edgecolor()[-1]
            loop.set_color([0.1,0.1,0.667,alpha])
            trans_y=0.0
            if 5_000<node_size[node]<=16_000:
                trans_y=0.1
            elif node_size[node]>16_000:
                trans_y=0.2
            t2 = (mpl
                  .transforms
                  .Affine2D()
                  .translate(0,trans_y)
                  + ax.transData
                 )
            loop.set_transform(t2)
    ax.set_gid(f'id_{i}')
    
def EdgeWeightsPctil(g,
                     value='similarity',
                     top_pct=50,
                    ):
    """
    Filters the edges to show only those which
    absolute weight values represents a specified
    top percentile. Returning a list of the edges
    to plot.
    
    Parameters
    ----------
    g: graph
       A networkx graph.
    value: str
       The label of the edges' weights.
    top_pct: float
       Top percentile that will be shown represented
       by the edges' weights.

    Returns
    -------
    list: Edges to plot
    
    """

    dict_edges = {(n0,n1):[np.abs(s)]
                  for n0, n1, s
                  in (g
                      .edges
                      .data(value)
                     )
                  if n0 != n1
                 }

    df = (pd
          .DataFrame(dict_edges)
          .T
         )
    pctil = (df
             .sum()
             * (1
                - (top_pct
                   / 100)
               )
            )
    df = (df
          .sort_values(by=0,
                      )
          .reset_index()
         )
    cumsum_col = 'cumsum'
    df[cumsum_col] = df.iloc[:,2].cumsum()
    idx = df[cumsum_col].searchsorted(pctil)[0]
    edges_pctil = (df
                   .iloc[idx:,:]
                   [::-1]
                   .set_index(['level_0',
                               'level_1',
                              ],
                             )
                   .iloc[:,0]
                   .to_dict()
                   .keys()
                  )

    return list(edges_pctil)
    
def SpectrumHSL(segments=11,
             name='spectrumHSL',
             hue_range=[0.0,
                        1.0,
                       ],
             sat_range=[1.00,
                        0.67,
                       ],
             lig_range=[0.1,
                        0.67,
                       ],
            ):
    '''
    Produces a colormap attempting to minimise the
    green region to have a more diverse colour palette,
    as well as a progressive reduction of the lightness, to
    differentiate the beginning and the end of the
    colormap. 
    
    Parameters
    -----------
    segments : Int, Default : 11
        Segments in which the hue space will be divided
    name : String, Default : 'spectrum'
        Name for the new colormap. 
    hue_range : Array-like, default: [0.0,1.0]
        Range of hue from which the new colormap
        will be sampled
    sat_range : Array-like, default: [0.67,0.75]
        Range of saturation from which the new
        colormap will be sampled
    lig_dange: Array-like, default: [0.5,0.9]
        Range of lightness from which the new 
        colormap will be sampled
    
    
    Returns
    -------
    `matplotlib.colors.ListedColormap`
    '''
            
    colors = mpl.colors
    lspace = np.linspace
    rgb_list = [cs
                .hls_to_rgb(h,l,s)
                             for h,s,l
                             in zip(lspace(*hue_range,
                                           segments,
                                          ),                                         
                                    lspace(*sat_range,
                                           segments
                                          ),
                                    lspace(*lig_range,
                                           segments
                                          ), 

                                   )
                            ]

    rgb_list = np.delete(rgb_list,
                         (segments//2)-1,
                         axis=0,
                        )
    spectrum = (colors
                .LinearSegmentedColormap
                .from_list(name=name,
                           colors=rgb_list,
                          )
               )

    return spectrum

spectrumHSL=SpectrumHSL()

def CSV2NetworkxGraphs(dict_mlvar):
    """Production of dictionary with
    `networkx.classes.graph.Graph` objects using CSV
    files resulting from R-mlVAR function.
    
    Parameters
    ----------
    dict_mlvar: dict
        Dictionary in which the keys are the type of 
        mlVAR (temporal, between or contamporaneous).
        These are the raw result of the mlVAR R library, 
        extracted from R object:
        `mlVAR$summary()[[<type of network>]]`
    
    Returns
    -------
    dict_graphs: dict
    
    """
    dict_graphs = {}
    for key, df_mlvar in dict_mlvar.items():
        constructs = sorted(df_mlvar
                      .iloc[:,:2]
                      .stack()
                      .value_counts()
                      .keys()
                     )
        constructs_combinations = df_mlvar.values

        edges = {(c0,c1):float(v)
                 for c0,c1,v
                 in constructs_combinations
                }

        nodes = {q:{}
                 for i, q 
                 in enumerate(constructs)
                }

        fontsize_node = 9.5

        max_distance = max([d for d in edges.values()])

        edges = {k:{'similarity':v,
                    'absolute_similarity': abs(v),
                   }
                 for k,v
                 in edges.items()
                }
        directional = (True
                       if key=='temporal'
                       else
                       False
                      )
        g = GraphBuilding(nodes,
                          edges,
                          directional=directional,
                         )

        betweenness = nx.betweenness_centrality(g,
                                                endpoints=True,
                                                weight='absolute_similarity',
                                               )
        closeness = nx.closeness_centrality(g,
                                            distance='absolute_similarity',
                                            wf_improved=True,
                                           )

        strength = {node:StrengthCentrality(g,
                                            node,
                                            'absolute_similarity',
                                           )
                    for node
                    in g.nodes()
                   }
        nodes_centralities = {n:{name:c[n] 
                                 for name, c
                                 in [['betweenness', betweenness],
                                     ['strength', strength],
                                     ['closeness', closeness],
                                    ]
                                } for n 
                              in nodes.keys()
                             }
        g = GraphBuilding(nodes_centralities,
                          edges,
                          directional=(True
                                       if key=='temporal' 
                                       else
                                       False
                                      ),
                         )
        key = (key
               if key!='between'
               else
               f'{key}-subjects'
              )
        dict_graphs[key] = g
    return dict_graphs