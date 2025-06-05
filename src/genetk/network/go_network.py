#!/usr/bin/env python

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.colors as mcolors
from collections import defaultdict
from itertools import combinations

class GONetwork:
    def __init__(self, enrichr_file, adj_pval_threshold=0.05):
        """
        Initialize GO network from Enrichr results file.
        
        Parameters:
        enrichr_file (str): Path to Enrichr results file
        adj_pval_threshold (float): Adjusted P-value threshold for filtering terms (default: 0.05)
        """
        self.enrichr_file = enrichr_file
        self.adj_pval_threshold = adj_pval_threshold
        self.data = self._load_data()
        self.G = self._build_network()
        
    def _load_data(self):
        """Load and process Enrichr results data."""
        data = pd.read_csv(self.enrichr_file, sep='\t')
        data = data.dropna(subset=['Genes'])
        data = data[data['Adjusted P-value'] <= self.adj_pval_threshold]
        data['Genes_list'] = data['Genes'].apply(lambda x: x.split(';'))
        return data
        
    def _jaccard_index(self, set1, set2):
        """Calculate Jaccard index between two gene sets."""
        intersection = len(set1.intersection(set2))
        union = len(set1.union(set2))
        return intersection / union if union > 0 else 0
        
    def _build_network(self):
        """Build networkx graph with Terms as nodes and Jaccard index as edge weights."""
        G = nx.Graph()
        
        # Add nodes with attributes
        for idx, row in self.data.iterrows():
            term = row['Term']
            G.add_node(term,
                      gene_set=row['Gene_set'],
                      overlap=row['Overlap'],
                      p_value=row['P-value'],
                      adjusted_p_value=row['Adjusted P-value'],
                      odds_ratio=row['Odds Ratio'],
                      combined_score=row['Combined Score'],
                      genes=set(row['Genes_list']))
        
        # Add edges with Jaccard index as weight
        terms = list(G.nodes())
        for term1, term2 in combinations(terms, 2):
            genes1 = G.nodes[term1]['genes']
            genes2 = G.nodes[term2]['genes']
            jaccard = self._jaccard_index(genes1, genes2)
            if jaccard > 0:
                G.add_edge(term1, term2, weight=jaccard)
        
        # Calculate degree weighted by Jaccard index
        for node in G.nodes():
            weighted_degree = G.degree(node, weight='weight')
            G.nodes[node]['weighted_degree'] = weighted_degree
            
        return G
    
    def filter_nodes_by_attribute(self, attribute='odds_ratio', threshold=None):
        """
        Filter nodes by attribute value.
        
        Parameters:
        attribute (str): Node attribute to filter by (default: 'odds_ratio')
        threshold (float): Threshold value (default: 25th percentile)
        
        Returns:
        GONetwork: New GONetwork instance with filtered nodes
        """
        if threshold is None:
            values = [self.G.nodes[node][attribute] for node in self.G.nodes()]
            threshold = np.percentile(values, 25)
        
        nodes_to_keep = [node for node in self.G.nodes() 
                        if self.G.nodes[node][attribute] >= threshold]
        
        filtered_G = self.G.subgraph(nodes_to_keep).copy()
        
        # Create new GONetwork instance with filtered graph
        new_network = GONetwork.__new__(GONetwork)
        new_network.enrichr_file = self.enrichr_file
        new_network.data = self.data[self.data['Term'].isin(nodes_to_keep)]
        new_network.G = filtered_G
        
        return new_network
    
    def filter_edges_by_weight(self, threshold=None):
        """
        Filter edges by weight.
        
        Parameters:
        threshold (float): Weight threshold (default: 25th percentile)
        
        Returns:
        GONetwork: New GONetwork instance with filtered edges
        """
        if threshold is None:
            weights = [data['weight'] for _, _, data in self.G.edges(data=True)]
            threshold = np.percentile(weights, 25)
        
        filtered_G = self.G.copy()
        edges_to_remove = [(u, v) for u, v, data in filtered_G.edges(data=True) 
                          if data['weight'] < threshold]
        filtered_G.remove_edges_from(edges_to_remove)
        
        # Create new GONetwork instance with filtered edges
        new_network = GONetwork.__new__(GONetwork)
        new_network.enrichr_file = self.enrichr_file
        new_network.data = self.data
        new_network.G = filtered_G
        
        return new_network
    
    def community(self, weight='weight'):
        """Find communities using Louvain algorithm."""
        lpc = nx.community.louvain_communities(self.G, weight=weight)
        community_index = {n: i for i, com in enumerate(lpc) for n in com}
        return community_index
    
    def bubble(self, name, color=None, community_colored=False, 
               text_annot=False, term_list=None, specify_color=None, 
               size_attribute='odds_ratio', color_attribute=None):
        """
        Create bubble plot visualization of GO network.
        
        Parameters:
        name (str): Name for output file
        color: Default color for nodes
        community_colored (bool): Color nodes by community
        text_annot (bool): Show text annotations
        term_list (list): List of terms to annotate
        specify_color (dict): Dictionary mapping nodes to colors
        size_attribute (str): Node attribute for sizing (default: 'odds_ratio')
        color_attribute (str): Node attribute for coloring (e.g., 'adjusted_p_value')
        """
        pos = nx.spring_layout(self.G, k=0.15, seed=4572321)
        
        # Node sizes based on specified attribute
        node_size = []
        for node in self.G.nodes():
            try:
                size_val = self.G.nodes[node][size_attribute]
                node_size.append(size_val * 10)
            except:
                node_size.append(10)
        
        # Node colors
        if community_colored:
            community_index = self.community(weight='weight')
            community_data = pd.DataFrame.from_dict(community_index, orient='index')
            community_data.to_csv(f'{name}_community.tsv', sep='\t', header=False, index=True)
            node_color = [community_index[n] for n in self.G.nodes()]
        elif specify_color:
            node_color = [specify_color.get(node, 'gray') for node in self.G.nodes()]
        elif color_attribute:
            # Color by -log10(adjusted p-value) or other attribute
            if color_attribute == 'adjusted_p_value':
                node_color = [-np.log10(self.G.nodes[node][color_attribute]) 
                             for node in self.G.nodes()]
            else:
                node_color = [self.G.nodes[node][color_attribute] 
                             for node in self.G.nodes()]
        else:
            node_color = [color] * len(self.G.nodes()) if color else ['blue'] * len(self.G.nodes())
        
        # Terms to annotate
        if term_list is None and text_annot:
            # Use top terms by weighted degree
            centrality = dict(self.G.degree(weight='weight'))
            terms_to_annotate = sorted(centrality, key=centrality.get, reverse=True)[:10]
        else:
            terms_to_annotate = term_list if term_list else []
        
        self.draw_params(pos, node_color, node_size, terms_to_annotate, 
                        color, name, text_annot, color_attribute, size_attribute)
    
    def draw_params(self, pos, node_color, node_size, terms_to_annotate, 
                   color, name, text_annot, color_attribute=None, size_attribute='odds_ratio'):
        """Draw network with specified parameters."""
        fig, ax = plt.subplots(figsize=(20, 15))
        
        if color_attribute == 'adjusted_p_value':
            # Use coolwarm colormap for p-values
            nx.draw_networkx(
                self.G,
                pos=pos,
                cmap=plt.get_cmap('coolwarm'),
                with_labels=False,
                node_color=node_color,
                node_size=node_size,
                edge_color="gainsboro",
                alpha=0.8,
                vmin=min(node_color) if node_color else 0,
                vmax=max(node_color) if node_color else 1
            )
            sm = plt.cm.ScalarMappable(cmap=plt.get_cmap('coolwarm'), 
                                     norm=plt.Normalize(vmin=min(node_color) if node_color else 0, 
                                                       vmax=max(node_color) if node_color else 1))
            sm._A = []
            cbar = plt.colorbar(sm, ax=ax)
            cbar.set_label('-log10(Adjusted P-value)', fontsize=14)
        else:
            # Regular drawing
            nx.draw_networkx(
                self.G,
                pos=pos,
                with_labels=False,
                node_color=node_color,
                node_size=node_size,
                edge_color="gainsboro",
                alpha=0.8
            )
        
        # Text annotations
        if text_annot and terms_to_annotate:
            for term in terms_to_annotate:
                if term in pos:
                    x, y = pos[term]
                    plt.text(x, y, term, fontsize=8, color=color or 'red', 
                            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7))
        
        # Title and legend
        font = {"color": "k", "fontweight": "bold", "fontsize": 20}
        ax.set_title(f"GO Network Analysis - {name}", font)
        
        # Legend text
        legend_color = color or "r"
        font["color"] = legend_color
        
        ax.text(
            0.80,
            0.10,
            "edge weight = Jaccard index",
            horizontalalignment="center",
            transform=ax.transAxes,
            fontdict=font,
        )
        ax.text(
            0.80,
            0.06,
            f"node size = {size_attribute}",
            horizontalalignment="center",
            transform=ax.transAxes,
            fontdict=font,
        )
        
        ax.margins(0.1, 0.05)
        fig.tight_layout()
        plt.axis("off")
        plt.savefig(f'figures/{name}_go_network.pdf', dpi=300, bbox_inches='tight')
        plt.show()