# Genetk

A toolkit for genetic analysis.

## Installation

```bash
uv pip install git+https://github.com/brainfo/genetk.git
```

## Features

- Use gseapy for gene set enrichment analysis
  - given gtf file, annotate and select only protein-coding genes for enrichment
  - The gseapy plotting utils support mapping q values to colors for dot plots but not barplots. They have barplots with fixed float (p or q values) as one axis, and at most color by database groups.
    - thus, here's a simple barplot function for flexibly choosing field as axis and mapping colors
  - Support both explorary saving out all sig terms directly from fresh analysis, or pre-saved tables
    - for pre-saved tables, the plot saving is boosted with plutils package
  - from the output gsea enriched table, use networkx for network analysis for the Terms. the edge is weighted by the Jaccard index of the Genes between Terms.
- Use magma to enrich for GWAS hits
- Use networkx for Network analysis given any edge list, e.g. from PPI analysis
  - In-degree, Out-degree, and Betweenness centrality
  - Find communities

## Example Usage

1. for enrichment analysis given DEGs

   Input: an excel of all DEG results with sheetnames, output **all** GO ORA and GSEA results tables and visualize **all** significant terms.

   default using databases 'MSigDB_Hallmark_2020', 'KEGG_2021_Human', and 'GO_Biological_Process_2023'
   
   ```python
   from genetk.enrich.go_enrichment import process_excel_sheets
process_excel_sheets("test/test.xlsx", gtf_file="test/genes.gtf", run_id="trial_001")
   ```

    Output:
      - ORA enrichment tables and one barplot for all enriched terms from given databases
     - GSEA enrichment tables from each database, GSEA plots with each significant term, one GSEA plot with all significant terms
     - All outputs are stored under `data/enrichr/<run_id>` and `figures/enrichr/<run_id>` so consecutive runs never clobber each other

2. Barplot of enriched terms (likely from 1. or any other dataframe)

   ```python
   import pandas as pd
   import numpy as np

   test_file = "test/GO_Biological_Process_2023.human.enrichr.reports.txt" 
   ```

   ```python
   df = pd.read_csv(test_file, sep='\t')
   genetk.enrich.go.barh(df, 'GO', 'down', 'exact', save = 'test.pdf', return_ax=False)
   ```


3. network analysis (Term-level) of enriched terms

   ```python
   from genetk.network.go_network import GONetwork
   go_net = GONetwork(test_file)

   filtered_net = go_net.filter_nodes_by_attribute('odds_ratio')
   edge_filtered_net = filtered_net.filter_edges_by_weight()

   communities = edge_filtered_net.community(resolution=0.5)
   num_communities = len(set(communities.values()))

   os.makedirs('figures', exist_ok=True)

   edge_filtered_net.bubble(
   name='test_community',
   community_colored=True,
   text_annot=True
   )

   ## output figure test/figures/test_community_go_network.pdf
   ```

