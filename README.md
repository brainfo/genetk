# Genetk

A toolkit for genetic analysis.

## Installation

```bash
pip install git+https://github.com/brainfo/plutils.git
pip install git+https://github.com/brainfo/genetk.git
```

## Features

- Use gseapy for gene set enrichment analysis
  - given gtf file, annotate and select only protein-coding genes for enrichment
  - The gseapy plotting utils support mapping q values to colors for dot plots but not barplots. They have barplots with fixed float (p or q values) as one axis, and at most color by database groups.
    - thus, here's a simple barplot function for flexibly choosing field as axis and mapping colors
  - Support both explorary saving out all sig terms directly from fresh analysis, or pre-saved tables
    - for pre-saved tables, the plot saving is boosted with plutils package 
- Use magma to enrich for GWAS hits
- Use networkx for Network analysis given any edge list
  - In-degree, Out-degree, and Betweenness centrality
  - Find communities

## Example Usage

```python
import pandas as pd
import numpy as np


data = {
    'Gene_set': ['GO_Biological_Process_2023'] * 5,
    'Term': [
        'Regulation Of Canonical Wnt Signaling Pathway (GO:0060828)',
        'Vascular Transport (GO:0010232)',
        'Positive Regulation Of Telomere Maintenance Via Telomere Lengthening (GO:1904358)'
    ],
    'Overlap': ['6/207', '4/83', '3/35'],
    'P-value': [
        0.0002375382634268994,
        0.00041336149822381393,
        0.0004244897820241493
    ],
    'Adjusted P-value': [
        0.0386193554324011,
        0.0386193554324011,
        0.0386193554324011
    ],
    'Old P-value': [0, 0, 0],
    'Old Adjusted P-value': [0, 0, 0],
    'Odds Ratio': [
        7.544967470340604,
        12.55506329113924,
        23.01388888888889
    ],
    'Combined Score': [
        62.96412550405352,
        97.81885907452366,
        178.6941623125691
    ],
    'Genes': [
        'KANK1;ZNRF3;WNK2;USP34;LGR5;RBMS3',
        'SLC13A3;ABCB1;SLC24A3;LRP2',
        'PKIB;HNRNPA2B1;MAP3K4'
    ]
}

df = pd.DataFrame(data)

genetk.enrich.go.barh(df, 'GO', 'down', 'exact', save = 'test.pdf', return_ax=False)
```
