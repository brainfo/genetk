import pandas as pd
import numpy as np

# Create example gene set enrichment analysis results
data = {
    'Gene_set': ['GO_Biological_Process_2023'] * 5,
    'Term': [
        'Regulation Of Canonical Wnt Signaling Pathway (GO:0060828)',
        'Vascular Transport (GO:0010232)',
        'Positive Regulation Of Telomere Maintenance Via Telomere Lengthening (GO:1904358)',
        'Regulation Of p38MAPK Cascade (GO:1900744)',
        'Regulation Of ERK1 And ERK2 Cascade (GO:0070372)'
    ],
    'Overlap': ['6/207', '4/83', '3/35', '3/37', '6/244'],
    'P-value': [
        0.0002375382634268994,
        0.00041336149822381393,
        0.0004244897820241493,
        0.0005008995516524148,
        0.0005697290270008587
    ],
    'Adjusted P-value': [
        0.0386193554324011,
        0.0386193554324011,
        0.0386193554324011,
        0.0386193554324011,
        0.0399328254379692
    ],
    'Old P-value': [0, 0, 0, 0, 0],
    'Old Adjusted P-value': [0, 0, 0, 0, 0],
    'Odds Ratio': [
        7.544967470340604,
        12.55506329113924,
        23.01388888888889,
        21.65795206971677,
        6.360051712992889
    ],
    'Combined Score': [
        62.96412550405352,
        97.81885907452366,
        178.6941623125691,
        164.5810512711608,
        47.51181041502257
    ],
    'Genes': [
        'KANK1;ZNRF3;WNK2;USP34;LGR5;RBMS3',
        'SLC13A3;ABCB1;SLC24A3;LRP2',
        'PKIB;HNRNPA2B1;MAP3K4',
        'DLG1;PHLPP1;MAP3K4',
        'SPRED2;DLG1;ERBB4;WNK2;PDGFD;ARRB1'
    ]
}

# Create DataFrame
gsea_results = pd.DataFrame(data)

# Example usage with genetk
if __name__ == '__main__':
    print("Example GSEA results DataFrame:")
    print(gsea_results)
    
    # Example of how to use with genetk plotting functions
    # from genetk import plot_gsea_results
    # plot_gsea_results(gsea_results, 
    #                  value_col='Adjusted P-value',
    #                  color_col='Combined Score',
    #                  top_n=5) 