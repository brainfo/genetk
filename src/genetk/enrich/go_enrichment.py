"""
This script performs enrichment analysis and GSEA on differential expression data.

Especially for explorations, to save out and plot all significant terms.

"""
import argparse
import os
import warnings
from datetime import datetime

import gseapy as gp
import pandas as pd
from gseapy import barplot

def gene_info(x):
    """Extract gene name and type from gtf attribute string"""
    g_name = list(filter(lambda x: 'gene_name' in x, x.split(";")))[0].split(" ")[2].strip('\"')
    g_type = list(filter(lambda x: 'gene_type' in x, x.split(";")))[0].split(" ")[2]
    return (g_name, g_type)

def load_gencode_annotation(gtf_file):
    """Load and process gencode annotation"""
    if gtf_file is None:
        return None, None
        
    gencode = pd.read_table(gtf_file, comment="#",
                        sep="\t", names=['seqname', 'source', 'feature', 'start', 'end', 
                                        'score', 'strand', 'frame', 'attribute'])
    
    gencode_genes = gencode[(gencode.feature == "transcript")][['attribute']].copy()
    gencode_genes["gene_name"], gencode_genes["gene_type"] = zip(
        *gencode_genes.attribute.apply(lambda x: gene_info(x)))
    pc_genes = gencode_genes.query("gene_type=='\"protein_coding\"'")
    pc_gene_set = set(pc_genes.gene_name)
    
    # Debug: Print some information about the protein coding genes
    print(f"\nDebug - Protein coding genes info:")
    print(f"Total number of protein coding genes: {len(pc_gene_set)}")
    print("First 5 protein coding genes:", list(pc_gene_set)[:5])
    return gencode_genes, pc_gene_set

def enrich_plots(
    genelist,
    name,
    gene_sets=['MSigDB_Hallmark_2020', 'KEGG_2021_Human', 'GO_Biological_Process_2023'],
    coding=False,
    pc_gene_set=None,
    output_dir='data/enrichr',
    figure_dir='figures/enrichr',
):
    """
    Perform enrichment analysis and create plots.
    
    Args:
        genelist (list): List of genes to analyze
        name (str): Base name for output files
        gene_sets (list): List of gene sets to use for enrichment
        coding (bool): If True, filter for protein-coding genes only
        pc_gene_set (set): Set of protein-coding genes if coding=True
        output_dir (str): Base directory where Enrichr outputs are written
        figure_dir (str): Base directory where plots are saved
    """
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(figure_dir, exist_ok=True)
    analysis_outdir = os.path.join(output_dir, name)
    os.makedirs(analysis_outdir, exist_ok=True)

    # Filter for coding genes if requested
    if coding and pc_gene_set is not None:
        original_len = len(genelist)
        print(f"\nDebug - Gene list info for {name}:")
        print(f"Original number of genes: {original_len}")
        print("First 5 genes in list:", genelist[:5])
        
        # Check for any overlap
        overlap = set(genelist).intersection(pc_gene_set)
        print(f"Number of genes that overlap with protein coding set: {len(overlap)}")
        if len(overlap) == 0:
            print("First 5 genes in input list:", genelist[:5])
            print("First 5 genes in protein coding set:", list(pc_gene_set)[:5])
            print("Are the gene names in the same format? (e.g., case sensitivity, version numbers)")
        
        genelist = [gene for gene in genelist if gene in pc_gene_set]
        if len(genelist) == 0:
            warnings.warn(f"No coding genes found after filtering for {name}. Original list had {original_len} genes.")
            return
        name = f"{name}_coding"
        print(f"Filtered to {len(genelist)} coding genes")
    
    enr_list = []
    for gene_st in gene_sets:
        try:
            print(f"enrichment start, genes in {gene_st}")
            enr = gp.enrichr(gene_list=genelist,
                 gene_sets=gene_st,
                 organism='human',
                 outdir=analysis_outdir,
                )
            enr_list.append(enr.results)
            print(f"enrichment success for {gene_st}")
        except:
            enr_list = enr_list
    if len(enr_list) > 0:
        enr_pd = pd.concat(enr_list, ignore_index=True)
        ## filter out any rows with a p-value of 1
        enr_pd = enr_pd[enr_pd['Adjusted P-value'] < 0.05]
        enr_pd.to_csv(os.path.join(output_dir, f'{name}.tsv'), sep='\t')
        top_enr_pd = enr_pd.head(15)
        barplot(top_enr_pd,
              column="Adjusted P-value",
              group='Gene_set',
              cutoff=0.05,
              size=10,
              figsize=(1.77,2.23*len(top_enr_pd.index)/10),
              ofname=os.path.join(figure_dir, f'{name}_.pdf'),
              color={'MSigDB_Hallmark_2020':'#4C72B0', 
                     'KEGG_2021_Human': '#DD8452',
                     'GO_Biological_Process_2023': '#55A868'})

def rnk_gsea(
    rnk,
    name,
    gene_set='KEGG_2021_Human',
    coding=False,
    pc_gene_set=None,
    output_dir='data/enrichr',
    figure_dir='figures/enrichr',
):
    """
    Perform GSEA analysis and create plots.
    
    Args:
        rnk (pd.DataFrame): DataFrame with gene names and log fold changes
        name (str): Base name for output files
        gene_set (str): Gene set to use for GSEA analysis
        coding (bool): If True, filter for protein-coding genes only
        pc_gene_set (set): Set of protein-coding genes if coding=True
        output_dir (str): Base directory where GSEA outputs are written
        figure_dir (str): Base directory where GSEA plots are saved
    """
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(figure_dir, exist_ok=True)
    gene_set_outdir = os.path.join(output_dir, gene_set)
    os.makedirs(gene_set_outdir, exist_ok=True)
    gsea_outdir = os.path.join(gene_set_outdir, name)
    os.makedirs(gsea_outdir, exist_ok=True)
    gene_set_fig_dir = os.path.join(figure_dir, gene_set)
    os.makedirs(gene_set_fig_dir, exist_ok=True)
    gsea_fig_path = os.path.join(gene_set_fig_dir, f'{name}__gsea.pdf')

    # Filter for coding genes if requested
    if coding and pc_gene_set is not None:
        original_len = len(rnk)
        print(f"\nDebug - GSEA gene list info for {name}:")
        print(f"Original number of genes: {original_len}")
        print("First 5 genes in list:", rnk.index[:5].tolist())
        
        # Check for any overlap
        overlap = set(rnk.index).intersection(pc_gene_set)
        print(f"Number of genes that overlap with protein coding set: {len(overlap)}")
        if len(overlap) == 0:
            print("First 5 genes in input list:", rnk.index[:5].tolist())
            print("First 5 genes in protein coding set:", list(pc_gene_set)[:5])
            print("Are the gene names in the same format? (e.g., case sensitivity, version numbers)")
        
        rnk = rnk[rnk.index.isin(pc_gene_set)]
        if len(rnk) == 0:
            warnings.warn(f"No coding genes found after filtering for {name} GSEA. Original list had {original_len} genes.")
            return
        name = f"{name}_coding"
        print(f"Filtered to {len(rnk)} coding genes")
    
    pre_res = gp.prerank(rnk=rnk,
                     gene_sets=gene_set,
                     threads=16,
                     min_size=5,
                     max_size=1000,
                     permutation_num=1000,
                     outdir=gsea_outdir,
                     seed=404,
                     verbose=True)

    terms = pre_res.res2d[pre_res.res2d['FDR q-val'] < 0.05].Term
    pre_res.plot(terms=terms,
                   show_ranking=True,
                   figsize=(1.77,2.23),
                   ofname=gsea_fig_path)

def create_gene_lists(data):
    """
    Create two gene lists based on differential expression criteria.
    
    Args:
        data (pd.DataFrame): DataFrame with differential expression results
        
    Returns:
        tuple: (upregulated_genes, downregulated_genes)
    """
    # Upregulated genes
    up_genes = data[
        (data['avg_log2FC'] > 1) & 
        (data['p_val_adj'] < 0.05) & 
        (abs(data['pct.1'] - data['pct.2']) > 0.1)
    ].index.tolist()
    
    # Downregulated genes
    down_genes = data[
        (data['avg_log2FC'] < (-1)) & 
        (data['p_val_adj'] < 0.05) & 
        (abs(data['pct.1'] - data['pct.2']) > 0.1)
    ].index.tolist()
    
    return up_genes, down_genes

def process_excel_sheets(
    excel_file,
    gtf_file=None,
    gene_sets=['MSigDB_Hallmark_2020', 'KEGG_2021_Human', 'GO_Biological_Process_2023'],
    coding=True,
    run_id=None,
):
    """
    Process each sheet in the Excel file and run enrichment analysis.
    
    Args:
        excel_file (str): Path to the Excel file
        gtf_file (str): Path to the GTF file for gene annotation
        gene_sets (list): List of gene sets to use for analysis
        coding (bool): Whether to filter for protein-coding genes
        run_id (str): Unique identifier to segregate outputs per run. If None,
            defaults to a timestamp.
    """
    if run_id is None:
        run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    safe_run_id = str(run_id)
    for sep in (os.sep, os.path.altsep):
        if sep:
            safe_run_id = safe_run_id.replace(sep, "_")

    run_output_dir = os.path.join('data', 'enrichr', safe_run_id)
    run_figure_dir = os.path.join('figures', 'enrichr', safe_run_id)
    os.makedirs(run_output_dir, exist_ok=True)
    os.makedirs(run_figure_dir, exist_ok=True)
    print(f"Run outputs will be stored under: {run_output_dir}")
    print(f"Run figures will be stored under: {run_figure_dir}")
    
    # Load gene annotation if provided
    _, pc_gene_set = load_gencode_annotation(gtf_file)
    
    # Read all sheets from Excel file
    excel_data = pd.read_excel(excel_file, sheet_name=None, index_col=0)
    
    for sheet_name, data in excel_data.items():
        print(f"\nProcessing sheet: {sheet_name}")
        print("First few rows of data:")
        print(data.head())
        print("\nColumn names:", data.columns.tolist())
        
        # Create gene lists
        up_genes, down_genes = create_gene_lists(data)
        
        # Run enrichment analysis for upregulated genes
        if up_genes:
            print(f"Running enrichment for {len(up_genes)} upregulated genes")
            enrich_plots(
                up_genes, 
                f"{sheet_name}_upregulated",
                gene_sets=gene_sets,
                coding=coding,
                pc_gene_set=pc_gene_set,
                output_dir=run_output_dir,
                figure_dir=run_figure_dir,
            )
        
        # Run enrichment analysis for downregulated genes
        if down_genes:
            print(f"Running enrichment for {len(down_genes)} downregulated genes")
            enrich_plots(
                down_genes,
                f"{sheet_name}_downregulated",
                gene_sets=gene_sets,
                coding=coding,
                pc_gene_set=pc_gene_set,
                output_dir=run_output_dir,
                figure_dir=run_figure_dir,
            )
        
        # Run GSEA analysis
        # Create rank file from log2FC
        rnk = data[['avg_log2FC']].copy()
        rnk = rnk.sort_values(by='avg_log2FC', ascending=False)
        
        # Run GSEA for each gene set
        for gene_set in gene_sets:
            print(f"Running GSEA for {gene_set}")
            rnk_gsea(
                rnk,
                sheet_name,
                gene_set=gene_set,
                coding=coding,
                pc_gene_set=pc_gene_set,
                output_dir=run_output_dir,
                figure_dir=run_figure_dir,
            )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run GO/KEGG enrichment and GSEA across Excel sheets."
    )
    parser.add_argument(
        "--excel-file",
        default="test/test.xlsx",
        help="Path to the Excel file containing differential expression results.",
    )
    parser.add_argument(
        "--gtf-file",
        default='test/genes.gtf',
        help="Optional path to a GTF file for retrieving protein-coding annotations.",
    )
    parser.add_argument(
        "--run-id",
        default=None,
        help="Unique identifier appended to output folders. Defaults to a timestamp.",
    )
    parser.add_argument(
        "--no-coding",
        action="store_true",
        help="Set to disable filtering for protein-coding genes.",
    )
    args = parser.parse_args()

    process_excel_sheets(
        args.excel_file,
        gtf_file=args.gtf_file,
        coding=not args.no_coding,
        run_id=args.run_id,
    )
