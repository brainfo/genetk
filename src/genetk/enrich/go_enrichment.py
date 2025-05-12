import pandas as pd
import gseapy as gp
from gseapy import barplot
import os
import plutils as plu
import warnings

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

def enrich_plots(genelist, name, gene_sets, coding=False, pc_gene_set=None):
    """
    Perform enrichment analysis and create plots.
    
    Args:
        genelist (list): List of genes to analyze
        name (str): Base name for output files
        gene_sets (list): List of gene sets to use for enrichment
        coding (bool): If True, filter for protein-coding genes only
        pc_gene_set (set): Set of protein-coding genes if coding=True
    """
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
                 outdir='data/enrichr',
                )
            enr_list.append(enr.results)
            print(f"enrichment success for {gene_st}")
        except:
            enr_list = enr_list
    if len(enr_list) > 0:
        enr_pd = pd.concat(enr_list, ignore_index=True)
        barplot(enr_pd,
              column="Adjusted P-value",
              group='Gene_set',
              size=10,
              top_term=5,
              figsize=(1.77,2.23),
              ofname=f'figures/enrichr/{name}_pcos_specific.pdf',
              color={'MSigDB_Hallmark_2020':'#4C72B0', 
                     'KEGG_2021_Human': '#DD8452',
                     'GO_Biological_Process_2023': '#55A868'})
        enr_pd.to_csv(f'data/enrichr/{name}.tsv', sep='\t')

def rnk_gsea(rnk, name, gene_set='KEGG_2021_Human', coding=False, pc_gene_set=None):
    """
    Perform GSEA analysis and create plots.
    
    Args:
        rnk (pd.DataFrame): DataFrame with gene names and log fold changes
        name (str): Base name for output files
        gene_set (str): Gene set to use for GSEA analysis
        coding (bool): If True, filter for protein-coding genes only
        pc_gene_set (set): Set of protein-coding genes if coding=True
    """
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
    
    # Create output directory for this gene set
    outdir = f'data/enrichr/{gene_set}'
    os.makedirs(outdir, exist_ok=True)
    
    pre_res = gp.prerank(rnk=rnk,
                     gene_sets=gene_set,
                     threads=16,
                     min_size=5,
                     max_size=1000,
                     permutation_num=1000,
                     outdir=outdir,
                     seed=404,
                     verbose=True)
    
    terms = pre_res.res2d.Term
    pre_res.plot(terms=terms[1:7],
                   show_ranking=True,
                   figsize=(1.77,2.23),
                   ofname=f'figures/enrichr/{gene_set}/{name}_pcos_specific_gsea.pdf')

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
        (data['avg_log2FC'] > 0.25) & 
        (data['p_val_adj'] < 0.05) & 
        (abs(data['pct.1'] - data['pct.2']) > 0.1)
    ].index.tolist()
    
    # Downregulated genes
    down_genes = data[
        (data['avg_log2FC'] < -0.25) & 
        (data['p_val_adj'] < 0.05) & 
        (abs(data['pct.1'] - data['pct.2']) > 0.1)
    ].index.tolist()
    
    return up_genes, down_genes

def process_excel_sheets(excel_file, gtf_file=None, gene_sets=None, coding=True):
    """
    Process each sheet in the Excel file and run enrichment analysis.
    
    Args:
        excel_file (str): Path to the Excel file
        gtf_file (str): Path to the GTF file for gene annotation
        gene_sets (list): List of gene sets to use for analysis
        coding (bool): Whether to filter for protein-coding genes
    """
    if gene_sets is None:
        gene_sets = ['MSigDB_Hallmark_2020', 'KEGG_2021_Human', 'GO_Biological_Process_2023']
    
    # Create output directories
    os.makedirs('data/enrichr', exist_ok=True)
    os.makedirs('figures/enrichr', exist_ok=True)
    
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
                pc_gene_set=pc_gene_set
            )
        
        # Run enrichment analysis for downregulated genes
        if down_genes:
            print(f"Running enrichment for {len(down_genes)} downregulated genes")
            enrich_plots(
                down_genes,
                f"{sheet_name}_downregulated",
                gene_sets=gene_sets,
                coding=coding,
                pc_gene_set=pc_gene_set
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
                pc_gene_set=pc_gene_set
            )

if __name__ == "__main__":
    wd = '/mnt/run/jh/projects/pcos_sc'
    os.chdir(wd)
    
    # Example usage
    excel_file = "data/degs/data/pcos_sex.xlsx"
    gtf_file = '/mnt/run/jh/reference/Homo_sapiens/NCBI/GRCh38/Annotation/Genes.gencode/genes.gtf'
    
    plu.reset_mpl_style(rcfile="/mnt/run/jh/reference/general.mplstyle")
    process_excel_sheets(excel_file, gtf_file=gtf_file)