import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
from collections import defaultdict
from gseapy import barplot


gene_sets=['MSigDB_Hallmark_2020','KEGG_2021_Human', 'GO_Biological_Process_2025', 'Reactome_Pathways_2024', 'Jensen_COMPARTMENTS', 'Human_Phenotype_Ontology', 'DGIdb_Drug_Targets_2024']
## from the intersect output get the shared across all
## get DEG list from file

def _plot(genelist, comparison, pattern):
    print(f"enrichment start {comparison} {pattern}")
    enr_list = []
    ## separately enrich and rbind
    for gset in gene_sets:
        try:
            enr = gp.enrichr(gene_list=genelist, # or "./tests/data/gene_list.txt",
                gene_sets=gset,
                organism='human', # don't forget to set organism to the one you desired! e.g. Yeast
                outdir='output/DEGs/enrichr', # don't write to disk
            )
            enr_list.append(enr.results)
            print(f"enrichment suscess {comparison} up {gset}")
        except:
            enr_list = enr_list
    if len(enr_list)>0:
        enr_pd = pd.concat(enr_list, ignore_index=True)
        barplot(enr_pd,
              column="Adjusted P-value",
              group='Gene_set', # set group, so you could do a multi-sample/library comparsion
              size=10,
              top_term=5,
              figsize=(3,5),
              ofname = f'figures/enrichr/{comparison}_{pattern}_snRNAseq.pdf',
              color = {'MSigDB_Hallmark_2020':'#4C72B0', 'KEGG_2021_Human': '#DD8452',
             'GO_Biological_Process_2023': '#55A868', 'Reactome_Pathways_2024': '#C44E52', 'Jensen_COMPARTMENTS': '#8172B2', 'Human_Phenotype_Ontology': '#9370DB', 'DGIdb_Drug_Targets_2024': '#8A3200'})
        enr_pd.to_csv(f'output/DEGs/enrichr/{comparison}_{pattern}_snRNAseq.tsv', sep='\t')

def enrich_plots(comparison, file):
    # data = pd.read_csv(file, sep=';', decimal=',')
    data = pd.read_csv(file, sep=',', decimal='.')
    cols = ["padj", "log2FoldChange"]
    data[cols] = data[cols].apply(pd.to_numeric, errors='coerce')
    genelists = defaultdict(list)
    data_sig = data.query("padj<0.05")
    print(data.head())
    genelists['up'],  genelists['down'] = data_sig.query("log2FoldChange>1")["Unnamed: 0"].to_list(), data_sig.query("log2FoldChange<(-1)")["Unnamed: 0"].to_list()
    for k,v in genelists.items():
        print(len(v))
        _plot(v, comparison, k)
    

if __name__ == "__main__":
    comparison = "placebo_control"
    file = f"output/DEGs/up_{comparison}_pc_snRNAseq.csv"
    enrich_plots(comparison, file)
    file = f"output/DEGs/down_{comparison}_pc_snRNAseq.csv"
    enrich_plots(comparison, file)