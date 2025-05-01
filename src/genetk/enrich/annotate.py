import pandas as pd

def gene_info(x):
    """Extract gene name and type from gtf attribute string"""
    g_name = list(filter(lambda x: 'gene_name' in x, x.split(";")))[0].split(" ")[2].strip('\"')
    g_type = list(filter(lambda x: 'gene_type' in x, x.split(";")))[0].split(" ")[2]
    return (g_name, g_type)

def get_pc_genes(gtf_file): 
    gencode = pd.read_table(gtf_file, comment="#",
                        sep="\t", names=['seqname', 'source', 'feature', 'start', 'end', 
                                        'score', 'strand', 'frame', 'attribute'])
    
    gencode_genes = gencode[(gencode.feature == "transcript")][['attribute']].copy()
    gencode_genes["gene_name"], gencode_genes["gene_type"] = zip(
        *gencode_genes.attribute.apply(lambda x: gene_info(x)))
    pc_genes = gencode_genes.query("gene_type=='\"protein_coding\"'")
    pc_gene_set = set(pc_genes.gene_name)
    return pc_gene_set

def annotate_genes(genes, pc_gene_set):
    """Annotate genes with gene name and type"""
    return [gene for gene in genes if gene in pc_gene_set]