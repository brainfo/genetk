## make geneset gmt file ##
## input file path ##
## input file format with field index and genes ##
## output field "geneset" \t "description" \t "genes" ##
## geneset is the index column from input, description is the index column value + input file names without suffix .tsv, genes is the genes column, aggregated by the index column, separated by space into one record ##

import pandas as pd
import os
import mygene

def gene_symbol_to_entrez_id(gene_symbols):
    mg = mygene.MyGeneInfo()
    gene_symbols_list = gene_symbols.tolist()
    results = mg.querymany(gene_symbols_list, scopes='symbol', fields='entrezgene', species='human')
    df_results = pd.DataFrame(results)
    return df_results

def mk_gmt(df, descript):
    ## drop na any, i.e., the rows that contains any nan
    df = df.dropna(subset=["index", "genes"])
    ## df aggreagate by "index" and genes to a str joined by space
    ## convert genes to entrez id
    entrez_df = gene_symbol_to_entrez_id(df["genes"])
    ## merge with df
    df = pd.merge(df, entrez_df, left_on="genes", right_on="query", how="left")
    ## drop the original genes column
    df = df.drop(columns=["genes"])
    ## rename the entrezgene column to genes
    df = df.rename(columns={"entrezgene": "genes"})
    df = df.dropna(subset=["index", "genes"])
    gmt_df = df.groupby("index")["genes"].apply(lambda x: " ".join(x)).reset_index()
    ## add description
    gmt_df["description"] = gmt_df["index"] + "_" + descript
    ## add geneset
    gmt_df["geneset"] = gmt_df["index"] + "_" + descript
    ## select columns
    gmt_df = gmt_df[["geneset", "description", "genes"]]
    return gmt_df

if __name__ == "__main__":
    ## read all files in the input path
    files = os.listdir("/mnt/data/hong/2022/DHJ1_human_obesity_placenta/output/DEGs/summary_Cscore")
    gmt_dfs = []
    ## read each file
    for file in files:
        df = pd.read_csv(f"/mnt/data/hong/2022/DHJ1_human_obesity_placenta/output/DEGs/summary_Cscore/{file}", sep="\t")
        descript = file.split(".")[0]
        ## make gmt
        gmt_dfs.append(mk_gmt(df, descript))
    ## pandas concat
    gmt = pd.concat(gmt_dfs)
    ## output
    gmt.to_csv("genese_obesity_placenta.gmt", sep="\t", index=False)