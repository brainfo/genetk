import gzip
import os
import pandas as pd
import shutil
import re
import logging

root_folder = '/mnt/data/hong/2022/DHJ1_human_obesity_placenta/data/gwas/egg-consortium.org'
## set logging path, and logging level
logging.basicConfig(filename=f'{root_folder}/mk_db.log', level=logging.DEBUG)
target_folder = '/mnt/storage/hong/2024/egg-consortium'

p = ['p', 'p-value', 'pval', 'p.value', 'p_value']
n=['n', 'nsamples', 'totalsamplesize']
rsid=['hm_rsid', 'rs_id','rsid','snp', 'markername', 'variant_id']
chromosome=['chromosome', 'chrom', 'chromosome_name', 'chr']
position=['position', 'pos', 'bp', 'basepair', 'base_pair_location', 'base_pair_location_start']


def unzip_file(file_path, target_folder):
    """unzip the gz file and move the gz file to a target folder"""
    file_name = file_path.replace('.gz', '')
    with gzip.open(file_path,"rb") as f_in, open(file_name,"wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
        ## move the zip file to a target folder
        shutil.move(file_path, os.path.join(target_folder, file_path.split("/")[-1]))

def mk_pval_file(df, p, n, rsid, filename):
    """
    return a file with fields rs_id, p, and n.
    """
    col_names = df.columns.tolist()
    rsid_cols = [col for col in col_names if col.lower().strip() in rsid]
    p_cols = [col for col in col_names if col.lower().strip() in p]
    n_cols = [col for col in col_names if col.lower().strip() in n]
    ## check if the rsid_cols, p_cols, and n_cols are unique
    if len(rsid_cols) != 1 or len(p_cols) != 1:
        ## not raise error, but print a warning and continue
        logging.warning(f"Warning: There should be only one rsid_col, p_col, and n_col")
        ## -[] what if not unique?
        return
    ## get the column names
    rsid_col = rsid_cols[0]
    p_col = p_cols[0]
    ## if n_cols is empty, use 42212
    ## simplified ifelse statement
    n_col = n_cols[0] if len(n_cols) > 0 else n
    ## create df2use with the column names ['rsid', 'p', 'n']
    df2use = pd.DataFrame({"rsid": df[rsid_col], "p": df[p_col], "n": df[n_col]})
    df2use['p'] = pd.to_numeric(df2use['p'], errors='coerce')
    df2use['n'] = pd.to_numeric(df2use['n'], errors='coerce')
    ## remove any nan row
    df2use = df2use[~df2use.isna().any(axis=1)]
    ## specify the column types using pandas for dataframe but not per column
    df2use = df2use.astype({'rsid': 'str', 'p': 'float', 'n': 'int'})
    ## write the file
    filename_new = f'{filename}_p.txt'
    df2use.to_csv(filename_new, index=False, sep='\t')

## for os walk 
# -[] maybe modulize to handle the breaks of the chain of programs

def mk_loc_file(df, rsid,  chromosome, position, filename):
    """
    return a file with fields rs_id, chromosome,position.
    """
    ## get the column names
    col_names = df.columns.tolist()
    rsid_cols = [col for col in col_names if col.lower().strip() in rsid]
    chromosome_cols = [col for col in col_names if col.lower().strip() in chromosome]
    ## using any(substring in string for substring in substring_list)
    position_cols = [col for col in col_names if col.lower().strip() in position]
    if len(rsid_cols) != 1 or len(chromosome_cols) != 1 or len(position_cols) != 1:
        ## not raise error, but print a warning and continue
        logging.warning(f"Warning: we found {len(rsid_cols)} rsid_col, {len(chromosome_cols)} chromosome_col, and {len(position_cols)} position_col")
        ## -[] what if not unique?
        return
    ## get the column names
    rsid_col = rsid_cols[0]
    chromosome_col = chromosome_cols[0]
    position_col = position_cols[0]
    ## create df2use with the column names ['rsid', 'chromosome', 'position']
    df2use = pd.DataFrame({"rsid": df[rsid_col], "chromosome": df[chromosome_col], "position": df[position_col]})
    ## coerce conversion for each column, pd.to_numeric, then astype
    df2use['chromosome'] = pd.to_numeric(df2use['chromosome'], errors='coerce')
    df2use['position'] = pd.to_numeric(df2use['position'], errors='coerce')
    ## remove any nan row
    df2use = df2use[~df2use.isna().any(axis=1)]
    df2use = df2use.astype({'rsid': 'str', 'chromosome': 'int', 'position': 'int'})
    filename_new = f'{filename}_loc.txt'
    df2use.to_csv(filename_new, index=False, sep='\t')

def run_single_gwas(root, file, target_folder, p, n, rsid, chromosome, position):
    unzip_file(os.path.join(root, file), target_folder)
    file_path = os.path.join(root, file.replace('.gz', ''))
    df = pd.read_csv(file_path, sep='\t')
    filename = file_path.split("/")[-1]
    mk_pval_file(df, p, n, rsid, filename)
    mk_loc_file(df, rsid,  chromosome, position, filename)
    ## remove the unzipped file
    os.remove(os.path.join(root, file.replace('.gz', '')))

if __name__ == "__main__":
    for root, dirs, files in os.walk(root_folder):
        for file in files:
            if file.endswith('.gz'):
                ## logging
                logging.info(f"Processing {os.path.join(root, file)}")
                ## these are defaults
                run_single_gwas(root, file, target_folder, p = p, n=n, rsid=rsid, chromosome=chromosome, position=position)