## loop over all loc files in the folder
# for file in *_loc.txt; do
#     snploc=$file
#     ncbi37=/mnt/data/hong/reference/gwas/magma/NCBI37.3.gene.loc
#     magma --annotate \
#       --snp-loc ${snploc} \
#       --gene-loc ${ncbi37} \
#       --out ${file%.loc} 
# done

## single file
file='EGG-GWAS-BL.txt_loc.txt'
snploc=$file
ncbi37=/mnt/data/hong/reference/gwas/magma/NCBI37.3.gene.loc
magma --annotate \
  --snp-loc ${snploc} \
  --gene-loc ${ncbi37} \
  --out ${file%.loc} 