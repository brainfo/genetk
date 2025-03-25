geneset=genese_obesity_placenta.gmt

for file in *genes.raw; do 
magma \
    --gene-results $file \
    --set-annot ${geneset} \
    --out ${file%.genes.raw}
done