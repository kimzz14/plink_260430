prefix=$1
threadN=24

plink2 \
    --pmerge-list pfile.list \
    --threads ${threadN} \
    --sample-inner-join \
    --allow-extra-chr \
    --make-pgen \
    --sort-vars \
    --out result/${prefix}.pfile \
    1>    result/${prefix}.pfile.log2 \
    2>    result/${prefix}.pfile.err
