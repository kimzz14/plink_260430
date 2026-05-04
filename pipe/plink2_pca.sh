prefix=$1
threadN=24

plink2 \
    --pfile result/${prefix}.pfile \
    --threads ${threadN} \
    --allow-extra-chr \
    --pca 10 \
    --out result/${prefix}.pca \
    1>    result/${prefix}.log2 \
    2>    result/${prefix}.pca.err
