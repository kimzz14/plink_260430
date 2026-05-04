prefix=$1
threadN=24

plink2 \
    --pfile result/${prefix}.pfile \
    --threads ${threadN} \
    --allow-extra-chr \
    --blocks no-pheno-req \
    --out result/${prefix}.ld_blokcs \
    1>    result/${prefix}.ld_blokcs.log \
    2>    result/${prefix}.ld_blokcs.err
