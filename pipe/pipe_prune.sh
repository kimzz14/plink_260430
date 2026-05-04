prefix=$1
threadN=24

plink2 \
    --pfile result/${prefix}.pfile \
    --threads ${threadN} \
    --allow-extra-chr \
    --indep-pairwise 1000 50 0.1 \
    --out result/${prefix} \
    1>    result/${prefix}.prune.log \
    2>    result/${prefix}.prune.err

plink2 \
    --pfile result/${prefix}.pfile \
    --allow-extra-chr \
    --extract result/${prefix}.prune.in \
    --make-pgen \
    --out result/${prefix}.prune.pfile \
    1>    result/${prefix}.prune.pfile.log2 \
    2>    result/${prefix}.prune.pfile.err
