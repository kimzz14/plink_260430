threadN=24

plink2 \
    --pfile result/pooled.HaplotypeCaller.ChrAll.all.pfile \
    --threads ${threadN} \
    --pheno pheno.tsv \
    --pheno-name HD.STD \
    --covar result/pooled.HaplotypeCaller.ChrAll.all.prune.pca.eigenvec \
    --covar-name PC1,PC2 \
    --glm recessive hide-covar cols=+a1freq,+nobs,+beta,+se,+p \
    --out gwas4/gwas \
    1>    gwas4/gwas.log \
    2>    gwas4/gwas.err
