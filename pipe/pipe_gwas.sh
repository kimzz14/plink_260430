threadN=24

plink2 \
    --pfile result/pooled.HaplotypeCaller.ChrAll.all.prune.pfile \
    --threads ${threadN} \
    --pheno pheno.tsv \
    --pheno-name HD.STD \
    --covar result/pooled.HaplotypeCaller.ChrAll.all.prune.pca.eigenvec \
    --covar-name PC1,PC2 \
    --glm hide-covar cols=+a1freq,+nobs,+beta,+se,+p \
    --out gwas/gwas \
    1>    gwas/gwas.log \
    2>    gwas/gwas.err
