prefix=$1
threadN=24

plink2 \
    --vcf vcf/${prefix}.vcf.gz \
    --threads ${threadN} \
    --double-id \
    --allow-extra-chr \
    --set-all-var-ids '@:#' \
    --make-pgen \
    --out result/${prefix}.pfile \
    1>    result/${prefix}.pfile.log2 \
    2>    result/${prefix}.pfile.err

#--keep keep_samples.txt \