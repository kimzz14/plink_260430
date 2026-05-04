import pandas as pd
import matplotlib.pyplot as plt


# Input files
eigenvec_file = "result/pooled.HaplotypeCaller.ChrAll.all.prune.pca.eigenvec"
eigenval_file = "result/pooled.HaplotypeCaller.ChrAll.all.prune.pca.eigenval"
metadata_file = "metadata.tsv"

# Output files
out_png = "PCA_PC1_PC2_by_country.png"
out_pdf = "PCA_PC1_PC2_by_country.pdf"

# 1. Read PCA result
pca = pd.read_csv(eigenvec_file, sep=r"\s+")

# PLINK2 eigenvec may have '#FID' instead of 'FID'
if "#FID" in pca.columns:
    pca = pca.rename(columns={"#FID": "FID"})

# 2. Read eigenvalues
eig = pd.read_csv(eigenval_file, header=None)[0]
var_exp = eig / eig.sum() * 100

# 3. Read metadata
meta = pd.read_csv(metadata_file, sep="\t")
print(meta)

# Required columns check
required_meta_cols = {"IID", "Country"}
missing = required_meta_cols - set(meta.columns)
if missing:
    raise ValueError(f"metadata.tsv is missing required columns: {missing}")


# 4. Merge PCA and metadata
dat = pca.merge(meta, on="IID", how="left")

# Check unmatched samples
n_unmatched = dat["Country"].isna().sum()
if n_unmatched > 0:
    print(f"Warning: {n_unmatched} samples do not have Country metadata.")
    dat["Country"] = dat["Country"].fillna("Unknown")

# 5. Print country sample counts
print("\nSample counts by country:")
print(dat["Country"].value_counts())

# 6. Plot PC1 vs PC2 by country
plt.figure(figsize=(7, 6))

for country, sub in dat.groupby("Country"):
    plt.scatter(
        sub["PC1"],
        sub["PC2"],
        s=35,
        alpha=0.8,
        label=f"{country} (n={len(sub)})"
    )

plt.xlabel(f"PC1 ({var_exp.iloc[0]:.2f}%)")
plt.ylabel(f"PC2 ({var_exp.iloc[1]:.2f}%)")
plt.title("PCA by country")
plt.legend(frameon=False, fontsize=8, loc="best")

plt.tight_layout()
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
plt.close()

print(f"\nSaved: {out_png}")
print(f"Saved: {out_pdf}")