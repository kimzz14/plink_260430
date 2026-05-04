import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


gwas_file = "gwas.HD.STD.glm.linear"
out_prefix = "trait_gwas"

CHUNKSIZE = 500_000


# -----------------------------
# 0. Header 확인
# -----------------------------
header = pd.read_csv(gwas_file, sep=r"\s+", nrows=0).columns.tolist()

chrom_col = "#CHROM" if "#CHROM" in header else "CHROM"

required_cols = [chrom_col, "POS", "P"]
missing_cols = [c for c in required_cols if c not in header]

if missing_cols:
    raise ValueError(f"Missing required columns: {missing_cols}")

usecols = [chrom_col, "POS", "P"]

print("Columns loaded:", usecols)


# -----------------------------
# 1. 첫 번째 pass:
#    유효 row 수, chromosome min/max 계산
# -----------------------------
chrom_min = {}
chrom_max = {}
n_valid = 0

for chunk in pd.read_csv(
    gwas_file,
    sep=r"\s+",
    usecols=usecols,
    chunksize=CHUNKSIZE
):
    if chrom_col != "CHROM":
        chunk = chunk.rename(columns={chrom_col: "CHROM"})

    chunk = chunk.dropna(subset=["P", "POS"])
    chunk = chunk[chunk["P"] > 0]

    if chunk.empty:
        continue

    chunk["CHROM"] = chunk["CHROM"].astype(str)

    n_valid += len(chunk)

    grouped = chunk.groupby("CHROM")["POS"].agg(["min", "max"])

    for chrom, row in grouped.iterrows():
        if chrom not in chrom_min:
            chrom_min[chrom] = row["min"]
            chrom_max[chrom] = row["max"]
        else:
            chrom_min[chrom] = min(chrom_min[chrom], row["min"])
            chrom_max[chrom] = max(chrom_max[chrom], row["max"])

if n_valid == 0:
    raise ValueError("No valid GWAS rows found after filtering.")

print(f"Valid rows: {n_valid:,}")


# -----------------------------
# 2. chromosome 순서 설정
# -----------------------------
wheat_chr_order = [
    "1A", "1B", "1D",
    "2A", "2B", "2D",
    "3A", "3B", "3D",
    "4A", "4B", "4D",
    "5A", "5B", "5D",
    "6A", "6B", "6D",
    "7A", "7B", "7D"
]

observed_chroms = set(chrom_min.keys())

chr_order = [c for c in wheat_chr_order if c in observed_chroms]

if len(chr_order) == 0:
    chr_order = sorted(
        observed_chroms,
        key=lambda x: int(x) if str(x).isdigit() else str(x)
    )

print("Chromosomes:", chr_order)


# -----------------------------
# 3. Manhattan plot용 offset 계산
# -----------------------------
chrom_offsets = {}
tick_positions = []
tick_labels = []

current_offset = 0

for chrom in chr_order:
    chrom_offsets[chrom] = current_offset

    cmin = chrom_min[chrom]
    cmax = chrom_max[chrom]

    tick_positions.append(current_offset + (cmin + cmax) / 2)
    tick_labels.append(chrom)

    current_offset += cmax


# -----------------------------
# 4. 염색체별 2색 반복 지정
# -----------------------------
alt_colors = ["#4C72B0", "#DD8452"]   # 파랑, 주황
chrom_colors = {
    chrom: alt_colors[i % 2]
    for i, chrom in enumerate(chr_order)
}


# -----------------------------
# 5. 두 번째 pass:
#    Manhattan plot을 chunk 단위로 그림
#    QQ plot용 P-value만 저장
# -----------------------------
pvals = np.empty(n_valid, dtype=np.float64)
p_index = 0

plt.figure(figsize=(14, 6))

for chunk in pd.read_csv(
    gwas_file,
    sep=r"\s+",
    usecols=usecols,
    chunksize=CHUNKSIZE
):
    if chrom_col != "CHROM":
        chunk = chunk.rename(columns={chrom_col: "CHROM"})

    chunk = chunk.dropna(subset=["P", "POS"])
    chunk = chunk[chunk["P"] > 0]

    if chunk.empty:
        continue

    chunk["CHROM"] = chunk["CHROM"].astype(str)
    chunk["minus_log10_p"] = -np.log10(chunk["P"].to_numpy())

    # QQ plot용 P-value 저장
    this_n = len(chunk)
    pvals[p_index:p_index + this_n] = chunk["P"].to_numpy()
    p_index += this_n

    # Manhattan plot: chromosome별로 chunk 안에서만 그림
    for chrom in chr_order:
        sub = chunk[chunk["CHROM"] == chrom]

        if sub.empty:
            continue

        x = sub["POS"].to_numpy() + chrom_offsets[chrom]
        y = sub["minus_log10_p"].to_numpy()

        plt.scatter(
            x,
            y,
            s=4,
            alpha=0.8,
            color=chrom_colors[chrom]
        )

# 혹시 남는 공간이 있으면 자름
pvals = pvals[:p_index]

# threshold
bonf = 0.05 / len(pvals)
suggestive = 1 / len(pvals)

plt.axhline(-np.log10(bonf), linestyle="--", linewidth=1)
plt.axhline(-np.log10(suggestive), linestyle=":", linewidth=1)

plt.xticks(tick_positions, tick_labels, rotation=90)
plt.xlabel("Chromosome")
plt.ylabel("-log10(P)")
plt.title("Manhattan plot")

plt.tight_layout()
plt.savefig(f"{out_prefix}_manhattan.png", dpi=300, bbox_inches="tight")
plt.close()


# -----------------------------
# 6. QQ plot
# -----------------------------
observed = -np.log10(np.sort(pvals))
expected = -np.log10(
    np.arange(1, len(observed) + 1) / (len(observed) + 1)
)

plt.figure(figsize=(6, 6))
plt.scatter(expected, observed, s=5, alpha=0.7)

max_val = max(expected.max(), observed.max())
plt.plot([0, max_val], [0, max_val], linestyle="--", linewidth=1)

plt.xlabel("Expected -log10(P)")
plt.ylabel("Observed -log10(P)")
plt.title("QQ plot")

plt.tight_layout()
plt.savefig(f"{out_prefix}_qq.png", dpi=300, bbox_inches="tight")
plt.close()


print("Saved:")
print(f"{out_prefix}_manhattan.png")
print(f"{out_prefix}_qq.png")