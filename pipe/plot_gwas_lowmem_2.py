import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


gwas_file = "gwas.HD.STD.glm.linear"
out_prefix = "trait_gwas"
genes_file = "known_genes.tsv"   # 유전자 annotation 파일. 없으면 None으로 설정 가능

CHUNKSIZE = 1_000_000
SEP = r"\s+"


# -----------------------------
# 0. Header 확인
# -----------------------------
header = pd.read_csv(gwas_file, sep=SEP, nrows=0).columns.tolist()
chrom_col = "#CHROM" if "#CHROM" in header else "CHROM"

required_cols = [chrom_col, "POS", "P"]
missing_cols = [c for c in required_cols if c not in header]
if missing_cols:
    raise ValueError(f"Missing required columns: {missing_cols}")

usecols = [chrom_col, "POS", "P"]
print("Columns loaded:", usecols)


# -----------------------------
# 1. 첫 번째 pass:
#    chromosome별 min/max와 valid row 수 계산
# -----------------------------
chrom_min = {}
chrom_max = {}
n_valid = 0

for chunk in pd.read_csv(
    gwas_file,
    sep=SEP,
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

    grouped = chunk.groupby("CHROM", sort=False)["POS"].agg(["min", "max"])

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
# 2. chromosome 순서
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
# 3. offset, tick, color 설정
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

alt_colors = ["#4C72B0", "#DD8452"]
chrom_colors = {
    chrom: alt_colors[i % 2]
    for i, chrom in enumerate(chr_order)
}


# -----------------------------
# 4. 두 번째 pass:
#    Manhattan plot + QQ용 p-value 저장
# -----------------------------
pvals = np.empty(n_valid, dtype=np.float64)
p_index = 0
global_ymax = 0.0

fig, ax = plt.subplots(figsize=(14, 6))

for chunk in pd.read_csv(
    gwas_file,
    sep=SEP,
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

    p = chunk["P"].to_numpy(dtype=np.float64)
    minus_log10_p = -np.log10(p)

    this_ymax = np.nanmax(minus_log10_p)
    if this_ymax > global_ymax:
        global_ymax = this_ymax

    this_n = len(p)
    pvals[p_index:p_index + this_n] = p
    p_index += this_n

    for chrom, sub in chunk.groupby("CHROM", sort=False):
        if chrom not in chrom_offsets:
            continue

        idx = sub.index.to_numpy()
        loc = chunk.index.get_indexer(idx)

        x = sub["POS"].to_numpy(dtype=np.int64) + chrom_offsets[chrom]
        y = minus_log10_p[loc]

        ax.scatter(
            x,
            y,
            s=4,
            alpha=0.8,
            color=chrom_colors[chrom],
            rasterized=True
        )

pvals = pvals[:p_index]


# -----------------------------
# 5. threshold 선
# -----------------------------
bonf = 0.05 / len(pvals)
suggestive = 1 / len(pvals)

ax.axhline(-np.log10(bonf), linestyle="--", linewidth=1)
ax.axhline(-np.log10(suggestive), linestyle=":", linewidth=1)


# -----------------------------
# 6. 알려진 유전자 annotation 추가
# -----------------------------
if genes_file is not None and os.path.exists(genes_file):
    print(genes_file)
    genes = pd.read_csv(genes_file, sep="\t")

    required_gene_cols = ["Gene", "CHROM", "START", "END"]
    missing_gene_cols = [c for c in required_gene_cols if c not in genes.columns]
    if missing_gene_cols:
        raise ValueError(f"Missing required columns in {genes_file}: {missing_gene_cols}")

    genes = genes.copy()
    genes["CHROM"] = genes["CHROM"].astype(str)
    genes["START"] = pd.to_numeric(genes["START"], errors="coerce")
    genes["END"] = pd.to_numeric(genes["END"], errors="coerce")
    genes = genes.dropna(subset=["START", "END"])
    genes["START"] = genes["START"].astype(np.int64)
    genes["END"] = genes["END"].astype(np.int64)

    # 현재 Manhattan plot에 존재하는 chromosome만 유지
    genes = genes[genes["CHROM"].isin(chr_order)].copy()

    if not genes.empty:
        # 유전자 표시를 위해 y축 위쪽 여유 확보
        gene_track_y = global_ymax * 1.05
        gene_label_y = global_ymax * 1.09
        ax.set_ylim(0, global_ymax * 1.18)

        for _, row in genes.iterrows():
            chrom = row["CHROM"]
            start = int(row["START"])
            end = int(row["END"])
            gene = str(row["Gene"])

            if chrom not in chrom_offsets:
                continue

            x_start = start + chrom_offsets[chrom]
            x_end = end + chrom_offsets[chrom]
            x_mid = (x_start + x_end) / 2

            # 유전자 구간선
            ax.hlines(
                y=gene_track_y,
                xmin=x_start,
                xmax=x_end,
                color="red",
                linewidth=2
            )

            # 중앙 표시선
            ax.vlines(
                x=x_mid,
                ymin=gene_track_y - global_ymax * 0.015,
                ymax=gene_track_y + global_ymax * 0.015,
                color="red",
                linewidth=1
            )

            # 유전자 이름
            ax.text(
                x_mid,
                gene_label_y,
                gene,
                rotation=90,
                ha="center",
                va="bottom",
                fontsize=8,
                color="red"
            )


# -----------------------------
# 7. Manhattan plot 저장
# -----------------------------
ax.set_xticks(tick_positions)
ax.set_xticklabels(tick_labels, rotation=90)
ax.set_xlabel("Chromosome")
ax.set_ylabel("-log10(P)")
ax.set_title("Manhattan plot")

fig.tight_layout()
fig.savefig(f"{out_prefix}_manhattan.png", dpi=300, bbox_inches="tight")
plt.close(fig)


# -----------------------------
# 8. QQ plot
# -----------------------------
observed = -np.log10(np.sort(pvals))
expected = -np.log10(
    np.arange(1, len(observed) + 1) / (len(observed) + 1)
)

fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(
    expected,
    observed,
    s=5,
    alpha=0.7,
    rasterized=True
)

max_val = max(expected.max(), observed.max())
ax.plot([0, max_val], [0, max_val], linestyle="--", linewidth=1)

ax.set_xlabel("Expected -log10(P)")
ax.set_ylabel("Observed -log10(P)")
ax.set_title("QQ plot")

fig.tight_layout()
fig.savefig(f"{out_prefix}_qq.png", dpi=300, bbox_inches="tight")
plt.close(fig)

print("Saved:")
print(f"{out_prefix}_manhattan.png")
print(f"{out_prefix}_qq.png")