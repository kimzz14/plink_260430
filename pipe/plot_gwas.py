import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


gwas_file = "gwas.HD.STD.glm.linear.Chr3D"
out_prefix = "trait_gwas"

# 1. GWAS 결과 읽기
df = pd.read_csv(gwas_file, sep=r"\s+")

# PLINK2는 #CHROM 또는 CHROM 컬럼명을 쓸 수 있음
if "#CHROM" in df.columns:
    df = df.rename(columns={"#CHROM": "CHROM"})

# 필요한 컬럼만 확인
required_cols = ["CHROM", "POS", "ID", "P"]
missing_cols = [c for c in required_cols if c not in df.columns]

if missing_cols:
    raise ValueError(f"Missing required columns: {missing_cols}")

# 2. P-value 정리
df = df.dropna(subset=["P"])
df = df[df["P"] > 0]
df["minus_log10_p"] = -np.log10(df["P"])

# 3. chromosome 순서 지정
# 밀 chromosome 이름이 1A,1B,1D,...,7D인 경우
wheat_chr_order = [
    "1A", "1B", "1D",
    "2A", "2B", "2D",
    "3A", "3B", "3D",
    "4A", "4B", "4D",
    "5A", "5B", "5D",
    "6A", "6B", "6D",
    "7A", "7B", "7D"
]

# CHROM을 문자열로 처리
df["CHROM"] = df["CHROM"].astype(str)

# 실제 데이터에 있는 chromosome만 사용
chr_order = [c for c in wheat_chr_order if c in set(df["CHROM"])]

# 만약 숫자 chromosome이면 자동 정렬
if len(chr_order) == 0:
    chr_order = sorted(df["CHROM"].unique(), key=lambda x: int(x) if x.isdigit() else x)

df["CHROM"] = pd.Categorical(df["CHROM"], categories=chr_order, ordered=True)
df = df.sort_values(["CHROM", "POS"])

# 4. Manhattan plot용 누적 position 계산
chrom_offsets = {}
current_offset = 0
tick_positions = []
tick_labels = []

for chrom in chr_order:
    chrom_df = df[df["CHROM"] == chrom]
    if chrom_df.empty:
        continue

    chrom_offsets[chrom] = current_offset
    chrom_min = chrom_df["POS"].min()
    chrom_max = chrom_df["POS"].max()

    tick_positions.append(current_offset + (chrom_min + chrom_max) / 2)
    tick_labels.append(chrom)

    current_offset += chrom_max

df["cum_pos"] = df.apply(
    lambda row: row["POS"] + chrom_offsets[str(row["CHROM"])],
    axis=1
)

# 5. Manhattan plot
plt.figure(figsize=(14, 6))

for i, chrom in enumerate(chr_order):
    chrom_df = df[df["CHROM"] == chrom]
    plt.scatter(
        chrom_df["cum_pos"],
        chrom_df["minus_log10_p"],
        s=4,
        alpha=0.8
    )

# Bonferroni threshold
bonf = 0.05 / len(df)
plt.axhline(-np.log10(bonf), linestyle="--", linewidth=1)

# suggestive threshold
suggestive = 1 / len(df)
plt.axhline(-np.log10(suggestive), linestyle=":", linewidth=1)

plt.xticks(tick_positions, tick_labels, rotation=90)
plt.xlabel("Chromosome")
plt.ylabel("-log10(P)")
plt.title("Manhattan plot")
plt.tight_layout()
plt.savefig(f"{out_prefix}_manhattan.png", dpi=300, bbox_inches="tight")
plt.close()

# 6. QQ plot
observed = -np.log10(np.sort(df["P"].values))
expected = -np.log10(np.arange(1, len(observed) + 1) / (len(observed) + 1))

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