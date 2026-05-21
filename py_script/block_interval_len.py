import sys
import pandas as pd

# ===================== 检查参数 =====================
if len(sys.argv) != 3:
    print("用法: python block_interval_len.py <block文件> <样本ID>")
    sys.exit(1)

block_file = sys.argv[1]
sample_id = sys.argv[2]

# ===================== 读取并过滤数据 =====================
df = pd.read_csv(
    block_file,
    sep=r"\s+",
    comment="#",
    names=["sample", "chromosome", "phase_set", "from", "to", "variants"]
)

df = df.sort_values(["chromosome", "from"]).reset_index(drop=True)

# ===================== 计算间隔 =====================
print("sampleID\tgap_length")  # 表头

for i in range(1, len(df)):
    prev_end = df.iloc[i-1]["to"]
    curr_start = df.iloc[i]["from"]
    gap = curr_start - prev_end
    print(f"{sample_id}\t{gap}")
