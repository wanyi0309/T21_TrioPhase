import pandas as pd
import gzip
import argparse
from collections import defaultdict

def main():
    # ===================== 命令行参数解析 =====================
    parser = argparse.ArgumentParser(description="Block and GT comparison script for phased VCF")
    parser.add_argument("--ori_block", required=True, help="Input original block.tsv file")
    parser.add_argument("--ori_vcf", required=True, help="Original phased VCF.gz")
    parser.add_argument("--shapeit_vcf", required=True, help="Shapeit VCF.gz")
    parser.add_argument("--out", required=True, help="Output block result file")
    args = parser.parse_args()

    BLOCK_FILE = args.ori_block
    ORI_VCF_GZ = args.ori_vcf
    SHAPEIT_VCF = args.shapeit_vcf
    OUTPUT_FILE = args.out

    # ===================== 读取 block 文件 =====================
    print("Reading block file...")
    block_df = pd.read_csv(BLOCK_FILE, sep="\t")
    block_map = defaultdict(str)
    for _, row in block_df.iterrows():
        chrom = row["chromosome"]
        start = row["from"]
        end = row["to"]
        block_label = f"{chrom}:{start}-{end}"
        for pos in range(start, end + 1):
            block_map[(chrom, pos)] = block_label

    # ===================== 读取原始 VCF =====================
    print("Reading original VCF...")
    with gzip.open(ORI_VCF_GZ, "rt") as f:
        lines = [line.strip() for line in f if not line.startswith("#")]

    ori_rows = []
    for line in lines:
        parts = line.split("\t")
        chrom = parts[0]
        pos = int(parts[1])
        fmt = parts[8].split(":")
        sample = parts[9].split(":")
        gt = dict(zip(fmt, sample)).get("GT", ".")
        if gt in ("0|1", "1|0"):
            ori_rows.append([chrom, pos, gt])

    ori_df = pd.DataFrame(ori_rows, columns=["CHROM", "POS", "ori_GT"])

    # ===================== 读取 shapeit VCF =====================
    print("Reading shapeit VCF...")
    with gzip.open(SHAPEIT_VCF, "rt") as f:
        lines = [line.strip() for line in f if not line.startswith("#")]

    shapeit_rows = []
    for line in lines:
        parts = line.split("\t")
        chrom = parts[0]
        pos = int(parts[1])
        fmt = parts[8].split(":")
        sample = parts[9].split(":")
        fmt_dict = dict(zip(fmt, sample))
        gt = fmt_dict.get("GT", ".")
        if gt in ("0|1", "1|0"):
            shapeit_rows.append([chrom, pos, gt])

    shapeit_df = pd.DataFrame(shapeit_rows, columns=["CHROM", "POS", "shapeit_GT"])

    # ===================== 合并 =====================
    final_df = ori_df.merge(shapeit_df, on=["CHROM", "POS"], how="left")
    final_df["ori_block"] = final_df.apply(
        lambda row: block_map.get((row["CHROM"], row["POS"]), "."), axis=1
    )

    # ===================== GT 比对 =====================
    def compare_gt(ori_gt, shapeit_gt):
        if pd.isna(shapeit_gt) or shapeit_gt == ".":
            return "."
        return "T" if ori_gt == shapeit_gt else "F"

    final_df["GT_match"] = final_df.apply(
        lambda row: compare_gt(row["ori_GT"], row["shapeit_GT"]), axis=1
    )

    # ===================== 交叉表 =====================
    cross_table = pd.crosstab(final_df['ori_block'], final_df['GT_match'], margins=True)

    # ===================== block 标记 =====================
    def label_block(row):
        t = row['T']
        f = row['F']
        if t > 0 and f == 0:
            return 'T'
        elif f > 0 and t == 0:
            return 'F'
        elif t == 0 and f == 0:
            return 'N'
        else:
            return 'MIX'

    cross_table['block_link_relationship'] = cross_table.apply(label_block, axis=1)
    block_result = cross_table[cross_table.index != 'All'].copy()

    # ===================== 输出 =====================
    block_result.to_csv(OUTPUT_FILE, sep="\t")
    print(f"Done! Result saved to: {OUTPUT_FILE}")

if __name__ == "__main__":
    main()