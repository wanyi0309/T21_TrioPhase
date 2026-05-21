import pandas as pd
import gzip
import argparse

# ===================== 命令行传参 =====================
parser = argparse.ArgumentParser(description="Build final child haplotype VCF")
parser.add_argument("--keep", required=True, help="filterSNP.add_A1A3hapinfo.txt")
parser.add_argument("--stat", required=True, help="step3_filter.stat")
parser.add_argument("--vcf", required=True, help="original vcf.gz")
parser.add_argument("--child", required=True, help="child sample ID (e.g., T21E11)")
parser.add_argument("--out", required=True, help="output VCF")
args = parser.parse_args()

# ===================== 读取文件 =====================
keep_df = pd.read_csv(args.keep, sep="\t", dtype={"chrom": str, "pos": int})
stat_df = pd.read_csv(args.stat, sep="\t", dtype=str)

# 构建 SNP → block 映射
snp_block_map = {}
snp_chdgt_map = {}
for _, r in keep_df.iterrows():
    chrom, pos = r["chrom"], r["pos"]
    snp_block_map[(chrom, pos)] = r["imp_step2_newblock"]
    snp_chdgt_map[(chrom, pos)] = str(r["chd_gt"])

# 构建 block → (final_pat, final_mat0, final_mat1)
block_hap_map = {}
for _, r in stat_df.iterrows():
    block = r["imp_step2_newblock"]
    p = r["final_pat"]
    m0 = r["final_mat0"]
    m1 = r["final_mat1"]
    block_hap_map[block] = (p, m0, m1)

# ===================== 处理 VCF =====================
with gzip.open(args.vcf, "rt") as fin:
    header_lines = []
    sample_col_idx = None
    out_header = None

    for line in fin:
        line = line.rstrip("\n")
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            parts = line.split("\t")
            # 找到孩子样本所在列
            sample_col_idx = parts.index(args.child)
            # 构建输出表头
            out_header = (
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + args.child
            )
            break

# ===================== 处理 VCF =====================
with gzip.open(args.vcf, "rt") as fin, open(args.out, "w", encoding="utf-8") as fout:
    fout.write("##fileformat=VCFv4.2\n")
    fout.write(out_header + "\n")

    for line in fin:
        line = line.rstrip("\n")
        if line.startswith("#"):
            continue

        parts = line.split("\t")
        chrom = parts[0]
        pos = int(parts[1])
        vid = parts[2]
        ref = parts[3]
        alt = parts[4]
        
        # ==========================================
        qual = "."
        filt = "."
        info = "."
        # ================================================================

        key = (chrom, pos)
        if key not in snp_block_map:
            continue
        
        if pos<13000000:
            continue

        # 获取 block 和单倍型索引
        block = snp_block_map[key]
        if block not in block_hap_map:
            continue
        ip, im0, im1 = block_hap_map[block]

        # 预期的孩子GT（来自filterSNP文件）
        expected_chd_gt = snp_chdgt_map[key]
        
        # ===================== 【新增】从VCF中取出真实GT并校验 =====================
        try:
            sample_fields = parts[sample_col_idx].split(":")
            real_gt = sample_fields[0]  # 取第一个字段，即GT
        except:
            real_gt = "."

        # 比对：不一样就警告
        if real_gt != expected_chd_gt:
            print(f"警告：{chrom}:{pos} GT不匹配 | 原VCF={real_gt} | 预期={expected_chd_gt}")
        # ==========================================================================

        # 原始孩子GT（使用filter文件中的，不受VCF影响）
        chd_gt = expected_chd_gt
        hap = chd_gt.split("|")
        if len(hap) != 3:
            continue

        # 转数字索引
        try:
            ip = int(float(ip))
            im0 = int(float(im0))
            im1 = int(float(im1))
        except:
            continue

        # 新GT：pat | mat0 | mat1
        new_gt = f"{hap[ip]}|{hap[im0]}|{hap[im1]}"

        # 输出基础 VCF
        outline = (
            f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{qual}\t{filt}\t{info}\tGT\t{new_gt}"
        )
        fout.write(outline + "\n")

print("完成！输出文件：", args.out)
