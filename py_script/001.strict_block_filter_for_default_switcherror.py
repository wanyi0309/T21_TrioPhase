######保留了strict中涉及父母纯合子代杂合的blcok，且已将只有一个SNP的blcok去掉，将有switch error的block切分为断开的block

# A1 严格判断正确
# 必须满足：2 个母亲等位 + 1 个父亲等位才会标记父源位置

# switch 切割逻辑
# 无 switch + ≥2 SNP → 保留完整 block
# 有 switch → 切割
# 切割后 **<2 SNP 的片段自动丢弃 **
# 第一段 start = 原始 block start
# 最后一段 end = 原始 block end

import sys
import gzip
import pandas as pd

def load_blocks(block_file):
    blocks = []
    with open(block_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            chrom = parts[1]
            start = int(parts[3])
            end = int(parts[4])
            blocks.append((chrom, start, end))
    return blocks

def find_block(chrom, pos, gt, blocks):
    if "|" not in gt:
        return "NA"
    for b_chrom, b_start, b_end in blocks:
        if b_chrom == chrom and b_start <= pos <= b_end:
            return f"{chrom}:{b_start}-{b_end}"
    return "NA"

# ===================== 【标准化GT】 =====================
def normalize_gt(gt):
    alleles = sorted(gt.replace("|", "/").split("/"))
    return "/".join(alleles)

# ===================== 【父母子分型】 =====================
def get_type(gt_mat, gt_pat, gt_child):
    m = normalize_gt(gt_mat)
    p = normalize_gt(gt_pat)
    c = normalize_gt(gt_child)

    if m == "0/0" and p == "1/1" and c == "0/0/1":
        return "A1"
    if m == "1/1" and p == "0/0" and c == "0/1/1":
        return "A1"

    if m == "1/1" and p == "0/1" and c == "0/1/1":
        return "A3"
    if m == "0/0" and p == "0/1" and c == "0/0/1":
        return "A3"

    if m == "0/1" and p == "0/0" and c == "0/0/1":
        return "B2"
    if m == "0/1" and p == "1/1" and c == "0/1/1":
        return "B2"

    if m == "0/1" and p == "0/1" and (c == "0/1/1" or c == "0/0/1"):
        return "C"

    return "F"

# ======================================================================
# 【核心】直接返回 带列名的 pandas DataFrame
# ======================================================================
def vcf_to_dataframe(vcf_path, child_name, mother_name, father_name, child_blk, mat_blk, pat_blk):
    rows = []
    child_blocks = load_blocks(child_blk)
    mat_blocks = load_blocks(mat_blk)
    pat_blocks = load_blocks(pat_blk)

    f = gzip.open(vcf_path, 'rt') if vcf_path.endswith(".gz") else open(vcf_path, 'r')
    child_idx = mat_idx = pat_idx = None

    for line in f:
        line = line.strip()
        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):
            header = line.split("\t")
            child_idx = header.index(child_name)
            mat_idx = header.index(mother_name)
            pat_idx = header.index(father_name)
            continue

        if child_idx is None:
            continue

        parts = line.split("\t")
        chrom = parts[0]
        pos = parts[1]
        gt_child = parts[child_idx].split(":")[0]
        gt_mat = parts[mat_idx].split(":")[0]
        gt_pat = parts[pat_idx].split(":")[0]

        if "|" not in gt_child:
            type_label = "N"
        else:
            type_label = get_type(gt_mat, gt_pat, gt_child)

        b_child = find_block(chrom, int(pos),gt_child, child_blocks)
        b_mat = find_block(chrom, int(pos),gt_mat, mat_blocks)
        b_pat = find_block(chrom, int(pos),gt_pat, pat_blocks)

        rows.append([chrom, pos, gt_child, gt_mat, gt_pat, b_child, b_mat, b_pat, type_label])

    f.close()

    # ===================== 【关键：带列名】 =====================
    columns = [
        "chrom", "pos", "chd_gt", "mat_gt", "pat_gt",
        "chd_block", "mat_block", "pat_block", "type"
    ]
    
    df = pd.DataFrame(rows, columns=columns)
    df["pos"] = df["pos"].astype(int)
    return df

# ======================================================================
# A1 严格判断 + 父源位置
# ======================================================================
def add_paternal_position(df):
    """
    给 df 增加一列：paternal_in_child → 孩子三倍体中，父源allele在第几位
    只有 A1 类型会计算，其他为 NaN
    """
    def get_paternal_pos(row):
        if row["type"] != "A1":
            return None
        
        chd_gt = row["chd_gt"]
        mat_gt = row["mat_gt"]
        pat_gt = row["pat_gt"]

        # 拆分成等位基因列表
        chd = chd_gt.split("|")  # eg['0','0','1']
        mat = mat_gt.split("/")
        pat = pat_gt.split("/")

        mat_allele = mat[0]
        pat_allele = pat[0]

        # 严格判断：孩子必须 2个来自母亲，1个来自父亲
        cnt_mat = chd.count(mat_allele)
        cnt_pat = chd.count(pat_allele)

        if cnt_mat == 2 and cnt_pat == 1:
            return chd.index(pat_allele)
        else:
            return None

    # 增加到原 df
    df["paternal_in_child"] = df.apply(get_paternal_pos, axis=1)
    return df


def split_blocks_by_switch(df_valid, df_all):
    block_results = []

    for block_name, group in df_valid.groupby("chd_block"):
        group = group.sort_values("pos").reset_index(drop=True)
        chrom = group["chrom"].iloc[0]

        # 原始 block 全部信息
        original_block_df = df_all[df_all["chd_block"] == block_name]
        original_start = original_block_df["pos"].min()
        original_end   = original_block_df["pos"].max()
        original_snp_count = len(original_block_df)

        unique_pat = group["paternal_in_child"].unique()
        no_switch = len(unique_pat) == 1

        # ==============================
        # 情况1：无 switch → 完整保留
        # ==============================
        if no_switch:
            if original_snp_count >= 2:
                # 多加一个标记：is_switch_cut = False
                block_results.append((
                    chrom, original_start, original_end,
                    unique_pat[0], f"{chrom}:{original_start}-{original_end}",
                    False  # 没有被切割过
                ))
            continue

        # ==============================
        # 情况2：有 switch → 切割
        # ==============================
        current_pat = group["paternal_in_child"].iloc[0]
        block_start = original_start

        for i in range(1, len(group)):
            pat = group["paternal_in_child"].iloc[i]

            if pat != current_pat:
                block_end = group["pos"].iloc[i-1]
                segment_df = original_block_df[
                    (original_block_df["pos"] >= block_start) &
                    (original_block_df["pos"] <= block_end)
                ]
                if len(segment_df) >= 2:
                    # 被切割过 → True
                    block_results.append((
                        chrom, block_start, block_end, current_pat,
                        f"{chrom}:{block_start}-{block_end}",
                        True  # 被switch切割过
                    ))

                current_pat = pat
                block_start = group["pos"].iloc[i]

        # ==============================
        # 最后一段
        # ==============================
        block_end = original_end
        segment_df = original_block_df[
            (original_block_df["pos"] >= block_start) &
            (original_block_df["pos"] <= block_end)
        ]
        if len(segment_df) >= 2:
            block_results.append((
                chrom, block_start, block_end, current_pat,
                f"{chrom}:{block_start}-{block_end}",
                True  # 被切割过
            ))

    # 列名也要加一个：is_switch_cut
    df_blocks = pd.DataFrame(block_results, columns=[
        "chrom", "start", "end", "paternal_idx", "new_block", "is_switch_cut"
    ])
    
    # 排序修复（数字排序）
    df_blocks["start"] = df_blocks["start"].astype(int)
    df_blocks["end"] = df_blocks["end"].astype(int)
    df_blocks = df_blocks.sort_values(["chrom","start"], ascending=True).reset_index(drop=True)

    return df_blocks


# ======================================================================
# 【最终输出格式】你要的表格
# ======================================================================
def build_final_output(df_new_blocks, df_all, sample_name):
    out_rows = []

    for _, r in df_new_blocks.iterrows():
        chrom = r["chrom"]
        s = r["start"]
        e = r["end"]
        is_switch = r["is_switch_cut"]

        # 计数：这段区间内 phased 的 SNP 数量（chd_block != NA）
        cnt = len(df_all[
            (df_all["chrom"] == chrom) &
            (df_all["pos"] >= s) &
            (df_all["pos"] <= e) &
            (df_all["chd_block"] != "NA")
        ])

        phase_set = s
        out_rows.append([
            sample_name, chrom, phase_set, s, e, cnt, is_switch
        ])

    df_out = pd.DataFrame(out_rows, columns=[
        "#sample", "chromosome", "phase_set", "from", "to", "variants", "is_switch_cut"
    ])
    return df_out


# ======================================================================
# 主函数
# ======================================================================
def main():
    if len(sys.argv) != 8:
        print("usage：python3 vcf_gt_block_merger.py <vcf.gz> <child_name> <mat_name> <pat_name> <child_block.file> <mat_block.file> <pat_block.file>")
        sys.exit(1)
       
    vcf_path, child_name, mother_name, father_name, child_blk, mat_blk, pat_blk = sys.argv[1:8]
    print("input：")
    print("VCF         :", vcf_path)
    print("child       :", child_name)
    print("mother      :", mother_name)
    print("father      :", father_name)

    df = vcf_to_dataframe(vcf_path, child_name, mother_name, father_name, child_blk, mat_blk, pat_blk)
    df_addpatpos = add_paternal_position(df)
    df_valid = df_addpatpos[df_addpatpos["paternal_in_child"].notna()].copy()
    df_new_blocks = split_blocks_by_switch(df_valid, df_addpatpos)

    # ===================== 最终输出 =====================
    df_final = build_final_output(df_new_blocks, df_addpatpos, child_name)
    #print("\n【最终输出表格】")
    #print(df_final.to_string(index=False))

    # 保存成文件
    out_fn = f"strick_blocks_for_default_switch_check.tsv"
    df_final.to_csv(out_fn, sep="\t", index=False)
    print(f"\n输出文件：{out_fn}")

    return df_final

if __name__ == "__main__":
    df = main()