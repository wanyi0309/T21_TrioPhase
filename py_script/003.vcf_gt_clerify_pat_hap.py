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
def vcf_to_dataframe(vcf_path, child_name, mother_name, father_name, child_blk, mat_blk, pat_blk, child_strict_blk):
    rows = []
    child_blocks = load_blocks(child_blk)
    mat_blocks = load_blocks(mat_blk)
    pat_blocks = load_blocks(pat_blk)
    child_strict_blocks = load_blocks(child_strict_blk)

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

        b_child = find_block(chrom, int(pos), gt_child, child_blocks)
        b_mat = find_block(chrom, int(pos), gt_mat, mat_blocks)
        b_pat = find_block(chrom, int(pos), gt_pat, pat_blocks)
        b_child_strict = find_block(chrom, int(pos), gt_child, child_strict_blocks)

        rows.append([chrom, pos, gt_child, gt_mat, gt_pat, b_child, b_mat, b_pat, b_child_strict, type_label])

    f.close()

    # ===================== 【关键：带列名】 =====================
    columns = [
        "chrom", "pos", "chd_gt", "mat_gt", "pat_gt",
        "chd_block", "mat_block", "pat_block", "chd_strict_block", "type"
    ]
    
    df = pd.DataFrame(rows, columns=columns)
    df["pos"] = df["pos"].astype(int)
    
    # ======================
    # 【新增】标记 valid_chd_block
    # ======================
    block_size = df.groupby("chd_block").size()
    df["valid_chd_block"] = df["chd_block"].map(block_size >= 2).map({True: "T", False: "F"})
    
    return df


# ======================================================================
# A1 严格判断 + 父源位置
# ======================================================================
def add_paternal_position(df):
    """
    给 df 增加一列：paternal_in_child → 孩子三倍体中，父源allele在第几位
    只有 A1 \A3类型会计算，其他为 NaN
    """
    def get_paternal_pos(row):
        # 只允许 A1 / A3 计算父源位置
        if row["type"] not in ("A1", "A3"):
            return pd.NA
        
        chd_gt = row["chd_gt"]
        mat_gt = row["mat_gt"]
        pat_gt = row["pat_gt"]

        # 拆分成等位基因列表
        chd = chd_gt.split("|")  # eg['0','0','1']
        mat = mat_gt.split("/")
        pat = pat_gt.split("/")

        mat_allele = mat[0]
        #pat_allele = pat[0]
        pat_allele = pd.NA
        if row["type"] == "A1":
            # A1: 父是纯合 0/0 或 1/1
            pat_allele = pat[0]

        elif row["type"] == "A3":
            # A3: 父是杂合 0/1 → 孩子有的、母亲没有的那个就是父源！
            maternal_set = set(mat)
            child_set = set(chd)
            candidate = list(child_set - maternal_set)
            if len(candidate) == 1:
                pat_allele = candidate[0]

        if pd.isna(pat_allele):
            return pd.NA
        

        # 严格判断：孩子必须 2个来自母亲，1个来自父亲
        cnt_mat = chd.count(mat_allele)
        cnt_pat = chd.count(pat_allele)

        if cnt_mat == 2 and cnt_pat == 1:
            return chd.index(pat_allele)
        else:
            return pd.NA

    # 增加到原 df
    df["paternal_in_child"] = df.apply(get_paternal_pos, axis=1)
    return df


# ======================================================================
# switch 切割 + strict block 优化断点 + 边界安全检查
#  2 个关键要求：
# 1. 区块起止坐标来自 chd_block 本身
# 2. 断点严格不越界：严格和前后SNP位置比较
# ======================================================================
def check_switch_error(df_addpatpos):
    #【关键】只选 valid=T 且 A1 有值的
    df_valid = df_addpatpos[
        (df_addpatpos["valid_chd_block"] == "T") &
        (df_addpatpos["paternal_in_child"].notna()) &
        (df_addpatpos["type"] == "A1")
    ].copy()
    
    df_full = df_addpatpos.copy()

    df_full["step1_switch_check"] = pd.NA
    df_full["step1_newchd_block"] = pd.NA

    if len(df_valid) == 0:
        return df_full

    # --------------------------
    # 1. 区分 pass / switch block
    # --------------------------
    block_consistency = df_valid.groupby("chd_block")["paternal_in_child"].nunique()
    pass_blocks = set(block_consistency[block_consistency == 1].index)
    switch_blocks = set(block_consistency[block_consistency > 1].index)

    # --------------------------
    # 2. 处理 PASS block（整块不变）
    # --------------------------
    for blk in pass_blocks:
        mask = df_full["chd_block"] == blk
        df_full.loc[mask, "step1_switch_check"] = "step1_pass"
        df_full.loc[mask, "step1_newchd_block"] = blk

    # --------------------------
    # 3. 处理 SWITCH block（核心切割）
    # --------------------------
    for blk in switch_blocks:
        sub = df_valid[df_valid["chd_block"] == blk].copy()
        sub = sub.sort_values("pos").reset_index(drop=True)
        if len(sub) < 2:
            continue

        # ======================
        # 修正 1：从 chd_block 字符串解析真实起止
        # ======================
        chrom, block_range = blk.split(":")
        block_start, block_end = block_range.split("-")
        orig_start = int(block_start)
        orig_end = int(block_end)

        segments = []
        current_pat = sub.iloc[0]["paternal_in_child"]
        seg_start = orig_start

        for i in range(1, len(sub)):
            prev_row = sub.iloc[i-1]
            curr_row = sub.iloc[i]
            prev_pat = prev_row["paternal_in_child"]
            curr_pat = curr_row["paternal_in_child"]

            if prev_pat != curr_pat:
                prev_pos = prev_row["pos"]
                curr_pos = curr_row["pos"]

                # ======================
                # 计算 seg_end（前一段结束）
                # ======================
                prev_strict = prev_row["chd_strict_block"]
                if pd.notna(prev_strict) and prev_strict != "NA":
                    _, s_range = prev_strict.split(":")
                    s_start, s_end = s_range.split("-")
                    strict_end = int(s_end)

                    # 修正 2：安全判断：>= prev_pos 才能用
                    if prev_pos <= strict_end < curr_pos:
                        seg_end = strict_end
                    else:
                        seg_end = prev_pos
                else:
                    seg_end = prev_pos

                # 保存前一段
                segments.append((seg_start, seg_end, prev_pat))

                # ======================
                # 计算 seg_start（后一段开始）
                # ======================
                curr_strict = curr_row["chd_strict_block"]
                if pd.notna(curr_strict) and curr_strict != "NA":
                    _, s_range = curr_strict.split(":")
                    s_start, s_end = s_range.split("-")
                    strict_start = int(s_start)

                    # 安全判断：<= curr_pos 才能用
                    if prev_pos < strict_start <= curr_pos:
                        new_seg_start = strict_start
                    else:
                        new_seg_start = curr_pos
                else:
                    new_seg_start = curr_pos

                seg_start = new_seg_start
                current_pat = curr_pat

        # 最后一段
        segments.append((seg_start, orig_end, current_pat))

        # ======================
        # 回填到整个 df_full
        # ======================
        for s, e, pat in segments:
            new_blk = f"{chrom}:{s}-{e}"
            mask = (
                (df_full["chd_block"] == blk) &
                (df_full["pos"] >= s) &
                (df_full["pos"] <= e)
            )
            df_full.loc[mask, "step1_switch_check"] = "step1_switch"
            df_full.loc[mask, "step1_newchd_block"] = new_blk
            
            
    # --------------------------
    # 【最终新增】给每个新block填充统一的pat位置
    # --------------------------
    block_pat_map = {}
    for blk in df_full["step1_newchd_block"].dropna().unique():
        blk_rows = df_full[df_full["step1_newchd_block"] == blk]
        pat_vals = blk_rows[blk_rows["type"] == "A1"]["paternal_in_child"].dropna().unique()
        if len(pat_vals) == 1:
            block_pat_map[blk] = pat_vals[0]

    df_full["step1_newchd_block_pat"] = df_full["step1_newchd_block"].map(block_pat_map)

    return df_full



# ======================================================================
# 主函数
# ======================================================================
def main():
    if len(sys.argv) != 9:
        print("usage：python3 script.py <vcf.gz> <child_name> <mat_name> <pat_name> <child.blk> <mat.blk> <pat.blk> <child_strict.blk>")
        print("example：python3 script.py sample.vcf.gz child mat pat child.blk mat.blk pat.blk child_strict.blk")
        sys.exit(1)

    # 读取输入参数
    vcf_path, child_name, mother_name, father_name, child_blk, mat_blk, pat_blk, child_strict_blk = sys.argv[1:9]

    print("===== 输入参数 =====")
    print("VCF        :", vcf_path)
    print("child      :", child_name)
    print("mother     :", mother_name)
    print("father     :", father_name)
    print("child_blk  :", child_blk)
    print("mat_blk    :", mat_blk)
    print("pat_blk    :", pat_blk)
    print("strict_blk :", child_strict_blk)

    # 1. 读取VCF并注释block、type、strict_block
    print("\n正在处理VCF...")
    df = vcf_to_dataframe(vcf_path, child_name, mother_name, father_name, child_blk, mat_blk, pat_blk, child_strict_blk)

    # 2. 计算A1位点的父源位置
    print("正在计算父源位置...")
    df_addpatpos = add_paternal_position(df)
    
    
    df_valid = df_addpatpos[df_addpatpos["paternal_in_child"].notna()].copy()
    
    df_final  = check_switch_error(df_addpatpos)

    # 3. 输出结果文件
    out_file = f"annotated.tsv"
    df_final.to_csv(out_file, sep="\t", na_rep="NA",index=False)

    print("\n处理完成！")
    print("输出文件:", out_file)
    print("列数:", len(df_final.columns))
    print("总行数:", len(df_final))

if __name__ == "__main__":
    main()