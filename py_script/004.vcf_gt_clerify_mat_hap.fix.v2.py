import pandas as pd
#import matplotlib.pyplot as plt
import os
import numpy as np

# ===========================================================================================
# ===========================================================================================
# step1结果筛选为step2做准备 
# ===========================================================================================
# ===========================================================================================

# 读step1之后的annotate文件
df = pd.read_csv(
    "annotated.tsv",
    sep="\t",          # 必须是 tab 分隔
    na_values="NA",    # 把 NA 识别成缺失值
    dtype={"pos": int} # 保证 pos 是数字
)

# 将step1_newchd_block和mat_block拼到一起，为了后面按chd大block中的mat小block分组做match
df["step1_newchd_block_in_mat"] = (
    df["step1_newchd_block"].astype(str) + "_" + df["mat_block"].astype(str)
)

# 筛选A1>=5 B2+C>=5 的block。A1>=5确保step1判断pat的准确性，B2+C>=5 确保有足够的SNP进行step2
step1_newchd_block_stat = df.groupby("step1_newchd_block")["type"].value_counts().unstack(fill_value=0)
filtered_blocks = step1_newchd_block_stat[
    (step1_newchd_block_stat["A1"] >= 5) & 
    (step1_newchd_block_stat["B2"] + step1_newchd_block_stat["C"] >= 5)
].index.tolist() 

# ===========================================================================================
# ====== 【重要新增】在筛选B2/C之前，先对所有合格block的全部SNP做全局投票 ======
# ===========================================================================================
# 切割GT（原有函数）
def split_gt(gt):
    return str(gt).replace("|", "/").split("/")

# 获取mat在chd的列信息（原有函数）
def get_maternal_cols(pat_col):
    all_cols = {0, 1, 2}
    mat_cols = sorted(all_cols - {int(pat_col)})
    return mat_cols[0], mat_cols[1]

# ====== 【全新函数】对 A1/A3/B/C 所有类型做投票 ======
# 规则：
# A1/A3：母本纯合 → 去掉父本等位后，孩子剩余两个单倍型必须也纯合 → 合格=T，不合格=F
# B/C：沿用原有逻辑，返回 A/B/F
# 其他：返回 F
def vote_all_type_row(row):
    chd = split_gt(row["chd_gt"])
    mat = split_gt(row["mat_gt"])
    p   = int(row["step1_newchd_block_pat"])
    a, b = get_maternal_cols(p)
    ca, cb = chd[a], chd[b]
    m0, m1 = mat[0], mat[1]
    snp_type = str(row["type"])

    # ---------------- A1 / A3 逻辑 ----------------
    if snp_type in ["A1", "A3"]:
        # 母本必须纯合
        if m0 != m1:
            return "F"
        # 孩子去掉父本后，剩余两个必须相同（纯合）
        if ca == cb  and ca == m0:
            return "T"
        else:
            return "F"

    # ---------------- B / C 逻辑（和原vote_b2_row完全一致） ----------------
    elif snp_type in ["B", "C", "B2"]:
        if {m0, m1} != {"0", "1"}:
            return "F"
        if ca == cb:
            return "F"
        scoreA = (ca == m0) + (cb == m1)
        scoreB = (ca == m1) + (cb == m0)
        if scoreA > scoreB:
            return "A"
        elif scoreB > scoreA:
            return "B"
        else:
            return "F"

    # 其他类型
    else:
        return "F"

# ====== 生成【全部合格block】的临时表，计算全局vote ======
all_valid_df = df[df["step1_newchd_block"].isin(filtered_blocks)].copy()
all_valid_df["global_vote"] = all_valid_df.apply(vote_all_type_row, axis=1)

# ====== 把 global_vote 合并回原始 df ======
df = df.merge(
    all_valid_df[["chrom", "pos", "global_vote"]],
    on=["chrom", "pos"],
    how="left"
)
# ===========================================================================================
# ====== 【新增结束】 ======
# ===========================================================================================

# step2主要使用B2+C的SNP。同时这些SNP要mat也是phased
filter_df = df[
    df["step1_newchd_block"].isin(filtered_blocks) &  # 保留目标block
    df["type"].isin(["B2", "C"])                      # 只保留B2和C
].copy()
filter_df = filter_df.sort_values(["chrom", "pos"]).reset_index(drop=True)
filter_df = filter_df[filter_df["mat_block"].notna()] # 同时这些SNP要mat也是phased

stat = filter_df.groupby("step1_newchd_block_in_mat")["type"].value_counts().unstack(fill_value=0)
print(stat)

# 先统计每个小block的B2+C数量
small_stat = filter_df.groupby("step1_newchd_block_in_mat")["type"].apply(
    lambda x: x.isin(["B2", "C"]).sum()
)
# 保留B2+C >=5的小block
filter_df = filter_df[filter_df["step1_newchd_block_in_mat"].isin(small_stat[small_stat >= 5].index)]


# ===========================================================================================
# ===========================================================================================
# step2中的几个小函数
# ===========================================================================================
# ===========================================================================================

# 切割GT
def split_gt(gt):
    return str(gt).replace("|", "/").split("/")

# 获取mat在chd的列信息
def get_maternal_cols(pat_col):
    """给定父源在 chd_gt 的列号 p, 返回另外两列 (a, b), a<b"""
    all_cols = {0, 1, 2}
    mat_cols = sorted(all_cols - {int(pat_col)})
    return mat_cols[0], mat_cols[1]

# 对单个 B2/C 位点投票
def vote_b2_row(row):
    """
    对单个 B2 位点投票: 比较 (chd_a, chd_b) 与 (mat0, mat1)
    返回:
       'A'  : 支持配置 chd_a<->mat0, chd_b<->mat1
       'B'  : 支持配置 chd_a<->mat1, chd_b<->mat0
       'F'  : tie / 信息不足
    """
    chd = split_gt(row["chd_gt"])
    mat = split_gt(row["mat_gt"])
    p   = int(row["step1_newchd_block_pat"])
    a, b = get_maternal_cols(p)

    ca, cb = chd[a], chd[b]
    m0, m1 = mat[0], mat[1]

    # mat 必须是 0/1 (B2 保证)
    if {m0, m1} != {"0", "1"}:
        return "F"
    # 母亲两个等位一样就无法区分
    if ca == cb:
        return "F"

    # 配置 A: ca<->m0, cb<->m1 的匹配数
    # 对于scoreA 如果ca=m0 (ca == m0)结果是true (ca == m0) + (cb == m1)=2
    # 对于scoreA 如果ca!=m0 (ca == m0)结果是false (ca == m0) + (cb == m1)=0
    # 对于scoreB 亦然 
    scoreA = (ca == m0) + (cb == m1)
    scoreB = (ca == m1) + (cb == m0)
    if scoreA > scoreB:  return "A"
    if scoreB > scoreA:  return "B"
    return "F"

# 滑窗判断A比例 窗口大小自定义
# votes是block中SNP经过vote_b2_row判断后的一个list
def rolling_majority(votes, w):
    """votes 是 'A'/'B'/'T' 列表, 返回每个位点的 A 比例 (只在 A+B 中算)"""
    n = len(votes)
    out = np.full(n, np.nan) # 先建一个全是空值的结果数组
    arr = np.array(votes) # 转成 numpy 数组方便计算
    for i in range(n):
        lo = max(0, i - w // 2) # 窗口左边界：不能小于 0
        hi = min(n, i + w // 2 + 1)  # 窗口右边界：不能超过总长度
        win = arr[lo:hi]  # 取出窗口内的所有投票
        a = int((win == "A").sum()) # 窗口里有多少 A
        b = int((win == "B").sum()) # 窗口里有多少 B
        if a + b >= 3: # 只有 有效投票数(A+B) 达到阈值（至少3个），才计算比例
            out[i] = a / (a + b)
    return out

# ===========================================================================================
# ===========================================================================================
# step2中的几个阈值或参数
# ===========================================================================================
# ===========================================================================================
WINDOW = 5          # 滑动窗口大小
#MIN_B2_FOR_CALL = 5  # 一个 block / 一个子段至少需要多少个 B2 才敢下结论
#SWITCH_MARGIN = 0.5  # 滑窗比例在 >0.5 算稳定 A, <0.5 算稳定 B, 中间算不确定
MIN_SEG_LEN = 3      # 一个block中至少要有几个 连续A或B的SNP才能被认为是一个小switch,否则视为噪声
F_RATIO_THRESHOLD = 0.8 #block中F的比例


# ===========================================================================================
# ===========================================================================================
# 主循环: 对每个 step1_newchd_block 做 step2
# ===========================================================================================
# ===========================================================================================
results_per_snp = []   # 回填 SNP 级别信息
block_summaries = []   # block 级别汇总
vote_storage = {} # 全局存储 vote

unique_blocks = sorted(filter_df["step1_newchd_block_in_mat"].unique())

# =============================================================

for blk in unique_blocks:
    blk_df_full = filter_df[filter_df["step1_newchd_block_in_mat"] == blk].copy()
    blk_df_full = blk_df_full.sort_values("pos").reset_index(drop=True)
    
    n_b2c = len(blk_df_full)
    n_b2= len(blk_df_full[blk_df_full["type"]=="B2"])
    n_c= len(blk_df_full[blk_df_full["type"]=="C"])
    blk_df_full["vote"] = blk_df_full.apply(vote_b2_row, axis=1)
    votes = blk_df_full["vote"].tolist()
    
    # 把当前 block 的所有 vote 存起来（所有 A/B/F 都存，绝不丢）
    for _, row in blk_df_full.iterrows():
        vote_storage[(row["chrom"], row["pos"])] = row["vote"]

    nA = votes.count("A")
    nB = votes.count("B")
    nF = votes.count("F")
    total = len(votes)

    # ====================== 1. F ≥ 80% 直接丢弃整个 block ======================
    #F_RATIO_THRESHOLD = 0.8
    if  nF / total >= F_RATIO_THRESHOLD:
        block_summaries.append({
            "step1_newchd_block": blk,
            "n_snp_total": len(blk_df_full),
            #"n_b2c": n_b2c,
            "n_b2": n_b2,
            "n_c": n_c,
            "decision": "SKIP_F_TOO_MANY",
            "configA_votes": nA,
            "configB_votes": nB,
            "tie_votes": nF,
            "ratioA_percent": round(nA/total*100,2) if total>0 else 0,
            "ratioB_percent": round(nB/total*100,2) if total>0 else 0,
            "ratioF_percent": round(nF/total*100,2) if total>0 else 0,
            "n_switch": 0,
            "step2_segments": "NA",
        })
        
        # ===================== 修复：保留所有位点，只把 step2 设为 NA，不跳过 =====================
        sub = blk_df_full.copy()
        # 所有 step2 注释清空为 NA，但是 vote 完全保留
        # sub["step2_newblock"] = np.nan
        # sub["step2_mat_config"] = np.nan
        sub["step2_newblock"] = pd.Series(
            pd.NA,
            index=sub.index,
            dtype="string"
        )

        sub["step2_mat_config"] = pd.Series(
            pd.NA,
            index=sub.index,
            dtype="string"
        )
        
        sub["chd_col_from_mat0"] = np.nan
        sub["chd_col_from_mat1"] = np.nan
        sub["chd_col_from_pat"] = np.nan
        results_per_snp.append(sub)
        continue
    # ===========================================================================

    # 正常处理
    # 只提取 A/B 和它们在 blk_df_full 中的原始索引
    valid_indices = []
    valid_votes = []
    for idx, v in enumerate(votes):
        if v in ("A", "B"):
            valid_indices.append(idx)
            valid_votes.append(v)
    
    valid_ratioA = rolling_majority(valid_votes, WINDOW)
    
    valid_labels = ["U"] * len(valid_votes)
    for i in range(len(valid_votes)):
        rat = valid_ratioA[i]
        # if not np.isnan(rat):
        #     if rat > SWITCH_MARGIN:
        #         valid_labels[i] = "A"
        #     elif rat < 1 - SWITCH_MARGIN:
        #         valid_labels[i] = "B"
        if rat > 0.6:
            valid_labels[i] = "A"
        elif rat < 0.4:
            valid_labels[i] = "B"
        else:
            valid_labels[i] = valid_votes[i]  # 模糊区间，信任原始投票
        
                
    # 判断是否有足够长的A/B段（基于纯净标签）
    nA_valid = sum(1 for x in valid_labels if x == "A")
    nB_valid = sum(1 for x in valid_labels if x == "B")
    has_A_region = nA_valid >= MIN_SEG_LEN
    has_B_region = nB_valid >= MIN_SEG_LEN

    # 正向填充（前向填充）：把中间的 U 用最近的 A/B 填满
    if has_A_region and has_B_region:
        filled = valid_labels.copy() # 1. 复制一份标签，不改动原始数据
#         last = None # 2. 记录“上一个有效标记”（A/B），一开始是空的
#         # 3. 从头到尾 正着走一遍（从左到右）
#         for i in range(len(filled)):
#             # 4. 如果当前是 A 或 B → 有效标记
#             if filled[i] in ("A", "B"):
#                 last = filled[i]
#             # 5. 如果当前是 U（不确定）→ 用“上一个有效标记”填充
#             else:
#                 # 如果前面有过 A/B → 填 last
#                 # 如果前面从来没有 → 保持 U 不变
#                 filled[i] = last if last is not None else filled[i]
        
#         # 1. 清空记录：准备从右边开始重新扫描
#         last = None
#         # 2. 关键！！！
#     # range(len(filled) - 1, -1, -1) 意思是：从 最后一个位置 → 第 0 个位置 倒着走！！！反向走！！！
#         for i in range(len(filled) - 1, -1, -1):
#             if filled[i] in ("A", "B"):
#                 last = filled[i]
#             elif filled[i] == "U":
#                 filled[i] = last if last is not None else "A"
        
        segs = []
        # 连续片段切割（把一长串 A/B 切成一段一段）
        cur_label = filled[0]
        cur_start = 0
        for i in range(1, len(filled)):
            if filled[i] != cur_label:
                segs.append((cur_start, i - 1, cur_label))
                cur_start = i
                cur_label = filled[i]
        segs.append((cur_start, len(filled) - 1, cur_label))

        # 合并太短的片段（去噪声）
        merged = []
        for s, e, lab in segs:
            length = e - s + 1
            # 少于 MIN_SEG_LEN 个 SNP 的段 = 噪声
            if length < MIN_SEG_LEN and merged:
                #如果这段太短 + 前面已经有段了 把这段短的，直接合并到上一段！
                ps, pe, pl = merged[-1]
                merged[-1] = (ps, e, pl)
            else:
                # 够长 → 保留
                merged.append((s, e, lab))
        segs = merged
        decision = "SWITCH"
        
    else:
        dominant = "A" if nA >= nB else "B"
        segs = [(0,len(valid_votes)-1, dominant)]
        decision = f"CONSISTENT_{dominant}"
    
    chrom = blk.split(":")[0]
    step2_segments_str = []

    
    for seg in segs:
        s_in_valid, e_in_valid, lab = seg
    
        # 转成原始 index
        raw_s = valid_indices[s_in_valid]
        raw_e = valid_indices[e_in_valid]
    
        # 原始 pos
        left_pos = blk_df_full["pos"].iloc[raw_s]
        right_pos = blk_df_full["pos"].iloc[raw_e]
        
        new_blk = f"{chrom}:{left_pos}-{right_pos}"
        step2_segments_str.append(f"{chrom}:{left_pos}-{right_pos}({lab},n_B2C={e_in_valid - s_in_valid + 1})")

        # 全部位点都拿出来，包括 F！！！
        sub = blk_df_full[(blk_df_full["pos"] >= left_pos) & (blk_df_full["pos"] <= right_pos)].copy()

        # 只给 A/B 加 step2 注释，F 不加（但保留 F 行）
        # sub["step2_newblock"] = np.nan
        # sub["step2_mat_config"] = np.nan
        sub["step2_newblock"] = pd.Series(
            pd.NA,
            index=sub.index,
            dtype="string"
        )

        sub["step2_mat_config"] = pd.Series(
            pd.NA,
            index=sub.index,
            dtype="string"
        )
        sub["chd_col_from_mat0"] = np.nan
        sub["chd_col_from_mat1"] = np.nan
        sub["chd_col_from_pat"] = np.nan

        # 给 A/B 赋值
        mask_ab = sub["vote"].isin(["A", "B"])
        pat_col = int(blk_df_full["step1_newchd_block_pat"].iloc[0])
        a_col, b_col = get_maternal_cols(pat_col)

        if lab == "A":
            sub.loc[mask_ab, "chd_col_from_mat0"] = a_col
            sub.loc[mask_ab, "chd_col_from_mat1"] = b_col
        else:
            sub.loc[mask_ab, "chd_col_from_mat0"] = b_col
            sub.loc[mask_ab, "chd_col_from_mat1"] = a_col

        sub.loc[mask_ab, "step2_newblock"] = new_blk
        sub.loc[mask_ab, "step2_mat_config"] = lab
        sub.loc[mask_ab, "chd_col_from_pat"] = pat_col

        results_per_snp.append(sub)
    # =================================================================================

    block_summaries.append({
        "step1_newchd_block": blk,
        "n_snp_total": len(blk_df_full),
        #"n_b2": n_b2c,
        "n_b2": n_b2,
        "n_c": n_c,
        "decision": decision,
        "configA_votes": nA,
        "configB_votes": nB,
        "tie_votes": nF,
        "ratioA_percent": round(nA/total*100,2) if total>0 else 0,
        "ratioB_percent": round(nB/total*100,2) if total>0 else 0,
        "ratioF_percent": round(nF/total*100,2) if total>0 else 0,
        "n_switch": max(0, len(segs) - 1),
        "step2_segments": " | ".join(step2_segments_str),
    })
    

if results_per_snp:
    snp_df = pd.concat(results_per_snp, ignore_index=True)

    # 只合并 step2 注释，不含 vote
    merge_cols = [
        "chrom", "pos",
        "step2_newblock", "step2_mat_config",
        "chd_col_from_mat0", "chd_col_from_mat1", "chd_col_from_pat"
    ]
    out_df = df.merge(snp_df[merge_cols], on=["chrom", "pos"], how="left")

    # ===================== 真正正确的 vote 回填 =====================
    out_df["vote"] = out_df.set_index(["chrom", "pos"]).index.map(vote_storage)
    # =================================================================

    # step2_filter 标注
    out_df["step2_used_block"] = out_df["step1_newchd_block"].isin(filtered_blocks).map({True: "T", False: "F"})

    out_df.to_csv("annotated.step2.tsv", sep="\t", index=False, na_rep="NA")
    #out_df.to_csv("test1", sep="\t", index=False, na_rep="NA")
    
summary_df = pd.DataFrame(block_summaries)
summary_df.to_csv("annotated.step2.block.tsv", sep="\t", index=False, na_rep="NA")
