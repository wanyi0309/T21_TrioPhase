import pandas as pd
import numpy as np
import re
import argparse

# ===================== 命令行传参配置 =====================
def main():
    parser = argparse.ArgumentParser(description="SNP筛选与单倍型区块整合统计脚本")
    parser.add_argument("--infile", type=str, required=True, help="输入tsv文件: annotated.step2.tsv")
    parser.add_argument("--out_snp", type=str, required=True, help="过滤后SNP结果输出路径")
    parser.add_argument("--out_stat", type=str, required=True, help="block统计结果输出路径")
    parser.add_argument("--chr_name", type=str, default="chr21", help="染色体名称，默认chr21")
    parser.add_argument("--mat_link_file", type=str, required=True, help="母亲block关联关系文件路径")
    args = parser.parse_args()

    df = pd.read_csv(
        args.infile,
        sep="\t",          # 必须是 tab 分隔
        na_values="NA",    # 把 NA 识别成缺失值
        dtype={"pos": int} # 保证 pos 是数字
    )

    # 基础过滤
    ## 过滤初始SNPtype判断为F的（这里F的主要是chd SNP不符合家系遗传规律）
    df1 = df[df["type"]!="F"]
    ## 过滤["valid_chd_block"]!="F" （这类主要是block中只有一个SNP）
    df2 = df1[df1["valid_chd_block"]!="F"]
    ## ["step2_used_block"]=="T" 保留判断了pat的A1>=5&(B2+C>=5)的valid_block
    df3=df2[df2["step2_used_block"]=="T"]
    ## 过滤掉["global_vote"]为"F"的（主要是在确定pat的信息后，剩下的两个来自mat的allele不符合遗传规律）
    df4=df3[df3["global_vote"]!="F"]
    ## 过滤掉["vote"]不等于["step2_mat_config"]的SNP（这一部分是指在判断chd∩mat_block的两个mathap时，判断时我们接受一串连续的hap组合中有1-2个小的跳转错误。判断完后，统计发现这部分的SNP比较少，考虑到保留准确度更高的，因而过滤掉）
    # 条件1：两列值相等
    cond1 = df4["vote"] == df4["step2_mat_config"]
    # 条件2：两列 同时 为缺失值（NA）
    cond2 = df4["vote"].isna() & df4["step2_mat_config"].isna()
    # 满足任意一个条件就保留
    df5 = df4[cond1 | cond2]

    ## 通过上面的过滤，基本已经去除掉可信度不高的SNP
    # 添加A1、A3类SNP信息，更新最碎的block的边界，过滤只有A1、A3的block（其实就是step2中的SKIP_F_TOO_MANY block）

    # 需要填充的 5 个目标列
    cols = [
        "step2_newblock",
        "step2_mat_config",
        "chd_col_from_mat0",
        "chd_col_from_mat1",
        "chd_col_from_pat"
    ]

    df5 = df5.copy() 

    # 先复制原值到 imp_ 开头的列
    for c in cols:
        df5[f"imp_{c}"] = df5[c].copy()

    # 按 pos 排序（必须按位置排序才能找最近）
    df5 = df5.sort_values(["step1_newchd_block", "pos"]).reset_index(drop=True)

    # ---------------------
    for c in cols:
        df5[f"imp_{c}"] = df5.groupby("step1_newchd_block", group_keys=False)[f"imp_{c}"].apply(
            lambda x: x.fillna(method="ffill").fillna(method="bfill")
        )
    # 先往下填空缺 → 再往上填空缺 最终 = 同组内离得最近的有值 SNP 来填充
    # ffill = forward fill = 向前填充 / 向下填充
    # bfill = backward fill = 向后填充 / 向上填充
    # 这里的前后，是按照 pos 从小到大排序后 的顺序

    # 对每一组 imp_step2_newblock
    # 找到这一组里所有 SNP 中最小的 pos → 作为新 block 的 start
    # 找到这一组里所有 SNP 中最大的 pos → 作为新 block 的 end
    # 把 block 名字统一改成：
    # chr21: 最小 pos - 最大 pos
    # --------------------------
    # 1. 按 imp_step2_newblock 分组，计算每个 block 的真实 min_pos 和 max_pos
    # --------------------------
    block_range = df5.groupby("imp_step2_newblock")["pos"].agg(
        block_start="min",  # 这个block里最小的pos
        block_end="max"     # 这个block里最大的pos
    ).reset_index()

    # 拼接成新的 block 字符串：chrXX:start-end
    block_range["new_block_name"] = (
        f"{args.chr_name}:" 
        + block_range["block_start"].astype(str) 
        + "-" 
        + block_range["block_end"].astype(str)
    )

    # 建立映射关系
    block_map = block_range.set_index("imp_step2_newblock")["new_block_name"]

    # --------------------------
    # 2. 替换掉原来的 imp_step2_newblock，变成统一放大后的区间
    # --------------------------
    df5["imp_step2_newblock"] = df5["imp_step2_newblock"].map(block_map)

    stat = df5[df5["valid_chd_block"]=="T"].groupby("step1_newchd_block")["type"].value_counts().unstack(fill_value=0)
    stat.head()

    # 1. 先算出每个 block：B2 + C 是否等于 0
    stat["B2_plus_C"] = stat["B2"] + stat["C"]

    # 2. 找出需要保留的 block（B2+C > 0）
    keep_blocks = stat[stat["B2_plus_C"] > 0].index.tolist()

    # 3. 过滤 df5，只保留这些 block
    df5_filtered = df5[df5["step1_newchd_block"].isin(keep_blocks)].copy()

    df5_filtered.to_csv(args.out_snp, sep="\t", index=False, na_rep="NA")
    
    # 2. 定义需要检查的列
    cols_to_check = [
        'chd_block', 
        'mat_block', 
        'imp_chd_col_from_mat0', 
        'imp_chd_col_from_mat1', 
        'imp_chd_col_from_pat'
    ]

    # 3. 按 imp_step2_newblock 分组
    results = []

    grouped = df5_filtered.groupby("imp_step2_newblock")

    for name, group in grouped:
        group_b2c = group[group["type"].isin(["B2", "C"])].copy()
        info = {
            "imp_step2_newblock": name,
            "n_snp": len(group),
            "n_a1":len(group[group["type"]=="A1"]),
            "n_A3":len(group[group["type"]=="A3"]),
             "n_B2":len(group[group["type"]=="B2"]),
             "n_C":len(group[group["type"]=="C"]),
            "step1_newchd_block": group['step1_newchd_block'].iloc[0] if 'step1_newchd_block' in group.columns else "NA"
        }
        
        is_pure = True
        for col in cols_to_check:
            # 【关键改进】：提取该列中非缺失的值
            # 我们排除掉 pd.isna 以及 字符串形式的 "NA", "nan"
            valid_series = group_b2c[col].replace(['NA', 'nan', 'None'], np.nan).dropna()
            
            # 获取去重后的有效值
            unique_vals = valid_series.unique()
            
            if len(unique_vals) == 0:
                # 如果全是 NA，标记为 NA
                info[col] = "NA"
            elif len(unique_vals) == 1:
                # 如果只有一个有效值，处理掉小数点并记录
                val = unique_vals[0]
                try:
                    # 尝试转为整数（去除 .0）
                    info[col] = str(int(float(val)))
                except:
                    info[col] = str(val)
            else:
                # 如果存在多个不同的有效值，标记冲突
                info[col] = f"MULTIPLE({list(unique_vals)})"
                is_pure = False
        
        info["is_consistent"] = "Yes" if is_pure else "No"
        results.append(info)

    # 4. 转换为 DataFrame 并排序
    summary_df = pd.DataFrame(results)

    def get_start(s):
        try: return int(str(s).split(':')[-1].split('-')[0])
        except: return 0

    summary_df['start_pos'] = summary_df['imp_step2_newblock'].apply(get_start)
    summary_df = summary_df.sort_values('start_pos').drop(columns=['start_pos'])

    
    # 1. 读取母亲block的link关系文件
    mat_link_df = pd.read_csv(args.mat_link_file, sep="\t")
    mat_map = mat_link_df.set_index("ori_block")[["T", "F"]].to_dict("index")

    # 2. 解析函数：按 T/F 数值求和判断
    def get_matblock_link(mat_block_val):
        if pd.isna(mat_block_val) or mat_block_val == "NA":
            return "NA"

        # 提取 block 列表
        if mat_block_val.startswith("MULTIPLE("):
            match = re.search(r"MULTIPLE\(\[(.*?)\]\)", mat_block_val)  # 用正则提取 MULTIPLE([ 内容 ]) 中间的内容
            if not match:
                return "NA"
            block_str = match.group(1).replace("'", "").strip() # 拿到中间的字符串，去掉单引号，去掉首尾空格
            block_list = [b.strip() for b in block_str.split(",")] # 按逗号分割成列表，每个元素再去空格
        else:
            block_list = [mat_block_val.strip()]

        # 求和 T 和 F
        sum_T = 0
        sum_F = 0
        for blk in block_list:
            if blk in mat_map:
                sum_T += mat_map[blk]["T"]
                sum_F += mat_map[blk]["F"]

        # 比较判断
        if sum_T > sum_F:
            return "T"
        elif sum_F > sum_T:
            return "F"
        else:
            return "NA"

    # 3. 加到 summary_df
    summary_df["matblock_link"] = summary_df["mat_block"].apply(get_matblock_link)

    def assign_final_cols(row):
        # 固定不变
        final_pat = row['imp_chd_col_from_pat']
        
        # 根据 matblock_link 判断是否交换母亲的两个列
        if row['matblock_link'] == 'T':
            final_mat0 = row['imp_chd_col_from_mat0']
            final_mat1 = row['imp_chd_col_from_mat1']
        elif row['matblock_link'] == 'F':
            # 交换！
            final_mat0 = row['imp_chd_col_from_mat1']
            final_mat1 = row['imp_chd_col_from_mat0']
        else:
            # 既不是T也不是F，保持不变
            final_mat0 = row['imp_chd_col_from_mat0']
            final_mat1 = row['imp_chd_col_from_mat1']
        
        return pd.Series([final_mat0, final_mat1, final_pat])

    # 新增三列
    summary_df[['final_mat0', 'final_mat1', 'final_pat']] = summary_df.apply(
        assign_final_cols, axis=1
    )

    summary_df.to_csv(args.out_stat, sep="\t", index=False, na_rep="NA")

if __name__ == "__main__":
    main()
