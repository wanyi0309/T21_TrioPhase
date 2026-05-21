import pandas as pd
import numpy as np
import argparse

# ===================== 分析函数（复用逻辑） =====================
def analyze_block_stats(df):
    # 【核心】按 step1_newchd_block 分组统计 type
    stat = df[df["valid_chd_block"]=="T"].groupby("step1_newchd_block")["type"].value_counts().unstack(fill_value=0)

    # 自然排序（sort -k1,1V）
    stat = stat.reset_index()
    stat[['chr', 'coords']] = stat['step1_newchd_block'].str.split(':', expand=True)
    stat['start'] = stat['coords'].str.split('-').str[0].astype(int)
    stat = stat.sort_values(['chr', 'start'])
    stat = stat.drop(['chr', 'coords', 'start'], axis=1).set_index('step1_newchd_block')

    # ===================== 添加计算列 =====================
    valid_types = ["A1", "A3", "B2", "C", "F"]
    # 确保列存在
    for col in valid_types:
        if col not in stat.columns:
            stat[col] = 0

    stat["all_total"] = stat[valid_types].sum(axis=1)
    stat["total_ABC"] = stat[["A1", "A3", "B2", "C"]].sum(axis=1)
    stat["valid_block"] = stat.apply(lambda row: "F" if row["A1"] == row["total_ABC"] else "T", axis=1)
    stat["total_BC"] = stat[["B2", "C"]].sum(axis=1)
    stat["A1>=5"] = stat["A1"].apply(lambda x: "T" if x >= 5 else "F")
    stat["A1>=5&B2+C>=5"] = stat.apply(lambda r: "T" if r["A1"]>=5 and r["total_BC"]>=5 else "F", axis=1)

    # ===================== pat 映射 =====================
    pat_mapping = df[df["valid_chd_block"]=="T"].groupby("step1_newchd_block")["step1_newchd_block_pat"].apply(
        lambda x: x.dropna().unique()
    )
    # 检查一个block是否对应多个pat
    bad_blocks = [block for block, vals in pat_mapping.items() if len(vals) > 1]
    if bad_blocks:
        print("❌ 发现block对应多个pat：", bad_blocks)
    pat_mapping = pat_mapping.str[0]

    # 加入pat列
    stat["pat"] = pat_mapping

    return stat

# ===================== 主函数 =====================
def main():
    parser = argparse.ArgumentParser(description="after step1 step1_newchd_block 统计脚本")
    parser.add_argument("--input", "-i", required=True, help="输入文件路径 annotated.tsv")
    args = parser.parse_args()

    # 读取数据
    df = pd.read_csv(
        args.input,
        sep="\t",
        na_values="NA",
        dtype={"pos": int}
    )

    # ===================== 1. 全部 SNP 分析 =====================
    print("🔹 处理全部数据...")
    stat_all = analyze_block_stats(df)
    stat_all.to_csv("afterstep1_block_stat_all.tsv", sep="\t", na_rep="NA")

    # ===================== 2. pos > 13000000 分析 =====================
    print("🔹 处理 pos > 13000000 数据...")
    df_up13M = df[df["pos"] > 13000000].copy()
    stat_up13M = analyze_block_stats(df_up13M)
    stat_up13M.to_csv("afterstep1_block_stat_up13M.tsv", sep="\t", na_rep="NA")

    print("\n✅ 分析完成！")
    print("📄 全部数据：afterstep1_block_stat_all.tsv")
    print("📄 pos>13M 数据：afterstep1_block_stat_up13M.tsv")

if __name__ == "__main__":
    main()