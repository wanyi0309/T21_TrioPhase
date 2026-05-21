import pandas as pd
import numpy as np
import argparse

def process_block_stats(df):
    # ===================== 分组统计 =====================
    ori = df.groupby("chd_block")["type"].value_counts().unstack(fill_value=0)
    valid_types = ["A1", "A3", "B2", "C", "F"]
    for col in valid_types:
        if col not in ori.columns:
            ori[col] = 0

    # ===================== 自然排序 (sort -k1,1V) =====================
    ori = ori.reset_index()
    ori[['chr', 'coords']] = ori['chd_block'].str.split(':', expand=True)
    ori['start'] = ori['coords'].str.split('-').str[0].astype(int)
    ori = ori.sort_values(['chr', 'start'])
    ori = ori.drop(['chr', 'coords', 'start'], axis=1).set_index('chd_block')

    # ===================== 添加计算列 =====================
    ori["all_total"] = ori[valid_types].sum(axis=1)
    ori["valid_chd_block"] = ori["all_total"].apply(lambda x: "T" if x > 1 else "F")
    ori["total_ABC"] = ori[["A1", "A3", "B2", "C"]].sum(axis=1)
    ori["total_BC"] = ori[["B2", "C"]].sum(axis=1)

    # 三列判断
    ori["A1>0"] = ori["A1"].apply(lambda x: "T" if x > 0 else "F")
    ori["A1>=5"] = ori["A1"].apply(lambda x: "T" if x >= 5 else "F")
    ori["A1>=0&B2+C>=0"] = ori.apply(lambda r: "T" if r["A1"]>=0 and r["total_BC"]>=0 else "F", axis=1)
    ori["A1>=5&B2+C>=5"] = ori.apply(lambda r: "T" if r["A1"]>=5 and r["total_BC"]>=5 else "F", axis=1)

    # ===================== invalid block 置 NA =====================
    F_mask = ori["valid_chd_block"] == "F"
    cols_to_null = [
        "total_ABC", "total_BC",
        "A1>0", "A1>=5", "A1>=0&B2+C>=0", "A1>=5&B2+C>=5"
    ]
    ori.loc[F_mask, cols_to_null] = pd.NA

    return ori

def main():
    # ===================== 1. 命令行参数 =====================
    parser = argparse.ArgumentParser(description="origin chd_block 统计与质控脚本")
    parser.add_argument("--input", "-i", default="annotated.tsv", help="输入文件路径")
    args = parser.parse_args()

    # ===================== 2. 读取数据 =====================
    df = pd.read_csv(
        args.input,
        sep="\t",
        na_values="NA",
        dtype={"pos": int}
    )

    # ===================== 3. 生成两份表格 =====================
    # 全量数据
    print("处理全量数据...")
    all_stats = process_block_stats(df)
    all_stats.to_csv("origin_type_stat_all.tsv", sep="\t", index=True, na_rep="NA")

    # pos > 13000000
    print("处理 pos > 13000000 数据...")
    df_up13M = df[df["pos"] > 13000000].copy()
    up13M_stats = process_block_stats(df_up13M)
    up13M_stats.to_csv("origin_type_stat_up13M.tsv", sep="\t", index=True, na_rep="NA")

    print("\n全部完成！")
    print("全量结果：origin_type_stat_all.csv")
    print("pos>13M 结果：origin_type_stat_up13M.csv")

if __name__ == "__main__":
    main()