#!/bin/bash
# README.sh - T21家系phase流程运行说明
# 运行前请确认：
# 1. 所有脚本中【需要修改的参数】（FAMILY_ID/样本名/BASE_DIR等）已修改
# 2. whatshap/shapeit4/gatk/bcftools/samtools等工具已安装并配置环境变量
# 3. 参考基因组/遗传图谱/人群面板/python脚本等文件路径修改确认

echo "========================================================="
echo "                       家系phase流程                     "
echo "                运行顺序：01 → 02 → 03 → 04              "
echo "========================================================="

# 步骤1：数据预处理 + 混合倍体jointcall
echo -e "\n==================== 运行 01.data_pre_and_jointcall.sh ===================="
bash 01.data_pre_and_jointcall.sh
if [ $? -ne 0 ]; then
    echo "ERROR: 01.data_pre_and_jointcall.sh 运行失败！"
    exit 1
fi

# 步骤2：家系phase + strict流程
echo -e "\n==================== 运行 02.phase.sh ===================="
bash 02.phase.sh
if [ $? -ne 0 ]; then
    echo "ERROR: 02.phase.sh 运行失败！"
    exit 1
fi

# 步骤3：ShapeIt4 母亲单倍型定型
echo -e "\n==================== 运行 03.shapeit_for_mat.sh ===================="
bash 03.shapeit_for_mat.sh
if [ $? -ne 0 ]; then
    echo "ERROR: 03.shapeit_for_mat.sh 运行失败！"
    exit 1
fi

# 步骤4：人工phase校正（hand phasing）
echo -e "\n==================== 运行 04.hand_phasing.sh ===================="
bash 04.hand_phasing.sh
if [ $? -ne 0 ]; then
    echo "ERROR: 04.hand_phasing.sh 运行失败！"
    exit 1
fi

echo -e "\n========================================================="
echo "                    所有流程运行完成！                     "
echo "========================================================="
