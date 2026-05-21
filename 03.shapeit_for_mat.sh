#!/bin/bash
set -eo pipefail

echo "========================================================="
echo "     ShapeIt4 母亲单倍型依据人群信息进行定型 + 区块判断   "
echo "                  所有参数都在脚本顶部                    "
echo "========================================================="

# ==============================================
# 【唯一需要修改：样本参数】
# 与前面 jointcall / phase 脚本完全一致
# ==============================================
FAMILY_ID="T21E6"
CHR="chr21"
MATERNAL_SAMPLE="T21E6_maternal_blood_2"  # 母亲样本ID

# 主目录（和前面脚本完全一致）
BASE_DIR="/data/work/06.chr21_phase/new_T21E6"

# 工具 & 数据目录（你提供的路径）
SCRIPT_DIR="${BASE_DIR}/all_script"
SHAPEIT_DATA_DIR="${SCRIPT_DIR}/shapeit_data"
PY_SCRIPT_DIR="${SCRIPT_DIR}/py_script"

# 输入（来自上一步 02.phase）
PHASE_DIR="${BASE_DIR}/02.phase"
MAT_PHASED_VCF="${PHASE_DIR}/${FAMILY_ID}.${CHR}.${MATERNAL_SAMPLE}.phased.clean.vcf.gz" ##这个vcf有更多的母本SNP信息，给shapeit更多vcf信息
MAT_BLOCK_TSV="${PHASE_DIR}/${FAMILY_ID}.family.phased.son_het.${MATERNAL_SAMPLE}.vcf.gz.block.tsv" ##这个vcf是后续手动phase参考用的vcf
MAT_FAMILY_VCF="${PHASE_DIR}/${FAMILY_ID}.family.phased.son_het.${MATERNAL_SAMPLE}.vcf.gz" ##这个vcf是后续手动phase参考用的vcf

# 输出目录
OUT_DIR="${BASE_DIR}/03.shapeit_for_mat"
mkdir -p ${OUT_DIR}
cd ${OUT_DIR}

# shapeit4 依赖文件
GMAP="${SHAPEIT_DATA_DIR}/${CHR}.b38.gmap.fixed.gz"
REFERENCE_PANEL="${SHAPEIT_DATA_DIR}/EAS.1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

# Python 脚本
PY_SCRIPT="${PY_SCRIPT_DIR}/002.mat_block_link_relationship.py"

# 线程
THREADS=4

# ==============================================
# 流程 1：运行 shapeit4
# ==============================================
echo -e "\n[1/2] 运行 ShapeIt4 母亲单倍型定型..."
SHAPEIT_OUTPUT="${MATERNAL_SAMPLE}.shapeit4.EASpopu.PS0.vcf.gz"

shapeit4 \
  --input ${MAT_PHASED_VCF} \
  --map ${GMAP} \
  --region ${CHR} \
  --output ${SHAPEIT_OUTPUT} \
  --thread ${THREADS} \
  --use-PS 0 \
  --reference ${REFERENCE_PANEL} \
  --sequencing

echo "ShapeIt4 完成！输出：${SHAPEIT_OUTPUT}"

# ==============================================
# 流程 2：运行 Python 脚本，生成 block 连锁关系
# ==============================================
echo -e "\n[2/2] 生成母亲 block 连锁关系文件..."

python ${PY_SCRIPT} \
  --ori_block ${MAT_BLOCK_TSV} \
  --ori_vcf ${MAT_FAMILY_VCF} \
  --shapeit_vcf ${SHAPEIT_OUTPUT} \
  --out mat_block_link_relationship.tsv

echo -e "\n========================================================="
echo "                    全部运行完成！"
echo " 输出文件：${OUT_DIR}/mat_block_link_relationship.tsv"
echo "========================================================="