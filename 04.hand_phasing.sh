#!/bin/bash
set -eo pipefail

echo "========================================================="
echo "             04. 人工pahse流程 (hand phasing)            "
echo "========================================================="

# ==============================================
# 【唯一需要改：样本参数，和前面完全一致】
# ==============================================
FAMILY_ID="T21E6"
CHR="chr21"
CHILD_SAMPLE="SC_T21E6_Midbrain_2"
MATERNAL_SAMPLE="T21E6_maternal_blood_2"
PATERNAL_SAMPLE="T21E6_paternal_blood_2"

# 主目录（和前面全套脚本一致）
BASE_DIR="/data/work/new_T21E6"

# 脚本目录（统一存放的位置）
SCRIPT_DIR="${BASE_DIR}/all_script"
PY_SCRIPT_DIR="${SCRIPT_DIR}/py_script"

# 前面流程的输出目录（自动读取）
PHASE_DIR="${BASE_DIR}/02.phase"
SHAPEIT_DIR="${BASE_DIR}/03.shapeit_for_mat"

# 本流程主目录
HAND_PHASE_DIR="${BASE_DIR}/04.hand_phasing"
STEP1_DIR="${HAND_PHASE_DIR}/01.first_step_clarify_pat"
STEP2_DIR="${HAND_PHASE_DIR}/02.second_step_clarify_mat"
STEP3_DIR="${HAND_PHASE_DIR}/03.expand_SNP_and_merge_block"

# 自动创建全部目录
mkdir -p ${STEP1_DIR} ${STEP2_DIR} ${STEP3_DIR}

# ==============================
# 输入文件（自动从上一步读取）
# ==============================
FAMILY_VCF="${PHASE_DIR}/${FAMILY_ID}.family.phased.son_het.vcf.gz"
CHILD_BLOCK_TSV="${PHASE_DIR}/${FAMILY_ID}.family.phased.son_het.${CHILD_SAMPLE}.vcf.gz.block.tsv"
MAT_BLOCK_TSV="${PHASE_DIR}/${FAMILY_ID}.family.phased.son_het.${MATERNAL_SAMPLE}.vcf.gz.block.tsv"
PAT_BLOCK_TSV="${PHASE_DIR}/${FAMILY_ID}.family.phased.son_het.${PATERNAL_SAMPLE}.vcf.gz.block.tsv"
STRICT_BLOCK_TSV="${PHASE_DIR}/strict/strick_blocks_for_default_switch_check.tsv"
MAT_LINK_FILE="${SHAPEIT_DIR}/mat_block_link_relationship.tsv"

# ==============================================
# 【STEP1】明确父型单倍型
# ==============================================
echo -e "\n========== [1/3] 第一步：明确父型单倍型 =========="
cd ${STEP1_DIR}

python ${PY_SCRIPT_DIR}/003.vcf_gt_clerify_pat_hap.py \
    ${FAMILY_VCF} \
    ${CHILD_SAMPLE} \
    ${MATERNAL_SAMPLE} \
    ${PATERNAL_SAMPLE} \
    ${CHILD_BLOCK_TSV} \
    ${MAT_BLOCK_TSV} \
    ${PAT_BLOCK_TSV} \
    ${STRICT_BLOCK_TSV}

# 统计
python ${PY_SCRIPT_DIR}/003_1.ori_SNP.stat.py --input annotated.tsv
python ${PY_SCRIPT_DIR}/003_2.afterstep1_SNP.stat.py --input annotated.tsv

# ==============================================
# 【STEP2】明确母型单倍型
# ==============================================
echo -e "\n========== [2/3] 第二步：明确母型单倍型 =========="
cd ${STEP2_DIR}

# 软链接上一步结果
cp ../01.first_step_clarify_pat/annotated.tsv ./

# 运行
python ${PY_SCRIPT_DIR}/004.vcf_gt_clerify_mat_hap.fix.v2.py

# ==============================================
# 【STEP3】SNP扩展 + 合并block + 生成最终VCF
# ==============================================
echo -e "\n========== [3/3] 第三步：SNP扩展 & 合并block & 生成最终vcf =========="
cd ${STEP3_DIR}

# 软链接上一步结果
cp ../02.second_step_clarify_mat/annotated.step2.tsv ./
cp -sf ../02.second_step_clarify_mat/annotated.step2.block.tsv ./

# 1. 合并预处理
python ${PY_SCRIPT_DIR}/005.merge_pre.py \
    --infile annotated.step2.tsv \
    --mat_link_file ${MAT_LINK_FILE} \
    --chr_name ${CHR} \
    --out_snp filterSNP.add_A1A3hapinfo.txt \
    --out_stat step3_filter.stat



# 2. 生成最终子代单倍型VCF
python ${PY_SCRIPT_DIR}/006.make_new_chd_vcf.py \
    --keep filterSNP.add_A1A3hapinfo.txt \
    --stat step3_filter.stat \
    --vcf ${FAMILY_VCF} \
    --child ${CHILD_SAMPLE} \
    --out final.child.hap.vcf

echo -e "\n========================================================="
echo "                04.hand_phasing 全部完成！   "
echo " 最终子代单倍型VCF：${STEP3_DIR}/final.child.hap.vcf"
echo "========================================================="
