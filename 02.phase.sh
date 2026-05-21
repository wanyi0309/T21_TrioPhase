#!/bin/bash
#set -euo pipefail
echo "========================================================="
echo "           一体化家系phase + strict流程 自动脚本          "
echo "                  所有参数都在脚本顶部                    "
echo "========================================================="

# ==============================================
# 【唯一需要你修改的地方：所有参数在这里】
# ==============================================
# 基本信息
FAMILY_ID="T21E6"               # 家系编号 T21E5 / T21E11
CHR="chr21"                     # 染色体
CHILD_SAMPLE="SC_T21E6_Midbrain_2"       # 子代样本名 ！！！与gvcf中sample的ID一致
MATERNAL_SAMPLE="T21E6_maternal_blood_2" # 母亲样本名 ！！！与gvcf中sample的ID一致
PATERNAL_SAMPLE="T21E6_paternal_blood_2" # 父亲样本名 ！！！与gvcf中sample的ID一致

# python脚本目录
SCRIPT_DIR="/data/work/06.chr21_phase/new_T21E6/all_script/"

# 参考基因组
REF_FA="/Files/chenjin/W202504210003221/reference/WGS_hg38.fa"

# ======================
# 和上一步 jointcall 脚本对齐
BASE_DIR="/data/work/06.chr21_phase/new_T21E6" 

# 输入 joint VCF（来自上一步01.data_pre_and_jointcall.sh输出）
INPUT_JOINT_VCF="${BASE_DIR}/01.jointcall/chr21_mixed_ploidy/${FAMILY_ID}.${CHR}.joint.snps.filtered.pass.vcf.gz"

# BAM 文件（和上一步 00.data 输出完全匹配）
MAT_BAM="${BASE_DIR}/00.data/mat_${CHR}_longreads.bam"
PAT_BAM="${BASE_DIR}/00.data/pat_${CHR}_longreads.bam"
CHILD_BAM="${BASE_DIR}/00.data/chd_${CHR}_longreads.bam"

# 线程
THREADS_DEFAULT=8
THREADS_STRICT=8

# block 过滤脚本（你自己路径）
SCRIPT_BLOCK_FILTER="${SCRIPT_DIR}/py_script/001.strict_block_filter_for_default_switcherror.py"

# bcftools / bgzip 路径
BCFTOOLS="/share/app/bcftools/1.11/bin/bcftools"
BGZIP="/share/app/htslib/bin/bgzip"

# work dir（自动创建）
WORK_DIR="${BASE_DIR}/02.phase"


# ==============================================
# 以下参数无需修改 自动运行流程
# ==============================================
PREFIX="${FAMILY_ID}.${CHR}"
SNPS_ONLY_VCF="${PREFIX}.snponly.vcf.gz"
SAMPLE_LIST="${SNPS_ONLY_VCF}.sample.list"

cd ${WORK_DIR}
mkdir -p strict logs

# ========== 1. 筛选双等位SNP ==========
echo -e "\n[1-1/7] 筛选只有两个allele的SNPs ..."
${BCFTOOLS} view -m2 -M2 -v snps -Oz -o ${SNPS_ONLY_VCF}.tmp.vcf.gz ${INPUT_JOINT_VCF}
${BCFTOOLS} index ${SNPS_ONLY_VCF}.tmp.vcf.gz  --tbi

# ========== 【关键新增】过滤：只保留父母子三个人都有GT（无缺失）的位点 ==========
echo -e "\n[1-2/7] 过滤掉任何基因型缺失的位点（不含 ./. 就保留） 即保证父母子在保留位点均有GT信息"
${BCFTOOLS} view ${SNPS_ONLY_VCF}.tmp.vcf.gz \
| grep -v "\./\." \
| ${BGZIP} -c > ${SNPS_ONLY_VCF}
${BCFTOOLS} index ${SNPS_ONLY_VCF} 
rm -f ${SNPS_ONLY_VCF}.tmp.vcf.gz ${SNPS_ONLY_VCF}.tmp.vcf.gz.tbi

# ========== 2. 拆分样本 ==========
echo -e "\n[2/7] 提取样本列表并拆分父母子 ..."
${BCFTOOLS} query -l ${SNPS_ONLY_VCF} > ${SAMPLE_LIST}
cat ${SAMPLE_LIST} | while read sample; do
  ${BCFTOOLS} view -s ${sample} ${SNPS_ONLY_VCF} -Oz -o ${SNPS_ONLY_VCF%.vcf.gz}.${sample}.vcf.gz
  ${BCFTOOLS} index ${SNPS_ONLY_VCF%.vcf.gz}.${sample}.vcf.gz
done

# ========== 3. phase 父母 + 子代polyphase ==========
echo -e "\n[3/7] Whatshap phase 父母 + 子代 ..."
whatshap phase -o ${PREFIX}.${MATERNAL_SAMPLE}.phased.vcf.gz \
  --reference ${REF_FA} ${SNPS_ONLY_VCF%.vcf.gz}.${MATERNAL_SAMPLE}.vcf.gz ${MAT_BAM} \
  --ignore-read-groups > logs/${MATERNAL_SAMPLE}.phase.log 2>&1
echo "mat phase done!!!"

whatshap phase -o ${PREFIX}.${PATERNAL_SAMPLE}.phased.vcf.gz \
  --reference ${REF_FA} ${SNPS_ONLY_VCF%.vcf.gz}.${PATERNAL_SAMPLE}.vcf.gz ${PAT_BAM} \
  --ignore-read-groups > logs/${PATERNAL_SAMPLE}.phase.log 2>&1
echo "pat phase done!!!"

whatshap polyphase ${SNPS_ONLY_VCF%.vcf.gz}.${CHILD_SAMPLE}.vcf.gz ${CHILD_BAM} \
  --ploidy 3 --ignore-read-groups --reference ${REF_FA} \
  -o ${PREFIX}.${CHILD_SAMPLE}.phased.vcf.gz -t ${THREADS_DEFAULT} \
  > logs/${CHILD_SAMPLE}.phase.log 2>&1
echo "chd phase done!!!"

# ========== 4. 清理VCF ==========
echo -e "\n[4/7] 清理VCF并建立索引 ..."
cat ${SAMPLE_LIST} | while read sample; do
  ${BCFTOOLS} annotate -x FORMAT/PL -Oz -o ${PREFIX}.${sample}.phased.clean.vcf.gz ${PREFIX}.${sample}.phased.vcf.gz
  ${BCFTOOLS} index ${PREFIX}.${sample}.phased.clean.vcf.gz
done

# ========== 5. 合并全家系 + 过滤 ==========
echo -e "\n[5/7] 合并全家系VCF ..."
${BCFTOOLS} merge -m all -Oz -o ${FAMILY_ID}.family.phased.vcf.gz \
  ${PREFIX}.${CHILD_SAMPLE}.phased.clean.vcf.gz \
  ${PREFIX}.${MATERNAL_SAMPLE}.phased.clean.vcf.gz \
  ${PREFIX}.${PATERNAL_SAMPLE}.phased.clean.vcf.gz

${BCFTOOLS} filter -i 'FMT/GT[0]!="0/0/0" && FMT/GT[0]!="1/1/1"' \
  -Oz -o ${FAMILY_ID}.family.phased.son_het.vcf.gz ${FAMILY_ID}.family.phased.vcf.gz
${BCFTOOLS} index ${FAMILY_ID}.family.phased.son_het.vcf.gz

# ========== 6. block 统计 ==========
echo -e "\n[6/7] 生成phase and filter block统计 ..."
cat ${SAMPLE_LIST} | while read sample; do
  ${BCFTOOLS} view -s ${sample} ${FAMILY_ID}.family.phased.son_het.vcf.gz -Oz -o tmp.${sample}.vcf.gz
  whatshap stats tmp.${sample}.vcf.gz --block-list ${FAMILY_ID}.family.phased.son_het.${sample}.vcf.gz.block.tsv > logs/block.stats.${sample}.txt 2>&1
  #rm -f tmp.${sample}.vcf.gz
done

# ========== 7. strict 流程 ==========
echo -e "\n[7/7] 运行 strict polyphase 流程 ..."
cd strict

whatshap polyphase \
  --ploidy 3 --threads ${THREADS_STRICT} \
  --block-cut-sensitivity 5 --min-overlap 3 --no-mav \
  --ignore-read-groups \
  -o ${PREFIX}.${CHILD_SAMPLE}.phased.strict.vcf.gz \
  --reference ${REF_FA} \
  ../${SNPS_ONLY_VCF%.vcf.gz}.${CHILD_SAMPLE}.vcf.gz ${CHILD_BAM} \
  > ../logs/${CHILD_SAMPLE}.strict.phase.log 2>&1

${BCFTOOLS} view -Oz -o ${PREFIX}.${CHILD_SAMPLE}.phased.strict.filter.vcf.gz \
  -T ../${FAMILY_ID}.family.phased.son_het.vcf.gz \
  ${PREFIX}.${CHILD_SAMPLE}.phased.strict.vcf.gz
${BCFTOOLS} index ${PREFIX}.${CHILD_SAMPLE}.phased.strict.filter.vcf.gz

whatshap stats ${PREFIX}.${CHILD_SAMPLE}.phased.strict.filter.vcf.gz \
  --block-list ${PREFIX}.${CHILD_SAMPLE}.phased.strict.filter.vcf.gz.block.tsv

${BCFTOOLS} annotate -x FORMAT/PL -Oz -o ${PREFIX}.${CHILD_SAMPLE}.phased.strict.filter.clean.vcf.gz \
  ${PREFIX}.${CHILD_SAMPLE}.phased.strict.filter.vcf.gz
${BCFTOOLS} index ${PREFIX}.${CHILD_SAMPLE}.phased.strict.filter.clean.vcf.gz

${BCFTOOLS} merge -m all -Oz -o ${FAMILY_ID}.family.phased.strict_chd.vcf.gz \
  ${PREFIX}.${CHILD_SAMPLE}.phased.strict.filter.clean.vcf.gz \
  ../${PREFIX}.${MATERNAL_SAMPLE}.phased.clean.vcf.gz \
  ../${PREFIX}.${PATERNAL_SAMPLE}.phased.clean.vcf.gz

python ${SCRIPT_BLOCK_FILTER} \
  ${FAMILY_ID}.family.phased.strict_chd.vcf.gz \
  ${CHILD_SAMPLE} \
  ${MATERNAL_SAMPLE} \
  ${PATERNAL_SAMPLE} \
  ${PREFIX}.${CHILD_SAMPLE}.phased.strict.filter.vcf.gz.block.tsv \
  ../${FAMILY_ID}.family.phased.son_het.${MATERNAL_SAMPLE}.vcf.gz.block.tsv \
  ../${FAMILY_ID}.family.phased.son_het.${PATERNAL_SAMPLE}.vcf.gz.block.tsv

cd ..

echo -e "\n========================================================="
echo "                      全部运行完成！                     "
echo "========================================================="