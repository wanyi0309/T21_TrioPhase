#!/bin/bash
set -eo pipefail

###########################################################
# 综合脚本：
# 00.data 提取chr21 BAM/VCF + 01.jointcall 子代三倍体gVCF + jointcall
###########################################################

##############################
# 修改输入参数
##############################
BASE_DIR="/data/work/06.chr21_phase/new_T21E6"
FAMILY_ID="T21E6"
CHROM="chr21"

# 原始输入文件路径
## 子代数据：二代短读WGS比对结果 + 三代长读WGS比对结果
CHD_SHORT_READS_BAM="/Files/chenjin/WGS_JIEGUO_BQSR/SC_T21E6_Midbrain_2_1/SC_T21E6_Midbrain_2.bqsr.bam"
CHD_LONG_READS_BAM="/Files/chenjin/cyclone_jieguo/SC_T21E6_Midbrain_2/SC_T21E6_Midbrain_2.merged.sorted.bam"

## 母本数据：三代长读WGS比对结果 + 二代WGS variant call的结果g.vcf
MAT_LONG_READS_BAM="/Files/chenjin/cyclone_jieguo/T21E6_maternal_blood/T21E6_maternal_blood.merged.sorted.bam"
MAT_GVCF_INPUT="/Files/chenjin/WGS_JIEGUO_BQSR/T21E6_maternal_blood_2/T21E6_maternal_blood_2.g.vcf.gz"

## 父本数据：三代长读WGS比对结果 + 二代WGS variant call的结果g.vcf
PAT_LONG_READS_BAM="/Files/chenjin/WGS_JIEGUO_BQSR/T21E6_paternal_blood_2/T21E6_paternal_blood_2.sort.bam"
PAT_GVCF_INPUT="/Files/chenjin/WGS_JIEGUO_BQSR/T21E6_paternal_blood_2/T21E6_paternal_blood_2.g.vcf.gz"

# 参考基因组
REF="/Files/chenjin/W202504210003221/reference/WGS_hg38.fa"

# 工具路径
SAMTOOLS="/share/app/samtools/1.11/bin/samtools"
BCFTOOLS="/share/app/bcftools/1.11/bin/bcftools"
GATK="gatk"
BGZIP="/share/app/htslib/bin/bgzip"
TABIX="/share/app/htslib/bin/tabix"

##############################
# 目录定义
##############################
DATA_DIR="${BASE_DIR}/00.data"
JOINTCALL_DIR="${BASE_DIR}/01.jointcall"
mkdir -p ${DATA_DIR} ${JOINTCALL_DIR}
cd ${BASE_DIR}

echo "====================================="
echo " 开始运行：${FAMILY_ID}"
echo " 工作目录：${BASE_DIR}"
echo "====================================="

###########################################################################
# 第一步：00.data 目录 —— 提取chr21 BAM + 父母gVCF
###########################################################################
echo -e "\n========== [1/3] 00.data：提取 chr21 BAM & GVCF =========="
cd ${DATA_DIR}

# 提取chr21 BAM
${SAMTOOLS} view -b ${MAT_LONG_READS_BAM} ${CHROM} > mat_${CHROM}_longreads.bam
${SAMTOOLS} view -b ${PAT_LONG_READS_BAM} ${CHROM} > pat_${CHROM}_longreads.bam
${SAMTOOLS} view -b ${CHD_LONG_READS_BAM} ${CHROM} > chd_${CHROM}_longreads.bam

# 建索引
${SAMTOOLS} index mat_${CHROM}_longreads.bam
${SAMTOOLS} index pat_${CHROM}_longreads.bam
${SAMTOOLS} index chd_${CHROM}_longreads.bam

# 提取父母 chr21 gVCF
${BCFTOOLS} view -r ${CHROM} ${MAT_GVCF_INPUT} -Oz -o mat_${CHROM}.g.vcf.gz
${BCFTOOLS} view -r ${CHROM} ${PAT_GVCF_INPUT} -Oz -o pat_${CHROM}.g.vcf.gz

# 建索引
${BCFTOOLS} index mat_${CHROM}.g.vcf.gz --tbi
${BCFTOOLS} index pat_${CHROM}.g.vcf.gz --tbi

###########################################################################
# 第二步：01.jointcall —— 子代三倍体 gVCF（ploidy 3）完全保留你原来脚本
###########################################################################
echo -e "\n========== [2/3] 01.jointcall：子代三倍体 HaplotypeCaller =========="
cd ${JOINTCALL_DIR}

TMP_DIR="${JOINTCALL_DIR}/tmp"
CHD_OUT_GVCF="${JOINTCALL_DIR}/chd_${CHROM}.g.vcf.gz"
mkdir -p ${TMP_DIR}

# 你原来的 gatk 命令一字不动
${GATK} --java-options "-Xms16g -Xmx48g" HaplotypeCaller \
    -R "${REF}" \
    -I "${CHD_SHORT_READS_BAM}" \
    -O "${CHD_OUT_GVCF}" \
    -L "${CHROM}" \
    -ERC GVCF \
    -ploidy 3 \
    --native-pair-hmm-threads 8 \
    --tmp-dir "${TMP_DIR}"

rm -rf ${TMP_DIR}
echo "子代三倍体gVCF完成：${CHD_OUT_GVCF}"

###########################################################################
# 第三步：jointcall
###########################################################################
echo -e "\n========== [3/3] 混合倍体 jointcall 流程 =========="

############################################
# 作用：
# 只针对 chr21 做 joint-calling
# 母本二倍体，父本二倍体，子代三倍体
#
# 输入文件必须分别是：
#   母本 chr21 的二倍体 gVCF
#   父本 chr21 的二倍体 gVCF
#   子代 chr21 的三倍体 gVCF
############################################

# 输入文件
MOTHER_GVCF="${DATA_DIR}/mat_${CHROM}.g.vcf.gz"
FATHER_GVCF="${DATA_DIR}/pat_${CHROM}.g.vcf.gz"
CHD_GVCF="${CHD_OUT_GVCF}"

# 输出目录
OUTPUT_DIR="${JOINTCALL_DIR}/${CHROM}_mixed_ploidy"
GATK_DB_BASE="${JOINTCALL_DIR}/01.gatk_db"
GATK_DB="${GATK_DB_BASE}/${CHROM}_mixed_ploidy"

# 过滤参数
FILTER_DP_MIN=20
FILTER_DP_MAX=400

# Java 和线程参数
GENOMICSDB_THREADS=16
GENOTYPE_JAVA_XMS="8g"
GENOTYPE_JAVA_XMX="24g"

# 临时目录
TMP_DIR_JOINT="${OUTPUT_DIR}/tmp_${CHROM}"

mkdir -p ${OUTPUT_DIR} ${TMP_DIR_JOINT} ${GATK_DB_BASE}

# 输入检查
for f in ${MOTHER_GVCF} ${FATHER_GVCF} ${CHD_GVCF}; do
    [ ! -s ${f} ] && echo "ERROR: 缺失文件 ${f}" && exit 1
done

# 1. GenomicsDBImport
# 如果数据库目录已存在，GATK 可能报错。
# 需要重跑时，先手动删除：
# rm -r "${GATK_DB}"
${GATK} GenomicsDBImport \
    -V ${MOTHER_GVCF} \
    -V ${FATHER_GVCF} \
    -V ${CHD_GVCF} \
    --reader-threads ${GENOMICSDB_THREADS} \
    --genomicsdb-workspace-path ${GATK_DB} \
    -L ${CHROM}

# 2. GenotypeGVCFs
${GATK} --java-options "-Xms${GENOTYPE_JAVA_XMS} -Xmx${GENOTYPE_JAVA_XMX}" GenotypeGVCFs \
    -R ${REF} \
    -V "gendb://${GATK_DB}" \
    -O ${OUTPUT_DIR}/${FAMILY_ID}.${CHROM}.joint.vcf.gz \
    -L ${CHROM} \
    --tmp-dir ${TMP_DIR_JOINT}

# 3. 提取SNP
${GATK} SelectVariants \
    -V ${OUTPUT_DIR}/${FAMILY_ID}.${CHROM}.joint.vcf.gz \
    -select-type SNP \
    -O ${OUTPUT_DIR}/${FAMILY_ID}.${CHROM}.joint.snps.vcf.gz

# 4. 硬过滤
${GATK} VariantFiltration \
    -V ${OUTPUT_DIR}/${FAMILY_ID}.${CHROM}.joint.snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -filter "DP < ${FILTER_DP_MIN}" --filter-name "DP_low_${FILTER_DP_MIN}" \
    -filter "DP > ${FILTER_DP_MAX}" --filter-name "DP_high_${FILTER_DP_MAX}" \
    -O ${OUTPUT_DIR}/${FAMILY_ID}.${CHROM}.joint.snps.filtered.vcf.gz

# 5. 只保留PASS
zcat ${OUTPUT_DIR}/${FAMILY_ID}.${CHROM}.joint.snps.filtered.vcf.gz | \
awk '/^#/ || $7=="PASS"' | \
${BGZIP} -c > ${OUTPUT_DIR}/${FAMILY_ID}.${CHROM}.joint.snps.filtered.pass.vcf.gz

# 6. 索引
${TABIX} -p vcf ${OUTPUT_DIR}/${FAMILY_ID}.${CHROM}.joint.snps.filtered.pass.vcf.gz

rm -rf ${TMP_DIR_JOINT}

echo -e "\n====================================="
echo " 运行完成！"
echo " 最终SNP文件：${OUTPUT_DIR}/${FAMILY_ID}.${CHROM}.joint.snps.filtered.pass.vcf.gz"
echo "====================================="