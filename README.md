# T21_TrioPhase

A shell/Python workflow for **trio-based haplotype phasing of chromosome 21 in trisomy 21 (T21)**.

This pipeline is designed for a T21 trio in which the child carries three copies of chromosome 21 while both parents are diploid. It is especially useful for cases inferred to originate from **maternal meiosis I nondisjunction**, where the goal is to resolve the child’s chr21 haplotypes and distinguish paternal and maternal haplotype contributions.

## Overview

The workflow contains four major steps:

```text
01.data_pre_and_jointcall.sh
    ↓
02.phase.sh
    ↓
03.shapeit_for_mat.sh
    ↓
04.hand_phasing.sh
    ↓
final.child.hap.vcf
```

## Repository structure

```text
T21_TrioPhase/
├── 01.data_pre_and_jointcall.sh      # chr21 extraction, child ploidy-3 gVCF, mixed-ploidy joint calling
├── 02.phase.sh                       # WhatsHap-based trio phasing and strict child polyphase
├── 03.shapeit_for_mat.sh             # ShapeIt4-assisted maternal haplotype anchoring
├── 04.hand_phasing.sh                # Rule-based paternal/maternal haplotype clarification and final VCF generation
├── total.sh                         # One-command wrapper: run 01 → 02 → 03 → 04
├── py_script/                        # Helper Python scripts used by steps 02–04
└── shapeit_data/                     # Genetic map and population reference panel information for ShapeIt4
```

## Requirements

Required command-line tools:

- `samtools`
- `bcftools`
- `bgzip` and `tabix` from HTSlib
- `GATK4`
- `WhatsHap`
- `ShapeIt4`
- `Python 3`

Required Python packages:

- `pandas`

Input data should include:

- Child short-read WGS BAM for ploidy-3 GATK variant calling
- Child long-read BAM for polyploid phasing
- Maternal long-read BAM
- Paternal long-read BAM
- Maternal diploid gVCF
- Paternal diploid gVCF
- Reference genome FASTA, for example hg38, with matching index files
- ShapeIt4 chr21 genetic map
- ShapeIt4 population reference panel, preferably matched to the study ancestry

> Important: chromosome naming must be consistent across all files, for example `chr21` rather than `21`.

## Quick start

### 1. Prepare the scripts

Clone the repository and place all scripts in your working directory:

```bash
git clone https://github.com/wanyi0309/T21_TrioPhase.git
cd T21_TrioPhase
```

Edit the user-defined parameters at the top of each shell script before running:

```bash
FAMILY_ID="T21E6"
CHR="chr21"
CHILD_SAMPLE="SC_T21E6_Midbrain_2"
MATERNAL_SAMPLE="T21E6_maternal_blood_2"
PATERNAL_SAMPLE="T21E6_paternal_blood_2"
BASE_DIR="/path/to/your/project"
SCRIPT_DIR="/path/to/T21_TrioPhase"
REF="/path/to/reference.fa"
```

Also check the paths to `samtools`, `bcftools`, `GATK`, `bgzip`, `tabix`, `WhatsHap`, and `ShapeIt4`.

### 2. Run the full pipeline

After all paths are configured, the full workflow can be launched with:

```bash
bash total.sh
```

This wrapper runs the four main scripts in order:

```bash
bash 01.data_pre_and_jointcall.sh
bash 02.phase.sh
bash 03.shapeit_for_mat.sh
bash 04.hand_phasing.sh
```

You can also run each step separately for debugging or stepwise inspection.

## Workflow details

### Step 01: data preparation and mixed-ploidy joint calling

Script:

```bash
01.data_pre_and_jointcall.sh
```

Main tasks:

1. Extract chr21 reads from parental and child long-read BAM files.
2. Extract parental chr21 gVCFs.
3. Run GATK `HaplotypeCaller` for the child on chr21 with `-ploidy 3`.
4. Joint-call maternal diploid gVCF, paternal diploid gVCF, and child triploid gVCF.
5. Select SNPs, apply hard filters, and retain only PASS variants.

Main output:

```text
01.jointcall/chr21_mixed_ploidy/${FAMILY_ID}.chr21.joint.snps.filtered.pass.vcf.gz
```

This file is used as the input VCF for the phasing step.

### Step 02: trio phasing and strict child polyphase

Script:

```bash
02.phase.sh
```

Main tasks:

1. Keep biallelic SNPs only.
2. Remove sites with missing genotypes in any trio member.
3. Split the VCF by sample.
4. Phase the diploid parents using `whatshap phase`.
5. Phase the triploid child using `whatshap polyphase --ploidy 3`.
6. Merge the phased trio VCFs.
7. Keep informative child-heterozygous sites.
8. Generate phase-block statistics.
9. Run a stricter child polyphase setting and filter/split blocks with potential switch errors.

Representative outputs:

```text
02.phase/${FAMILY_ID}.family.phased.vcf.gz
02.phase/${FAMILY_ID}.family.phased.son_het.vcf.gz
02.phase/${FAMILY_ID}.family.phased.son_het.${CHILD_SAMPLE}.vcf.gz.block.tsv
02.phase/${FAMILY_ID}.family.phased.son_het.${MATERNAL_SAMPLE}.vcf.gz.block.tsv
02.phase/${FAMILY_ID}.family.phased.son_het.${PATERNAL_SAMPLE}.vcf.gz.block.tsv
02.phase/strict/${FAMILY_ID}.chr21.${CHILD_SAMPLE}.phased.strict.filter.vcf.gz.block.tsv
```

### Step 03: ShapeIt4-assisted maternal haplotype anchoring

Script:

```bash
03.shapeit_for_mat.sh
```

Main tasks:

1. Run `shapeit4` on the maternal phased VCF using a chr21 genetic map and a population reference panel.
2. Compare the original maternal phase blocks with the ShapeIt4-phased output.
3. Generate a block-level relationship table indicating whether each maternal block is consistent with the population-assisted phase.

Main output:

```text
03.shapeit_for_mat/mat_block_link_relationship.tsv
```

This file is used in the final hand-phasing step to anchor and merge maternal haplotype blocks.

## ShapeIt4 reference data

The directory `shapeit_data/` is used to store ShapeIt4 dependency files.

Expected files include:

```text
chr21.b38.gmap.fixed.gz
EAS.1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
EAS.1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.csi
```

Example preparation of the chr21 genetic map:

```bash
tar -zxvf genetic_maps.b38.tar.gz
zcat chr21.b38.gmap.gz \
  | awk 'BEGIN{FS="\t";OFS="\t"} $2==21{$2="chr21"} 1' \
  | gzip > chr21.b38.gmap.fixed.gz
```

Example extraction of East Asian samples from the 1000 Genomes high-coverage phased panel:

```bash
awk '$7=="EAS"' 20130606_g1k_3202_samples_ped_population.txt \
  | awk '{print $2}' > EAS_samples.txt

bcftools view \
  -S EAS_samples.txt \
  -m2 -M2 -v snps \
  -Oz \
  -o EAS.1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
  1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz

bcftools index EAS.1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
```

Because the population reference panel is large, it may need to be generated locally rather than stored directly in the repository.

### Step 04: rule-based hand phasing

Script:

```bash
04.hand_phasing.sh
```

Main tasks:

1. Clarify paternal haplotype contribution using informative genotype patterns and child phase blocks.
2. Generate SNP-level statistics before and after paternal clarification.
3. Clarify maternal haplotype contribution.
4. Expand informative SNPs and merge blocks.
5. Generate the final child haplotype VCF.

Main output:

```text
04.hand_phasing/03.expand_SNP_and_merge_block/final.child.hap.vcf ##pat|mat0|mat1
```

This is the final phased child chr21 VCF produced by the pipeline. 

## Helper Python scripts

The `py_script/` directory contains helper scripts used by the main shell workflow.

| Script | Main purpose |
| --- | --- |
| `001.strict_block_filter_for_default_switcherror.py` | Filters strict child polyphase blocks and splits blocks with potential switch errors. |
| `002.mat_block_link_relationship.py` | Compares original maternal phase blocks with ShapeIt4 phasing and summarizes block-level consistency. |
| `003.vcf_gt_clerify_pat_hap.py` | Clarifies paternal haplotype contribution using trio genotypes and child phase blocks. |
| `003_1.ori_SNP.stat.py` | Summarizes SNP categories before paternal clarification. |
| `003_2.afterstep1_SNP.stat.py` | Summarizes SNP categories after paternal clarification. |
| `004.vcf_gt_clerify_mat_hap.fix.v2.py` | Clarifies maternal haplotype contribution after paternal haplotype assignment. |
| `005.merge_pre.py` | Prepares SNP/block information for final block merging and VCF generation. |
| `006.make_new_chd_vcf.py` | Generates the final child haplotype VCF. |
| `block_interval_len.py` | Utility script for phase-block interval inspection. |

## Output directory layout

After a complete run, the working directory is expected to contain:

```text
${BASE_DIR}/
├── 00.data/
│   ├── mat_chr21_longreads.bam
│   ├── pat_chr21_longreads.bam
│   ├── chd_chr21_longreads.bam
│   ├── mat_chr21.g.vcf.gz
│   └── pat_chr21.g.vcf.gz
├── 01.jointcall/
│   └── chr21_mixed_ploidy/
│       └── ${FAMILY_ID}.chr21.joint.snps.filtered.pass.vcf.gz
├── 02.phase/
│   ├── ${FAMILY_ID}.family.phased.vcf.gz
│   ├── ${FAMILY_ID}.family.phased.son_het.vcf.gz
│   ├── *.block.tsv
│   ├── logs/
│   └── strict/
├── 03.shapeit_for_mat/
│   └── mat_block_link_relationship.tsv
└── 04.hand_phasing/
    ├── 01.first_step_clarify_pat/
    ├── 02.second_step_clarify_mat/
    └── 03.expand_SNP_and_merge_block/
        └── final.child.hap.vcf
```

## Notes and troubleshooting

### Sample names must match the VCF header

`CHILD_SAMPLE`, `MATERNAL_SAMPLE`, and `PATERNAL_SAMPLE` must exactly match the sample IDs in the VCF files. Check sample names with:

```bash
bcftools query -l input.vcf.gz
```

### GenomicsDBImport cannot overwrite an existing workspace

If Step 01 is rerun, remove the old GATK database directory first:

```bash
rm -r ${BASE_DIR}/01.jointcall/01.gatk_db/chr21_mixed_ploidy
```

### Use consistent chromosome naming

If the reference genome and VCF use `chr21`, all BAM, VCF, genetic map, and reference-panel files should also use `chr21`. If the ShapeIt4 genetic map uses `21`, convert it to `chr21` before running Step 03.

### The population reference panel is not fully included

The 1000 Genomes chr21 population panel can be large. Generate it locally following the commands in the `shapeit_data/` section.

## Minimal command summary

```bash
# 1. Edit parameters in all shell scripts first
vim 01.data_pre_and_jointcall.sh
vim 02.phase.sh
vim 03.shapeit_for_mat.sh
vim 04.hand_phasing.sh

# 2. Run all steps
bash total.sh

# 3. Check final output
ls 04.hand_phasing/03.expand_SNP_and_merge_block/final.child.hap.vcf
```

## Citation


