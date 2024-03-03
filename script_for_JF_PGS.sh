#!/bin/bash

PHEN_NAME=BMI_10

working_dir=/your/working/directory
reference_panel=${working_dir}/reference_panel
clumping_output_dir=${working_dir}/clumped/${PHEN_NAME}

mkdir -p ${clumping_output_dir}
mkdir -p ${clumping_output_dir}/tmp



# Clump jackknifing gwas
for i in {1..10}
do

## specify the path of the GWAS summary statistics
GWAS_FILE=${working_dir}/data/UKBB/GWAS_results/${PHEN_NAME}_${i}_imputed.txt.gz

## extract the RSIDs and p-values from GWAS summary statistics for clumping
zcat ${GWAS_FILE} | awk '{print $1"\t"$2}' > ${clumping_output_dir}/tmp/${PHEN_NAME}_all_${i}_gwas_pval.txt

## perform GWAS clumping
plink \
    --bfile $reference_panel \
    --clump ${clumping_output_dir}/tmp/${PHEN_NAME}_all_${i}_gwas_pval.txt \
    --clump-field P_BOLT_LMM_INF \
    --clump-p1 5e-08 \
    --clump-p2 5e-08 \
    --clump-kb 1000 \
    --clump-r2 0.001 \
    --out ${clumping_output_dir}/tmp/${PHEN_NAME}_all_${i}_gwas

## extract the RSIDs of SNPs remained after clumping
awk '{print $3}' ${clumping_output_dir}/tmp/${PHEN_NAME}_all_${i}_gwas.clumped | sed '/^\s*$/d' > ${clumping_output_dir}/tmp/${PHEN_NAME}_all_${i}_gwas.clumped.snps

## extract the clumped SNPs from the GWAS summary statistics
zcat ${GWAS_FILE} | grep -wFf ${clumping_output_dir}/tmp/${PHEN_NAME}_all_${i}_gwas.clumped.snps > ${clumping_output_dir}/${PHEN_NAME}_all_${i}_gwas_instruments.txt

done


# Extracting genotype information

temp_geno_prefix=${PHEN_NAME}
bgen_pattern=/path/to/genotype/bgen/file/data.chrCHROM.bgen # genotype files by chromosome
bgen_index_pattern=/path/to/genotype/bgen/index/file/data.chrCHROM.bgen.bgi # genotype index files by chromosome
resource_dir=/path/to/SNP_extraction_files # a folder containing some necessary files (see below)
genotype_output_dir=${working_dir}/genotype/${PHEN_NAME} # save a temporary genotype file to avoid loading SNPs in PGS calculation

mkdir -p ${genotype_output_dir}

sort -u ${clumping_output_dir}/tmp/${PHEN_NAME}*.clumped.snps > ${clumping_output_dir}/tmp/${PHEN_NAME}_clumped_snp_ids.txt

for chrom in {1..22}; do
    chrom_padd=$(printf "%0*d\n" 2 $chrom)
    inbgen=${bgen_pattern/CHROM/$chrom_padd}
    inbgenidx=${bgen_index_pattern/CHROM/$chrom_padd}
    bgenix \
        -g $inbgen \
        -i $inbgenidx \
        -incl-rsids ${clumping_output_dir}/tmp/${PHEN_NAME}_clumped_snp_ids.txt \
        > ${genotype_output_dir}/${temp_geno_prefix}.${chrom_padd}.bgen
done

cmd=""
for chrom in {01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22}; do
        cmd="${cmd} ${genotype_output_dir}/${temp_geno_prefix}.${chrom}.bgen"
done

cat-bgen \
    -g ${cmd} \
    -clobber \
    -og ${genotype_output_dir}/${PHEN_NAME}_genotype.bgen

plink2 \
--bgen ${genotype_output_dir}/${PHEN_NAME}_genotype.bgen \
# specify all sample IDs
--sample ${resource_dir}/data.all.sample \
# remove related individuals
--remove ${resource_dir}/strict_genetic_removes_plink_new.txt \
--hard-call-threshold 0.4999 \
--make-bed \
--out ${genotype_output_dir}/tmp


# Calculate PGS

BFILE=${genotype_output_dir}/tmp
cut -f2 ${BFILE}.bim | uniq -d > ${working_dir}/data/UKBB/genotype/$PHEN_NAME/dups

OUTPUT_DIR=$working_dir/data/UKBB/prs/$PHEN_NAME
mkdir -p $OUTPUT_DIR

# Jackknifing prs for all (males and females)

for i in {1..10}
do

IND_LIST=/path/to/phenotype/file/${PHEN_LABEL}/${PHEN_NAME}_grp_${i}.txt # phenotype for individual in each block
SCORE_FILE=${clumping_output_dir}/${PHEN_NAME}_${i}_gwas_instruments.txt

plink --bfile ${BFILE} \
--exclude ${working_dir}/data/UKBB/genotype/$PHEN_NAME/dups \
--score ${SCORE_FILE} 1 5 11 header sum \
--keep ${IND_LIST} \
--out ${OUTPUT_DIR}/${PHEN_NAME}_${i}_prs

done

awk 'FNR>1 || NR==1' ${OUTPUT_DIR}/${PHEN_NAME}_*_prs.profile > ${OUTPUT_DIR}/${PHEN_NAME}_jackknifing_prs.profile

