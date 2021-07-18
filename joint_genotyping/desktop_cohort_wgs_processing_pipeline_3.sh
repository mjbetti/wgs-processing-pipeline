#Installation of the following commandline tools is required: fastqc, bwa, picard, gatk, java (OpenJDK), ggplot2, r-gplots, r-gsalib, samtools, vcftools, and tabix. The easiest way to install all of these is via the included wgs_pipeline_env.yml Anaconda file.

#Broad Google Cloud Bucket: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
#The GRCh38 Reference Sequence used here was downloaded from https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta

#dbSNP file downloaded from dbSNP FTP (https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/00-All.vcf.gz)
#Mills, 1000G, omni, and HapMap files downloaded from Broad FTP (https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0)

#Declare all files and directories that will be used for inputs and outputs. Clinvar reference files can be used for downstream annotations with tools like ANNOVAR, but we will not use them here.

#For Dante Labs data, the read group(s) of your sample can be determined by running the following command with your hg19 aligned BAM file:
#   samtools view -H sample.bam | grep '@RG'

#Modify the number of GVCF arguments to account for the number of individuals that will be jointly genotyped. This script has arguments for 5 GVCFs included by default

GVCF1=${1}
GVCF2=${2}
GVCF3=${3}
GVCF4=${4}
GVCF5=${5}
OUT_PREF=${6}
REF_GENOME=${7}
READ_GROUPS=${8}
INTER_DIR=${9}
DBSNP=${10}
MILLS=${11}
SNPS1000G=${12}
OMNI=${13}
HAPMAP=${14}
THREADS=${15}
RAM=${16}

#Perform genotyping on one or more samples pre-called with HaplotypeCaller
echo 'Genotyping pre-called samples...'
ulimit -c unlimited
export _JAVA_OPTIONS=-Djava.io.tmpdir=$TMP_DIR
gatk --java-options "-Xmx"$RAM\g \
GenotypeGVCFs \
-R $REF_GENOME \
-V $GVCF1 \
-V $GVCF2 \
-V $GVCF3 \
-V $GVCF4 \
-V $GVCF5 \
-O $INTER_DIR\/$OUT_PREF\.vcf

bgzip $INTER_DIR\/$OUT_PREF\.vcf

#Build a SNP recalibration model to score variant quality and then apply it to filter the SNPS in the generated VCF
#Referred to WGS recalibration code on (https://github.com/BD2KGenomics/gatk-whole-genome-pipeline/blob/master/HAPvariantCalling.sh)
echo 'Building SNP recalibration model...'
gatk VariantRecalibrator \
    -R $REF_GENOME \
    -V $INTER_DIR\/$OUT_PREF\.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 $OMNI \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 $SNPS1000G \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP\
    -an QD \
    -an DP \
    -an FS \
    -an ReadPosRankSum \
    --mode SNP \
    --output $INTER_DIR\/$OUT_PREF\.snp.recal \
    --tranches-file $INTER_DIR\/$OUT_PREF\.snp.tranches \
    --rscript-file $INTER_DIR\/$OUT_PREF\.snp.plots.R

echo 'Applying SNP recalibration model...'
gatk ApplyVQSR \
    -R $REF_GENOME \
    -V $INTER_DIR\/$OUT_PREF\.vcf.gz \
    --output $INTER_DIR\/$OUT_PREF\.recal.snp.vcf \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file $INTER_DIR\/$OUT_PREF\.snp.tranches \
    --recal-file $INTER_DIR\/$OUT_PREF\.snp.recal \
    -mode SNP

bgzip $INTER_DIR\/$OUT_PREF\.recal.snp.vcf

#Perform a second round of recalibration, this time focusing on the indels
echo 'Building indel recalibration model...'
gatk VariantRecalibrator \
    -R $REF_GENOME \
    -V $INTER_DIR\/$OUT_PREF\.recal.snp.vcf.gz \
    --resource:mills,known=true,training=true,truth=true,prior=12.0 $MILLS \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
    -an DP \
    -an FS \
    -an ReadPosRankSum \
    --mode INDEL \
    --output $INTER_DIR\/$OUT_PREF\.snp.indel.recal \
    --tranches-file $INTER_DIR\/$OUT_PREF\.snp.indel.tranches \
    --rscript-file $INTER_DIR\/$OUT_PREF\.snp.indel.plots.R

echo 'Applying indel recalibration model...'
gatk ApplyVQSR \
    -R $REF_GENOME \
    -V $INTER_DIR\/$OUT_PREF\.vcf.gz \
    --output $MAIN_DIR\/$OUT_PREF\.recal.snp.indel.vcf \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file $INTER_DIR\/$OUT_PREF\.snp.indel.tranches \
    --recal-file $INTER_DIR\/$OUT_PREF\.snp.indel.recal \
    -mode INDEL

echo 'Compressing and indexing SNP/indel VCF...'
bgzip $MAIN_DIR\/$OUT_PREF\.recal.snp.indel.vcf
tabix -f -p vcf $MAIN_DIR\/$OUT_PREF\.recal.snp.indel.vcf.gz

#Split the final VCF output into two separate SNP and indel VCF files
#echo 'Splitting VCF into separate SNP and indel files...'
#vcftools \
#	--gzvcf $MAIN_DIR\/$OUT_PREF\.recal.snp.indel.vcf.gz \
#	--remove-indels \
#	--recode \
#	--recode-INFO-all \
#	--stdout > \
#	$MAIN_DIR\/$OUT_PREF\.snp.vcf

#bgzip $MAIN_DIR\/$OUT_PREF\.snp.vcf
#tabix -f -p vcf $MAIN_DIR\/$OUT_PREF\.snp.vcf.gz

#vcftools \
#	--gzvcf $MAIN_DIR\/$OUT_PREF\.recal.snp.indel.vcf.gz \
#	--keep-only-indels \
#	--recode \
#	--recode-INFO-all \
#	--stdout > \
#	$MAIN_DIR\/$OUT_PREF\.indel.vcf

#bgzip $MAIN_DIR\/$OUT_PREF\.indel.vcf
#tabix -f -p vcf $MAIN_DIR\/$OUT_PREF\.indel.vcf