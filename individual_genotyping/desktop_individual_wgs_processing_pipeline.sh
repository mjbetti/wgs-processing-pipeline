#Installation of the following commandline tools is required: fastqc, bwa, picard, gatk, java (OpenJDK), ggplot2, r-gplots, r-gsalib, samtools, vcftools, and tabix. The easiest way to install all of these is via the included wgs_pipeline_env.yml Anaconda file.

#Broad Google Cloud Bucket: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
#The GRCh38 Reference Sequence used here was downloaded from https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta

#dbSNP file downloaded from dbSNP FTP (https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/00-All.vcf.gz)
#Mills, 1000G, omni, and HapMap files downloaded from Broad FTP (https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0)

#Declare all files and directories that will be used for inputs and outputs. Clinvar reference files can be used for downstream annotations with tools like ANNOVAR, but we will not use them here.

#For Dante Labs data, the read group(s) of your sample can be determined by running the following command with your hg19 aligned BAM file:
#   samtools view -H sample.bam | grep '@RG'
FASTQ1=$1
FASTQ2=$2
OUT_PREF=$3
MAIN_OUT_DIR=$4
REF_GENOME=$5
READ_GROUPS=$6
TMP_DIR=$7
INTER_DIR=$8
DBSNP=$9
MILLS=$10
SNPS1000G=$11
OMNI=$12
HAPMAP=$13
THREADS=$14
RAM=$15

#Index the downloaded refrence genome using BWA
echo 'Indexing reference genome...'
bwa index -a bwtsw $REF_GENOME

#Align the FASTQ files using BWA-MEM
echo 'Aligning reads to the reference genome...'
mkdir $INTER_DIR

bwa mem \
-t $THREADS \
-T 0 \
-R $READ_GROUPS \
$REF_GENOME \
$FASTQ1 \
$FASTQ2 | \
samtools view \
-Shb \
-o $INTER_DIR\/$OUT_PREF\.bam

#Sort the aligned BAM file using Picard tools
echo 'Sorting aligned BAM...'
mkdir $TMP_DIR
ulimit -c unlimited
picard SortSam "-Xmx"$RAM\g  \
    CREATE_INDEX=true \
    INPUT=$INTER_DIR\/$OUT_PREF\.bam \
    OUTPUT=$INTER_DIR\/$OUT_PREF\.sorted.bam \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=STRICT \
    TMP_DIR=$TMP_DIR

#Merge reads from multiple runs
echo 'Merging reads from multiple runs...'
ulimit -c unlimited
picard MergeSamFiles "-Xmx"$RAM\g \
    ASSUME_SORTED=false \
    CREATE_INDEX=true \
    INPUT=$INTER_DIR\/$OUT_PREF\.sorted.bam \
    MERGE_SEQUENCE_DICTIONARIES=false \
    OUTPUT=$INTER_DIR\/$OUT_PREF\.sorted.merged.bam \
    SORT_ORDER=coordinate \
    USE_THREADING=true \
    VALIDATION_STRINGENCY=STRICT

samtools sort -m $(($RAM / 4))\g -@ $THREADS \
    -o $INTER_DIR\/$OUT_PREF\.sorted.merged.sorted.bam \
    $INTER_DIR\/$OUT_PREF\.sorted.merged.bam

#Identify duplicate reads originating from the same single fragment of DNA, i.e. PCR artifacts
echo 'Identifying PCR duplicates...'
export _JAVA_OPTIONS=-Djava.io.tmpdir=$TMP_DIR
ulimit -c unlimited
gatk MarkDuplicatesSpark \
    -I $INTER_DIR\/$OUT_PREF\.sorted.merged.sorted.bam \
    -O $INTER_DIR\/$OUT_PREF\.sorted.merged.sorted.marked_duplicates.bam \
    -M $INTER_DIR\/$OUT_PREF\.sorted.merged.sorted.marked_duplicates_metrics.txt

#Generate an index and dictionary for the reference genome for use in gatk processing
#echo 'Creating reference genome index and dict for GATK...'
gatk CreateSequenceDictionary \
    --REFERENCE $REF_GENOME

samtools faidx $REF_GENOME

#Before generating the table, the VCF files will need corresponding indices generated
echo 'Indexing reference VCF files...'
gatk IndexFeatureFile --input $DBSNP
gatk IndexFeatureFile --input $MILLS
gatk IndexFeatureFile --input $SNPS1000G
gatk IndexFeatureFile --input $OMNI
gatk IndexFeatureFile --input $HAPMAP

#Generate a recalibration table for Base Quality Score Recalibration
echo 'Generating base quality score recalibration table...'
gatk BaseRecalibrator \
    -I $INTER_DIR\/$OUT_PREF\.sorted.merged.sorted.marked_duplicates.bam \
    -R $REF_GENOME \
    --known-sites $DBSNP \
    --known-sites $MILLS \
    --known-sites $SNPS1000G \
    -O $INTER_DIR\/$OUT_PREF\_recal_data.table

#Apply the calculated base quality score recalibration
echo 'Applying calculated recalibration scores...'
gatk ApplyBQSR \
    -R $REF_GENOME \
    -I $INTER_DIR\/$OUT_PREF\.sorted.merged.sorted.marked_duplicates.bam \
    --bqsr-recal-file $INTER_DIR\/$OUT_PREF\_recal_data.table \
    -O $MAIN_OUT_DIR\/$OUT_PREF\.sorted.merged.sorted.marked_duplicates.recalibrated.bam

#Evaluate the base quality score recalibration tables
echo 'Evaluating base quality score recalibration tables...'
gatk AnalyzeCovariates \
-bqsr $INTER_DIR\/$OUT_PREF\_recal_data.table \
-plots $INTER_DIR\/$OUT_PREF\_AnalyzeCovariates_mb.pdf

#Somatic short variant discovery (SNVs + Indels)
#Use HaplotypeCaller to call germline SNPs and indels via local re-assemply of haplotypes. The output will be an intermediate GVCF file, containing raw, unfiltered SNP and indel calls
echo 'Initial calling of SNPs and indels...'
ulimit -c unlimited
export _JAVA_OPTIONS=-Djava.io.tmpdir=$TMP_DIR
gatk --java-options "-Xmx"$RAM\g HaplotypeCaller \
    -R $REF_GENOME \
    -I $MAIN_OUT_DIR\/$OUT_PREF\.sorted.merged.sorted.marked_duplicates.recalibrated.bam \
    -O $INTER_DIR\/$OUT_PREF\.g.vcf \
    -ERC GVCF

bgzip $INTER_DIR\/$OUT_PREF\.g.vcf

#Perform genotyping on one or more samples pre-called with HaplotypeCaller
echo 'Genotyping pre-called samples...'
ulimit -c unlimited
export _JAVA_OPTIONS=-Djava.io.tmpdir=$TMP_DIR
gatk --java-options "-Xmx"$RAM\g \
GenotypeGVCFs \
-R $REF_GENOME \
-V $INTER_DIR\/$OUT_PREF\.g.vcf \
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
bgzip $MAIN_OUT_DIR\/$OUT_PREF\.recal.snp.indel.vcf
tabix -f -p vcf $MAIN_OUT_DIR\/$OUT_PREF\.recal.snp.indel.vcf.gz

#Split the final VCF output into two separate SNP and indel VCF files
echo 'Splitting VCF into separate SNP and indel files...'
vcftools \
	--gzvcf $MAIN_OUT_DIR\/$OUT_PREF\.recal.snp.indel.vcf.gz \
	--remove-indels \
	--recode \
	--recode-INFO-all \
	--stdout > \
	$MAIN_OUT_DIR\/$OUT_PREF\.snp.vcf

bgzip $MAIN_OUT_DIR\/$OUT_PREF\.snp.vcf
tabix -f -p vcf $MAIN_OUT_DIR\/$OUT_PREF\.snp.vcf.gz

vcftools \
	--gzvcf $MAIN_OUT_DIR\/$OUT_PREF\.recal.snp.indel.vcf.gz \
	--keep-only-indels \
	--recode \
	--recode-INFO-all \
	--stdout > \
	$MAIN_OUT_DIR\/$OUT_PREF\.indel.vcf

bgzip $MAIN_OUT_DIR\/$OUT_PREF\.indel.vcf
tabix -f -p vcf $MAIN_OUT_DIR\/$OUT_PREF\.indel.vcf