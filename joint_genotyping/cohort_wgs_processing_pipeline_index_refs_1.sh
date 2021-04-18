#Installation of the following commandline tools is required: fastqc, bwa, picard, gatk, java (OpenJDK), ggplot2, r-gplots, r-gsalib, samtools, vcftools, and tabix. The easiest way to install all of these is via the included wgs_pipeline_env.yml Anaconda file.

#Broad Google Cloud Bucket: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
#The GRCh38 Reference Sequence used here was downloaded from https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta

#dbSNP file downloaded from dbSNP FTP (https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/00-All.vcf.gz)
#Mills, 1000G, omni, and HapMap files downloaded from Broad FTP (https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0)

#Declare all files and directories that will be used for inputs and outputs. Clinvar reference files can be used for downstream annotations with tools like ANNOVAR, but we will not use them here.

#For Dante Labs data, the read group(s) of your sample can be determined by running the following command with your hg19 aligned BAM file:
#   samtools view -H sample.bam | grep '@RG'
REF_GENOME=$1
TMP_DIR=$2
INTER_DIR=$3
DBSNP=$4
MILLS=$5
SNPS1000G=$6
OMNI=$7
HAPMAP=$8

#Index the downloaded refrence genome using BWA
echo 'Indexing reference genome...'
bwa index -a bwtsw $REF_GENOME

#Create a directory for intermediate files to be stored
echo 'Making directory for intermediate files...'
mkdir $INTER_DIR

#Create a directory for temporary files to be stored
echo 'Making directory for temporary files...'
mkdir $TMP_DIR

#Generate an index and dictionary for the reference genome for use in gatk processing
#echo 'Creating reference genome index and dict for GATK...'
echo 'Creating index and dictionary for reference genome...'
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