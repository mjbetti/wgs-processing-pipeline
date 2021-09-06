#Declare the path of the target directory to which the reference files will be downloaded
DOWN_DIR=${1}

#Download the hg38 FASTA
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta \
$DOWN_DIR

#Download the set of Mills and 1000 Genomes gold standard indels
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
$DOWN_DIR

#Download the set of 1000 Genomes Phase 1 high-confidence SNPs
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
$DOWN_DIR

#Download the Omni SNP reference panel
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz \
$DOWN_DIR

#Download the HapMap reference panel
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz \
$DOWN_DIR

#Download the dbSNP b151 reference file
wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/00-All.vcf.gz \
-P $DOWN_DIR