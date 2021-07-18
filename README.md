# wgs-processing-pipeline
## Introduction
This pipeline is intended for use in the initial processing (alignment and variant calling) of whole genome sequencing data. Two sets of pipeline scripts are included, one for performing variant calling on a single sample and another set for performing joing variant calling on a set of multiple samples. The sequential workflow used here is based on guidelines outlined in the GATK Best Practices Workflow (https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows).

### Variant calling pipeline
The included ```desktop_individual_wgs_processing_pipeline_individual.sh``` or ```cluster_individual_wgs_processing_pipeline_individual.sh``` scripts (depending on one's computational setup) are the simplest ways to run the pipeline. Details can be found under the Included Scripts section, including documentation for additional scripts provided for joint genotyping (if processing multiple samples simultaneously).

and then proceeds to map the paired reads to the most current genome build (GRCh38/hg38) using ```bwa-mem```. After some downstream processing, variants (SNPs and indels) are called using GATK's ```HaplotypeCaller```. The resulting single VCF can optionally be split into separate SNP and indel files.

### Variant Annotation
The workflow for setting up and running GenomeChronicler is based off of the turotials presented by the Personal Genome Project UK (https://github.com/PGP-UK/GenomeChronicler) and Singularity (https://sylabs.io/guides/3.1/user-guide/). This tool takes in aligned reads (as a BAM file) and generates a PDF report highlighting the clinical significance of a panel of SNPs found in the proband's WGS data.

### Future Additions
Genetic ancestry estimates, copy number variant (CNV) calling, and fine-grain variant annotation are will hopefully be added to this workflow eventually.

## Getting started
### Installing dependencies
The easiest way to install all of the required tools is via a package manager such as Anaconda (https://docs.conda.io/en/latest/miniconda.html). For your convenience, a ready-to-use Anaconda environment can be installed using the ```wgs_pipeline_env.yml``` file. This environment contains all of the following dependencies:

* ```fastqc```
* ```bwa```
* ```picard```
* ```gatk4```
* ```java (OpenJDK)```
* ```R (r-base)```
* ```ggplot2```
* ```r-gplots```
* ```r-gsalib```
* ```samtools```
* ```vcftools```
* ```tabix```
* ```bcftools```

Assuming you have Anaconda already installed on your machine, the environment can be compiled with the following command:
```conda env create -f wgs_pipeline_env.yml```

### Downloading reference files
Most of the required reference files can be downloaded from the Broad Institute's Google Cloud Bucket (https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/). The dbSNP reference panel was downloaded from downloaded from https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/

In total, we will need the following files:
* GRCh38 (hg38) reference genome FASTA (```resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta```)
* Mills and 1000 Genomes gold standard indels (```resources-broad-hg38-v0-Mills_and_1000G_gold_standard.indels.hg38.vcf```)
* 1000 Genomes Phase 1 high-confidence SNPs (```resources-broad-hg38-v0-1000G_phase1.snps.high_confidence.hg38.vcf```)
* Omni reference panel (```resources-broad-hg38-v0-1000G_omni2.5.hg38.vcf```)
* HapMap reference panel (```resources-broad-hg38-v0-hapmap_3.3.hg38.vcf```)
* dbSNP reference panel (```00-All.vcf.gz```)

There are variable assignments in the script for ClinVar and dbSNP reference files, as well, but these are not actually used in this workflow.

## Usage
Optional: If you are unsure about or interested in assessing the quality of your raw reads prior to alignment, you can so do using ```fastqc```. a tool included in the Anaconda environment:
```
FASTQ1=/path/to/forward/read
FASTQ2=/path/to/reverse/read
FASTQC_OUT=/path/to/desired/output directory (will be created by fastqc)
fastqc -o $FASTQC_OUT $FASTQ1 $FASTQ2
```

### Included scripts
This repository contains several variations of the core pipeline script:
* The ```individual_genotyping``` directory contains scripts that should be used when working with sequencing data from a single individual, meaning that you will not be performing joing genotyping. The ```cluster_individual_wgs_processing_pipeline.sh``` script is optimized for running on a UNIX-based cluster, while the ```desktop_individual_wgs_processing_pipeline.sh``` should be used if you are running this pipeline on a desktop computer. It has only been tested on macOS and Linux, although it should also work in a UNIX-based terminal on Windows. If you are using a cluster, the pipeline script should be run interactively (as opposed to as a batch submission), as some of the GATK tools will otherwise not run properly.
  * You can run this script by specifying the required command line arguments, which will be read in by the parser. The cluster-optimized script can be run in exactly the same way, except that the RAM argument does not need to be specified.
    ```
    FASTQ1=/path/to/forward/read
    FASTQ2=/path/to/reverse/read
    OUT_PREF=string (desired prefix for all output files)
    MAIN_OUT_DIR=/path/to/desired/root/directory/for/outputs
    REF_GENOME=path/to/reference/fasta
    READ_GROUPS=string (information on required format found at https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups, use "samtools view -H aligned_hg19.bam | grep '@RG" to find read group names in an existing BAM file) 
    TMP_DIR=$MAIN_OUT_DIR\/path/to/directory/for/temp/files (will be created by the pipeline script)
    INTER_DIR=$MAIN_OUT_DIR\/path/to/store/intermediate/files/during/processing (will be created by the pipeline script)
    DBSNP=/path/to/dbsnp/reference/file
    MILLS=/path/to/mills/reference/file
    SNPS1000G=/path/to/1000_genomes/reference/file
    OMNI=/path/to/omni/reference/file
    HAPMAP=/path/to/hapmap/reference/file
    THREADS=int (number of CPU threads to use)
    RAM=int (amount of RAM to use in GB)
    
    desktop_individual_wgs_processing_pipeline.sh \
    $FASTQ1 \
  	$FASTQ2 \
  	$OUT_PREF \
  	$MAIN_OUT_DIR \
  	$REF_GENOME \
  	$READ_GROUPS \
  	$TMP_DIR \
  	$INTER_DIR \
  	$DBSNP \
  	$MILLS \
  	$SNPS1000G \
  	$OMNI \
  	$HAPMAP \
  	$THREADS \
  	$RAM

    ```
* The ```joint_genotyping``` directory, as its name suggests, contains pipeline scripts that should be used if you will be performing joing genotyping using multiple samples (i.e. sequencing data from more than one individual). Like the scripts for individual genotyping, the scripts beginning with "cluster" are optimized for running on a UNIX-based cluster, while the scripts designated "desktop" are intended for running on a desktop computer. The ```cohort_wgs_processing_pipeline_index_refs_1.sh``` script should be run first regardless of system, and then the remaining two scripts should be run based on whether you are using a cluster or desktop computer. As an example, the scripts should be run in the following sequence:
  * The initial ```cohort_wgs_processing_pipeline_index_refs_1.sh``` script will generate all of the required indices and dictionaries for all of your reference files. This is to avoid the process being repeated over and over again if you will be aligning and calling variants for multiple samples in parallel.
    ```
    REF_GENOME=path/to/reference/fasta
    TMP_DIR=$MAIN_OUT_DIR\/path/to/directory/for/temp/files (will be created by the pipeline script)
    INTER_DIR=$MAIN_OUT_DIR\/path/to/store/intermediate/files/during/processing (will be created by the pipeline script)
    DBSNP=/path/to/dbsnp/reference/file
    MILLS=/path/to/mills/reference/file
    SNPS1000G=/path/to/1000_genomes/reference/file
    OMNI=/path/to/omni/reference/file
    HAPMAP=/path/to/hapmap/reference/file

    cohort_wgs_processing_pipeline_index_refs_1.sh \
      $REF_GENOME \
      $TMP_DIR \
      $INTER_DIR \
      $DBSNP \
      $MILLS \
      $SNPS1000G \
      $OMNI \
      $HAPMAP
    ```
  * The second ```desktop_cohort_wgs_processing_pipeline_index_refs_2.sh``` script will run alignment and generation of the initial GVCF files. This script should be run individually (either in parallel or sequentially) for each individual sample. The ```cluster_cohort_wgs_processing_pipeline_index_refs_2.sh``` script optimized for use on a cluster would be run in the same way, except that the final RAM argument should not be specified.
    ```
    FASTQ1=/path/to/forward/read
    FASTQ2=/path/to/reverse/read
    MAIN_OUT_DIR=/path/to/desired/root/directory/for/outputs
    OUT_PREF=string (desired prefix for all output files)
    REF_GENOME=/path/to/reference/fasta
    READ_GROUPS=string (information on required format found at https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups)
    TMP_DIR=$MAIN_OUT_DIR\/path/to/directory/for/temp/files (will be created by the pipeline script)
    INTER_DIR=$MAIN_OUT_DIR\/path/to/store/intermediate/files/during/processing (will be created by the pipeline script)
    DBSNP=/path/to/dbsnp/reference/file
    MILLS=/path/to/mills/reference/file
    SNPS1000G=/path/to/1000_genomes/reference/file
    OMNI=/path/to/omni/reference/file
    HAPMAP=/path/to/hapmap/reference/file
    THREADS=int (number of CPU threads to use)
    RAM=int (amount of RAM to use in GB)
    
    desktop_cohort_wgs_processing_pipeline_index_refs_2.sh \
      $FASTQ1 \
      $FASTQ2 \
      $MAIN_OUT_DIR \
      $OUT_PREF \
      $REF_GENOME \
      $READ_GROUPS \
      $TMP_DIR \
      $INTER_DIR \
      $DBSNP \
      $MILLS \
      $SNPS1000G \
      $OMNI \
      $HAPMAP \
      $THREADS \
      $RAM
     ```
  * The third ```desktop_cohort_wgs_processing_pipeline_index_refs_3.sh``` script should be run once all GVCFs (for all samples) have been generated. This script performs joing variant calling, as well as downstream variant quality score recalibration for both SNPs and indels. By default, this script is written to take in GVCF files for 5 samples, but the argument parser can be modified to accept any number of samples. The ```cluster_cohort_wgs_processing_pipeline_index_refs_3.sh``` script optimized for use on a cluster would be run in the same way, except that the final RAM argument should not be specified.
    ```
    GVCF1=/path/to/sample1/gvcf
    GVCF2=/path/to/sample2/gvcf
    GVCF3=/path/to/sample3/gvcf
    GVCF4=/path/to/sample4/gvcf
    GVCF5=/path/to/sample5/gvcf
    OUT_PREF=string (desired prefix for all output files)
    REF_GENOME=/path/to/reference/fasta
    INTER_DIR=/path/to/store/intermediate/files/during/processing (will be created by the pipeline script)
    DBSNP=/path/to/dbsnp/reference/file
    MILLS=/path/to/mills/reference/file
    SNPS1000G=/path/to/1000_genomes/reference/file
    OMNI=/path/to/omni/reference/file
    HAPMAP=/path/to/hapmap/reference/file
    THREADS=int (number of CPU threads to use)
    RAM=int (amount of RAM to use in GB)
     
    desktop_cohort_wgs_processing_pipeline_index_refs_2.sh \
      $GVCF1 \
      $GVCF2 \
      $GVCF3 \
      $GVCF4 \
      $GVCF5 \
      $OUT_PREF \
      $REF_GENOME \
      $INTER_DIR \
      $DBSNP \
      $MILLS \
      $SNPS1000G \
      $OMNI \
      $HAPMAP \
      $THREADS \
      $RAM
    ```

Once the pipeline has completed running, the files you will likely want to retain for downstream analysis will be
* ```$MAIN_DIR\/$OUT_PREF\.sorted.merged.sorted.marked_duplicates.recalibrated.bam``` (BAM containing your aligned reads)
* ```$MAIN_DIR\/$OUT_PREF\.recal.snp.indel.vcf.gz``` (recalibrated VCF cotaining SNP and indel calls)

## Initial variant annotation using GenomeChronicler
This GenomeChronicler workflow is best thought of as an interesting exploratory analysis, generating a neat PDF report of GWAS associations with common variants present in the input BAM file. It is likely not applicable to most research contexts and should be considered separate from the main variant-valling workflow detailed above.

### Installing Golang and Singlularity (summarized from Golang and Singularity documentation)
The easiest way to run GenomeChronicler is via a Singularity container, which in turn requires a Golang installation. Running Singularity containers is very similar to Docker, and more information about Singularity can be found at the following link: https://sylabs.io/guides/3.1/user-guide/quick_start.html. Go installation documentation: https://golang.org/doc/install. 

GenomeChronicler only seems to work with human genome build GRCh38, so acces to a wider set of downstream analysis tools is a tangible benefit of realigning your reads to the latest build.

For use with an Ubuntu Linux installation, one first needs to insure that all required development tools and libraries are installed on the system:
```
sudo apt-get update && \
sudo apt-get install -y build-essential \
libseccomp-dev pkg-config squashfs-tools cryptsetup
```
Install Golang if not already on your system. If you are updating from a older version, remove ```/usr/local/go``` before reinstalling the latest version:
```
export VERSION=1.14.9 OS=linux ARCH=amd64  #alter these variables as required to align with user system
wget -O /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz https://dl.google.com/go/go${VERSION}.${OS}-${ARCH}.tar.gz && \
  sudo tar -C /usr/local -xzf /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz
```
Set up environment to use Go (add Go path to ```.bashrc``` file):
```
echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
  echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
  source ~/.bashrc
```
Install golangci-lint:
```
curl -sfL https://install.goreleaser.com/github.com/golangci/golangci-lint.sh |
  sh -s -- -b $(go env GOPATH)/bin v1.21.0
```
Clone the Singularity git repository from the source:
```
mkdir -p ${GOPATH}/src/github.com/sylabs && \
  cd ${GOPATH}/src/github.com/sylabs && \
  git clone https://github.com/sylabs/singularity.git && \
  cd singularity
```
To build a stable version of Singularity, check out a release tag before compiling:
```
git checkout v3.6.3
```
Build Singularity from source code:
```
cd ${GOPATH}/src/github.com/sylabs/singularity && \
  ./mconfig && \
  cd ./builddir && \
  make && \
  sudo make install
```
You can check the specific version of your Singularity install using the following command:
```
singularity version
```
Once Singularity is built, clone the compiled GenomeChronicler git repository into your home directory:
```
git clone https://github.com/PGP-UK/GenomeChronicler.git
```

### Installing and running GenomeChronicler
Download the pre-packaged GenomeChronicler image from SingularityHub:
```
cd GenomeChronicler
singularity pull shub://PGP-UK/GenomeChronicler
```
After cloning this repository, run the ```SetupMeFirst.sh``` script in your local system to retrieve the extra files needed to run the pipeline (around 10GB total, so too large to have included in the git repository).
```
~/GenomeChronicler/SetupMeFirst.sh
```
Because the input BAM appears to need to be in the same directory as the ```GenomeChronicler_latest.sif``` Singularity file, the aligned BAM file and its corresponding index (ending in ```.bai```) should next be copied to the ```GenomeChronicler/``` directory. Attempting to use a symbolic link instead to avoid duplicating files did not seem to work.

Once the BAM is copied over, navigate to the ```GenomeChronicler/``` directory and run the following command:
```
singularity run GenomeChronicler_latest.sif --bamFile=mb_hg38_60820188479382.sorted.merged.sorted.marked_duplicates.recalibrated.bam
```

This command should generate a new results directory inside of ```GenomeChronicler/results/``` containing the PDF variant annotations report, as well as corresponding ```.txt```, ```.xlsx```, and ```.vcf.gz``` files, such as the following path:
```
~/GenomeChronicler/results/results_mb_hg38_60820188479382.sorted.merged.sorted.marked_duplicatesibrated
```
