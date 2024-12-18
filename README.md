# Overview
This snakemake workflow will convert Illumina genotyping .idat files to .gtc, then to .vcf format. This workflow requires that you have the Illumina IAAP-CLI software installed, the appropriate Illumina manifest files (.bpm and .csv format) and clusterfile for the genotyping array downloaded from the illumina website, and a reference genome in fasta format. You will also need to have Conda or Mamba installed. Instructions for installing the Illumina IAAP-CLI software and reference genome GRCh38/hg38 were modified from this [gtc2vcf README](https://github.com/freeseek/gtc2vcf/blob/1898320dab37c9da4e355f0fa31d2ab28d3632d5/README.md#installation) created by GitHub user freeseek. 

## Dependencies:
  * Illumina IAAP-CLI
  * Conda/Mamba

## Required Support Files:
  * Illumina manifest.bpm
  * Illumina manifest.csv
  * Illumina clusterfile.egt
  * reference genome in fasta format

## Install Illumina IAAP-CLI 
Download the appropriate installer for the [Illumina IAAP-CLI software](https://support.illumina.com/downloads/iaap-genotyping-cli.html). If you need to install the software on a remote server, `scp` the installer to the host server, then continue with the installation:

```shell
# Create a bin directory if one does not already exist.
mkdir ~/bin
# Extract the installer from the tarball. 
tar -xzvf iaap-cli-linux-x64-1.1.0-sha.80d7e5b3d9c1fdfc2e99b472a90652fd3848bbc7.tar.gz -C ~/bin iaap-cli-linux-x64-1.1.0-sha.80d7e5b3d9c1fdfc2e99b472a90652fd3848bbc7/iaap-cli/ --strip-components=1
# Remove the tarballed installer
rm iaap-cli-linux-x64-1.1.0-sha.80d7e5b3d9c1fdfc2e99b472a90652fd3848bbc7.tar.gz
# Create an alias for the iaap executable.
echo "export PATH=$PATH:~/bin/iaap-cli/" >> ~/.bash_profile
source ~/.bash_profile
# Test that the tool installed properly.
iaap gencall -h
```
## Install Reference Genome
If you would prefer to install GRCh37/hg19 instead, an ftp link and guide on which reference genome build to use can be found [here](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use).

### Dependencies
 * samtools
 * bwa

If you do not already have these tools installed, you can install them with Conda/Mamba (replace `mamba` with `conda` if you don't have mamba installed):
```shell
# Not a requirement, but encouraged to create a new environment.
mamba create -n genetics -y
mamba activate genetics
# Install samtools and bwa.
mamba install bioconda::samtools
mamba install bioconda::bwa
```
Then run the following:
```shell
# Make directories for reference_genomes and build version.
mkdir ref_genomes/
mkdir ref_genomes/GRCh38/
cd ref_genomes/GRCh38/
# Download the reference genome.
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
# Unzip the reference genome.
gzip -d GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
# Index the reference genome in fasta format.
samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
# Index reference genome for alignment.
bwa index GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```
Deactivate conda/mamba environment:
```shell
mamba deactivate 
```

# Using this Workflow
## Install snakemake 
If you already have a conda environment for running snakemake, activate the environment and proceed to the next step. Otherwise, follow the below commands to install snakemake from the env.yml. This yml will also install the snakemake executor plugin 'cluster-generic' which allows snakemake to submit jobs to a cluster scheduler.

```shell
# Create environment from yml file
mamba env create -f env.yml
# Activate environment
mamba activate snakemake
```
## Edit Snakemake Config
The Snakemake configuration file (`workflow/config.yml`) will need to be modified so that the input paths, output filenames, and paramaters are specific to your computing environment and resources available. The provided samples_file should be a text file of the samples in the format of `<sentrix_chip_number>_<sentrix_chip_position>` with one sample per line that you want processed. The workflow will not work on samples with incomplete data (both red and green idat files must exist). 

If you don't have a meta-data file of the samples in the raw_data_directory beforehad to make the necessary samples_file, you can use the included `workflow/make_samples_file.smk` to generate two text files: `samples.txt` and `incomplete_samples.txt`. This accessory snakemake workflow uses the same config file as the main workflow, but the samples_file parameter isn't used and can thus can be set as an empty string or an arbitrary filename. The `samples.txt` is what you will provide as the samples_file to the main snakemake workflow. The `incomplete_samples.txt` will include the names of any samples that did not have both red and green idat files (if the file is empty there were no incomplete samples). 

If you need to create the samples_file for the main snakemake workflow, run the following:
```shell
snakemake -j 1 -s workflow/make_samples_file.smk
```
Then, add the path for the output `samples.txt` to the config_file and run the main workflow. 

## Run Snakemake Locally
The following command will run the main snakemake workflow with 1 job submitted at a time, a wait period of 30 seconds for 'missing' files to be generated, and will use conda (necessary for the rule `gtc_to_vcf`). In this version of the workflow, the rules `idat_to_gtc` and `gtc_to_vcf` take in a *directory* of samples and parallelizes by threading, so at this time there is no benefit to having snakemake submit more than 1 job at a time.

```shell
mamba activate snakemake
snakemake -j 1 \
 --latency-wait 30 \
 --use-conda
```
## Run Snakemake on Remote Cluster
This workflow can also be run on a cluster system. An example of a job script for submitting this workflow to a Sun Grid Engine (SGE) scheduler can be found in `submit_job.sh`. Using `nohup` before the snakemake command will ensure that a job will finish even if the ssh-connection is dropped. If the ssh-connection is dropped or the process that snakemake was running on reaches a time or memory limit and is terminated before all jobs in the workflow are completed, you can restart snakemake by unlocking the working directory (`snakemake --unlock`) then resubmitting the job script with `--rerun-incomplete` appended to the snakemake command.

As an additional note, make sure that you have specified parallelization by threading and the same number of threads in your cluster scheduling command as is 

# Planned Updates
 * The workflow will be updated so that the conversion steps work by *sample* and not by *directory* so that parallelization can be performed across cores vs threads.
