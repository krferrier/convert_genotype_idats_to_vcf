# Overview
This snakemake workflow will convert Illumina genotyping .idat files to .gtc, then to .vcf format. This workflow requires that you have the Illumina IAAP-CLI software installed, the appropriate Illumina manifest files (.bpm and .csv format) and clusterfile for the genotyping array downloaded from the illumina website, and a reference genome in fasta format. You will also need to have Conda or Mamba installed. Instructions for installing the Illumina IAAP-CLI software and reference genome GRCh38/hg38 were modified from this gtc2vcf [README.md](https://github.com/freeseek/gtc2vcf/blob/1898320dab37c9da4e355f0fa31d2ab28d3632d5/README.md#installation). 

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
### Dependencies

# Using this Workflow
## Install snakemake 
## Fix for Circular Dependency
## Run Snakemake Locally
## Run Snakemake on Remote Cluster
# Planned Updates
