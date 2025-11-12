# SatTCR pipeline

[Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for assembling TCR data, and perform a saturation analysis.

More details available on:

<https://ong-research.github.io/SatTCR/>

## Instructions

To download from the command line using the command:

```sh
git clone git@github.com:Ong-Research/SatTCR.git
```

The SatTCR pipeline requires:

1. Docker: <https://www.docker.com/>
2. Snakemake: <https://snakemake.readthedocs.io/en/stable/>

It uses Snakemake to schedule the jobs to run the pipeline, and every job is run in a different container. 

### Getting Docker images

To pull Docker images that are going to be utilized by the pipeline, using the following commands:

```sh
cd SatTCR

docker pull staphb/fastqc # FastQC image
docker pull  staphb/multiqc # MultiQC image
docker pull staphb/trimmomatic # Trimmomatic image
docker build -t tcr/sat - < Dockerfile # R and Quarto image
docker pull ghcr.io/milaboratory/mixcr/mixcr:latest # MIXCR image
```

For MIXCR to work, it is necessary to get a license from <https://mixcr.com/mixcr/getting-started/milm/> and save it into a file.

### Configuring the SatTCR pipeline

Create a comma-separated value (csv) with 2 columns:

- `sample_name` : The name of the sample
- `sample_file`: The prefix of the files until before the `_R1` and `_R2` parts, e.g. if the pair of RNA-seq files are data/sample1_R1_L001.fastq.gz and data/sample1_R2_L001.fastq.gz, then this column is `data/sample1`.

Edit the `config/config.yaml` file. This file is divided by pieces in order to easily configure running the pipeline:

General configuration parameters:
- `threads`: Max. # of parallel threads used per process.
- `samplefile`: Location of the file with the samples.
- `seed`: Seed number for random number generation and sequence sampling during saturation analysis.
- `run_*`: Logical indicators to determine if running a stage of the pipeline
- `suffix`: This is regarding to the `samplefile`. If the pair of RNA-seq files are `data/sample1_R1_L001.fastq.gz` and `data/sample1_R2_L001.fastq.gz`. The suffix would be the remaining part after the R1/R2 parts, i.e. `_L001.fastq.gz`.
	
Docker configuration parameters:
- `run_line`: This is the docker command used to run every rule.
- `fastqc`, `multiqc`, `trimmomatic`, `rquarto` and `mixcr` are the names of the images that were pulled before.

In general, it is not necessary to modify these parameters unless a different image name is used or a specific need to configure how docker runs in the user’s system.

Trimmomatic configuration parameters:
- `trimmer`: A vector with the `trimmomatic` configuration to use. More information is available in <http://www.usadellab.org/cms/?page=trimmomatic>. But the general idea is to remove the low-quality nucleotides at the end of the sequences, or very short sequences.

MIXCR configuration parameters:
- `params`: The configuration line used to control MIXCR behavior. We used the line below to assemble the clonotypes analyzed used in this manuscript `rna-seq –species dog -b imgt.202214-2 –rna`. MIXCR provides a comprehensive list of preset configuration in <https://mixcr.com/mixcr/reference/overview-built-in-presets/>.
- `license_file`: Location of the file with the license. The pipeline uses this file to run MIXCR in a docker container.

Saturation configuration parameters:
- `samples`: A vector with the sample keys for which the saturation analysis is going to be processed
- `block_size` or `nblocks`: Either the # of sequences that are going to be sampled by block or the # of blocks of sequences used to split the original sequence files. 
- `bootstrap_replicates`: The number of times that the block bootstrap sampling procedure is going to be repeated.
This rule is computationally intensive, because in total there are going to be sampled `n_blocks-1 x n_boot_reps` pairs of sequence files and then MIXCR is used for each pair of files.

### Running the pipeline


In the instructions below, the flag `-c{k}` stands for running the rule with `{k}` parallel threads.

- Quality control: `snakemake -c{k} qc`. The output of this rule are an html report generated with MultiQC and quality profiles generated with the R package dada2 (Callahan et al. 2016). Either one of these analyses will depict quality score summaries at each position of the sequence files.
- Trim sequences: `snakemake -c{k} trim`. The output of this rule are the trimmed versions for every raw sequence file.
- Clonotype assembly with MIXCR: `snakemake -c{k} mixcr`. The output of this rule is a tsv file according to the AIRR format (<https://docs.airr-community.org/en/stable/datarep/overview.html>) for every set of RNA-seq paired files.
- Block bootstrap sampling: `snakemake -c{k} saturation`. This rule generates `n_blocks-1 x n_boot_reps` pairs of compressed fastq files.
- Saturation analysis: `snakemake -c{k} saturation`
- Generate the report: `snakemake -c{k} report`. This rule produces an html report compiled by `quarto` summarizing the results of the analysis.

## References

- [Bolotin DA, Poslavsky S, Mitrophanov I, Shugay M, Mamedov IZ, Putintseva EV, Chudakov DM. _"MiXCR: software for comprehensive adaptive immunity profiling_". Nature methods (2015)](https://www.nature.com/articles/nmeth.3364)

- [Bolotin DA, Poslavsky S, Davydov AN, Frenkel FE, Fanchi L, Zolotareva OI, Hemmers S, Putintseva EV, Obraztsova AS, Shugay M, Ataullakhanov RI, Rudensky AY, Schumacher TN, Chudakov DM. _"Antigen receptor repertoire profiling from RNA-seq data"_. Nature Biotechnology 35, (2017)](https://www.nature.com/articles/nbt.3979)

- [Greiff V, Miho E, Menzel U, Reddy ST. _"Bioinformatic and Statistical Analysis of Adaptive Immune Repertoires"_. Trend in Immunology (2015)](https://www.sciencedirect.com/science/article/abs/pii/S1471490615002239)
