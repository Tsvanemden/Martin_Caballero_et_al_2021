# Snakemake pipeline for NGS data analysis in Martín Caballero et al., 2021

The snakemake pipeline will process the NGS data of Martín Caballero et al., 2021
In short, STAR will map the reads to the transriptome and RSEM will further process them.
With the enclosed R files the NGS related figures of Martín Caballero et al., 2021 can be reproduced,
although authors may choose to plot the data using a different software. 

## Setup

* clone this directory.

* create a conda environment specified in sra_star_rsem.yaml.

```
conda env create --name caballero_et_al_2021 --file caballero_et_al_2021.yaml
```

*activate the conda environment

```
conda activate caballero_et_al_2021
```


## Downlaod data from GEO

* navigate to this directory. 

* download the data and store it in data/sraFiles directory

```
prefetch <accession>
fastq-dump <accession> -O data/sraFiles/
```

### on local machine
with 16 cores

```
snakemake -np
snakemake --cores 16
```

## on the cluster
with many cores

```
snakemake -j 999 --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -n {cluster.n} --mem {cluster.mem}"
```


## Downlaod data from GEO

* Run the figures_Caballero_et_al_2021.R script in this directory to produce the figures

