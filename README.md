# Evaluation of MinHash Segment Length Distribution
This repository contains the evaluation workflow used for the Segment Length Distribution chapter of my PhD Thesis.
The generated plots are stored in the folder `results/plots`.

## Software Requirements
The pipelines require a snakemake version of 5.10.0 or above and an installation of the [conda](https://docs.conda.io/en/latest/miniconda.html) package manager to run (we recommend using the miniconda Python 3.7 installer).
Lower versions result in a crash due to a `WorkflowError`.
Additionally, `pandas>=1.0.0` is required.
Assuming conda is already installed (for installation instructions please refer to the [(mini)conda](https://docs.conda.io/en/latest/miniconda.html) website), the required version of snakemake can be installed into a conda environment as follows:
```bash
$ conda create -n snakemake-env "snakemake>=5.10.0" "pandas>=1.0.0" -c bioconda -c conda-forge
```
After installation, activate the environment using:
```bash
$ conda activate snakemake-env
(snakemake-env)$ snakemake --version
5.10.0
```

## System Requirements
This workflow requires between 10 and 20GB of free disk space, preferably on a fast storage medium like an SSD.

It will automatically download reference genomes amounting to about 5 GB of zipped `.fasta.gz` files.
The download paths for these files and their local file system paths are specified in the file `workflow/genomes.tsv`.

The prediction of distribution requires both a large amount of RAM and can take up several hours.
We recommend assigning a small number of jobs to this workflow to judge its impact and subsequently increase the amount to fit your available hardware.

## Configuration
The workflow can be customized using the configuration file `workflow/config.yaml`.
In this file, the parameter `tmp_prefix` should be set to point to a storage location for simulated genomes and intermediary files.
This should ideally be fast, local storage, like an SSD.
Per default, a folder named `intermediate_results` is created in the workflow folder
Additionally, the config file can be used to adjust the parameters of both the genome simulation and analysis.
For example the number of simulated genomes (defined by the parameter `nr` in the `synthetic_genomes` section), can be reduced to lower the storage requirements of the workflow.
Finally, the config contains the path of the genomes file.

In the genomes file, the download location for the reference genomes is specified.
Unchanged, these will be downloaded into folder named `genomes` located in the the workflow folder.

## Usage

The workflow is started with the following command, assuming the snakemake environment has already been activated:

```bash
(snakemake-env)$ snakemake --use-conda --jobs <number of cores>
```

As mentioned above, select a number of jobs that can be accommodated by your hardware.
