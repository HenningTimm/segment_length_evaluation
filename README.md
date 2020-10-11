# segment_length_evaluation
This repository contains the evalutation workflow used for the Segment Length Distribution chapter of my PhD Thesis.


## Requirements
The pipelines require a snakemake version of 5.10.0 or above and an installation of the [conda](https://docs.conda.io/en/latest/miniconda.html) package manager to run (we recommend using the miniconda Python 3.7 installer).
Lower versions result in a crash due to a `WorkflowError`.
Assuming conda is already installed (for installation instructions please refer to the [(mini)conda](https://docs.conda.io/en/latest/miniconda.html) website), the required version of snakemake can be installed into a conda environment as follows:
```bash
$ conda create -n snakemake-env "snakemake>=5.10.0" -c bioconda -c conda-forge
```
After installation, activate the environment using:
```bash
$ conda activate snakemake-env
(snakemake-env)$ snakemake --version
5.10.0
```
