# Requirements

## Input Data

MitoDate accepts multiple sequence alignments in FASTA format
(aligned DNA sequence data). Metadata and prior information
are provided in separate files. Compressed files (gzip or
bzip2) are supported and will be automatically uncompressed.

Input data must be organized as described in
[Sample Metadata](1_requirements_and_configuration.md#sample-metadata),
[FASTA Sequences](1_requirements_and_configuration.md#fasta-sequences),
[Priors File](1_requirements_and_configuration.md#priors-file).

## Computational Requirements

Computational requirements depend on dataset size (number of
sequences), alignment length, and MCMC chain length. Small
datasets (e.g., 10-20 sequences with 10-20 kb alignments) can
be processed on a single machine with 8 cores and 16 GB
memory in a few hours. Larger datasets require more
computational resources and benefit significantly from
multi-threading and hardware acceleration (BEAGLE library).

Nextflow must be available in your environment (install via
Conda: `conda env create -f environment.yaml`). The pipeline
runs software dependencies in containers (highly recommended)
or with Conda/Mamba. The pipeline supports Docker, Apptainer,
and Singularity containers. By default, container images are
downloaded and stored in the `work` directory. This can be
changed by setting environment variables:

- Apptainer: `export NXF_APPTAINER_CACHEDIR=/path/to/cache`
- Singularity: `export NXF_SINGULARITY_CACHEDIR=/path/to/cache`
- Docker: `export NXF_DOCKER_CACHEDIR=/path/to/cache`

Using a terminal multiplexer (e.g., tmux or screen) is
highly recommended to keep Nextflow running in the background.

Ensure sufficient storage space, as the pipeline generates
many large intermediate files in the `work` directory.
Temporary files can be deleted after successful completion,
see [Data Cleanup](4_data_cleanup.md).


# Configuration

## Sample Metadata

The sample metadata file is a TSV (tab-separated values) file
that provides information about each sample to be analyzed.
A template is available at `assets/README.md`.

The metadata file must contain the following columns:

- `Sample_ID`: Unique identifier for each sample
- `Species`: Species name
- `Origin`: Geographic origin (optional)
- `Group-By`: Taxonomic grouping for analysis
- `Calibrated_yBP`: Age in years before present, or `ND` if undated
- `TipDating`: Indicate `Y` if the sample is used for tip dating, or `N` if not.


## FASTA Sequences

A multiple sequence alignment in FASTA format with all
sequences aligned to the same length. Valid characters are
A, T, G, C, N and IUPAC degenerate DNA codes. Sequences must
have simple headers (e.g., `>Sample_ID`) and the header must
correspond exactly to the sample ID in the metadata.

## Priors File

A CSV (comma-separated values) file specifying the prior distributions used for tip dating individual samples.

Example format:
```
taxa,date,prior,param1,param2,offset
YUK001,40000,uniform,1000,1000000,0
YUK002,45000,uniform,1000,1000000,0
```

Column descriptions:
- `taxa`: Sample ID
- `date`: Approximate age of the sample in years before present used to initialize the Bayesian chain.
- `prior`: Type of prior distribution.
    - Options: `uniform`, `lognormal`, or `normal`.
- `param1`: First parameter of the prior distribution.
  - `uniform`: lower bound
  - `lognormal` / `normal`: mean (`mu`)
- `param2`: Second parameter of the prior distribution.
  - `uniform`: upper bound
  - `lognormal` / `normal`: standard deviation (`sigma`)
- `offset`: Offset applied to the distribution (default: 0)

> **Note:** Taxa names must exactly match the `Sample_ID` in the metadata file and the sequence headers in the FASTA alignment.


## Configuration File

The configuration file is a text file in Groovy format that
contains pipeline workflow steps and analysis parameters. A
template configuration file is available at `assets/custom.config`.
It is recommended to make a copy and modify it according to
your needs.

See [How to Run MitoDate](3_how_to_run.md) for detailed parameter descriptions
and use cases.

## Compute Resources Configuration

Configuration files specifying compute resources for MitoDate
processes are located in the `configs/` directory. Different
profiles are provided for various compute environments:

- `configs/dardel.config` - For the Dardel HPC cluster
- `configs/modules.config` - Standard module-level configuration
- `configs/nf-core-defaults.config` - nf-core standard settings

You can adjust the resources allocated to each process in
these configuration files. If a process fails due to
insufficient memory or time, update the relevant
configuration to allocate more resources.


To add a new compute resource configuration file, place it in
the `configs/` directory and include it in the configuration:

```
includeConfig 'configs/my_compute_environment.config'
```
