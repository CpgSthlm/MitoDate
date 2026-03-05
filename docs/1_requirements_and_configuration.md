# Requirements

## Input Data

MitoDate accepts multiple sequence alignments in FASTA format
(DNA sequence data). Metadata and prior information are provided
in separate files. Compressed files (gzip or bzip2) are supported
and will be automatically uncompressed.

Input data must be organized as described in
[Sample Metadata](1.-Pipeline-requirements-and-configuration#sample-metadata),
[FASTA Sequences](1.-Pipeline-requirements-and-configuration#fasta-sequences),
and [Priors File](1.-Pipeline-requirements-and-configuration#priors-file).

## Computational Requirements

Computational requirements depend on dataset size (number of
sequences), alignment length, and MCMC chain length. Small
datasets (e.g., 5-10 sequences with 10-20 kb alignments) can
be processed on a single machine with 8 cores and 16 GB memory
in a few hours. Larger datasets require more computational
resources and benefit significantly from multi-threading and
hardware acceleration (BEAGLE library).

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
see [Data Cleanup](4.-data-cleanup).


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


## FASTA Sequences

A multiple sequence alignment in FASTA format with all
sequences aligned to the same length. Valid characters are
A, T, G, C, N and IUPAC degenerate DNA codes. Sequences must
have simple headers (e.g., `>Sample_ID`).

## Priors File

A CSV (comma-separated values) file defining prior distributions
for analysis parameters (e.g., root divergence time, clock rate).

Example format:

```
Parameter,Distribution,Mean,StdDev,Lower,Upper
root_age,normal,50000,5000,0,100000
clock_rate,lognormal,0.001,0.0005,0.0001,0.01
```

## Configuration File

The configuration file is a text file in Groovy format that
contains pipeline workflow steps and analysis parameters. A
template configuration file is available at `assets/custom.config`.
It is recommended to make a copy and modify it according to
your needs.


See [How to Run MitoDate](3.-How-to-run-MitoDate) for detailed
parameter descriptions and use cases.

## Compute Resources Configuration

Configuration files specifying compute resources for MitoDate
processes are located in the `configs/` directory. Different
profiles are provided for various compute environments:

- `configs/dardel.config` - For the Dardel HPC cluster
- `configs/modules.config` - Standard module-level configuration
- `configs/nf-core-defaults.config` - nf-core standard settings

You can adjust the resources allocated to each process in
these configuration files. If a process fails due to
insufficient memory or time, update the relevant configuration
to allocate more resources.


To add a new compute resource configuration file, place it in
the `configs/` directory and include it in the configuration:

```
includeConfig 'configs/my_compute_environment.config'
```
