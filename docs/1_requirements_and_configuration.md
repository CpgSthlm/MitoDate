# Requirements

## Input Data

MitoDate accepts multiple sequence alignments in FASTA format.
Metadata and prior information are provided in separate files.
Compressed files (gzip or bzip2) are supported and will be
automatically uncompressed.

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
A template is available at `assets/template_meta.tsv`.

The metadata file must contain the following columns:

- `Sample_ID`: Unique identifier for each sample
- `Species`: Species name
- `Origin`: Geographic origin (optional)
- `Group-By`: Taxonomic grouping for analysis
- `Calibrated_yBP`: Age in years before present, or `ND` if undated
- `TipDating`: Indicate `Y` if the sample is used for tip dating, or `N` if not.

Example format:
```
Sample_ID	Species	Origin	Group-By	Calibrated_yBP	TipDating
Sample1	M.primigenius	WEU	Clade1	10000	N
Sample2	M.primigenius	NNA	Clade2	ND	Y
SAmple3	M.primigenius	MEU	Clade3	ND	Y
```


## FASTA Sequences

A multiple sequence alignment in FASTA format with all
sequences aligned to the same length. Valid characters are
A, T, G, C, N and IUPAC degenerate DNA codes. Sequences must
have simple headers (e.g., `>Sample_ID`) and the header must
correspond exactly to the sample ID in the metadata.

## Priors File

A CSV (comma-separated values) file specifying the prior distributions used for tip dating individual samples.
A template is available at `assets/template_priors.csv`.

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
- `offset`: Offset applied to the distribution (default: 0


Example format:
```
taxa,date,prior,param1,param2,offset
Sample1,40000,uniform,1000,1000000,0
Sample2,45000,uniform,1000,1000000,0
```

> **Note:** Taxa names must exactly match the `Sample_ID` in the metadata file and the sequence headers in the FASTA alignment.

> **Note:** Tip-dating priors are the initial assumptions of molecular ages,
> or time constraints, applied to the tips of the phylogeny. They play a
> crucial role in the analysis by incorporating existing knowledge and
> uncertainty into the Bayesian inference. Such priors can be derived from
> archaeological, geographical, or stratigraphic context.


## Configuration File

The configuration file is a text file in Groovy format that
contains pipeline workflow steps and analysis parameters. A
template configuration file is available at `assets/custom.config`.
It is recommended to make a copy and modify it according to
your needs.

Start by setting `run_mode` to one of `tip_dating`, `rerun_samples`,
or `joint_tree`. The remaining required fields depend on the selected
mode.

See [How to Run MitoDate](3_how_to_run.md) for detailed parameter descriptions
and use cases.

## Compute Resources Configuration

Configuration files specifying compute resources for MitoDate
processes are located in the `configs/` directory. Different
profiles are provided for various compute environments:

- `configs/dardel.config` - For the Dardel HPC cluster (PDC)
- `configs/mjolnir_copenhagen.config` - For the Mjolnir HPC cluster
- `configs/modules.config` - Module-level configuration
- `configs/nf-core-label.config` - nf-core standard resource labels

You can adjust the resources allocated to each process in
these configuration files. If a process fails due to
insufficient memory or time, update the relevant
configuration to allocate more resources.


To add a new compute resource configuration file, place it in
the `configs/` directory and include it in the configuration:

```
includeConfig 'configs/my_compute_environment.config'
```
