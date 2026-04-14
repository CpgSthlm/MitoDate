# Pipeline Steps

MitoDate consists of several workflow steps that can be enabled or
disabled in the config file. This section describes these workflow
steps in more detail.

Following are the main workflow steps of the pipeline:

- FASTA Processing
- XML Generation
- BEAST Analysis
- Result Parsing

## FASTA Processing

This workflow performs preprocessing on FASTA sequence alignments:

**Multiple Sample Dating** (`single_sample_dating = false`):
Reads the multiple sequence alignment and renames sequences based on
metadata information. New sequence headers include Sample ID, Species,
Origin, Group, and Age. Output is a single renamed FASTA file for
all samples.

**Single Sample Dating** (`single_sample_dating = true`):
Extracts only sequences with calibration dates (excluding samples
marked 'ND' as undated). Creates individual FASTA alignment files
for each dated sample, allowing independent analysis of each sample.

## XML Generation

This workflow generates BEAST configuration files that define the
complete phylogenetic model. The pipeline creates XML files containing:

- Sequence data and taxon information
- Taxon groups and age calibrations (from metadata)
- DNA substitution model (HKY or GTR)
- Molecular clock model (strict)
- Population model (SkyGrid)
- Prior distributions for parameters
- MCMC settings (chain length, logging frequency)

One XML file is generated per sample or partition combination. These
XML files serve as the complete specification for BEAST analysis.

## BEAST Analysis

This workflow performs Bayesian phylogenetic inference using
BEAST 1.10.4. For more details of BEAST, please check:
https://beast.community/about

- Runs MCMC sampling to explore parameter space
- Estimates divergence times (root age) for samples
- Estimates substitution rates and other evolutionary parameters
- Uses hardware acceleration (BEAGLE library) for speed
- Multi-threading support scales to available CPU cores

Output files:
- `*.log` - Parameter trace file with MCMC samples
- `*.trees` - Sampled phylogenetic trees in Newick format

**Computation time** depends on dataset size and chain length. For example, if 100 million iterations are used:

- Small datasets (5-20 sequences): 1-6 hours
- Medium datasets (50-100 sequences): 12-48 hours
- Large datasets (200+ sequences): 48+ hours

## Result Parsing

This workflow extracts and summarizes results from BEAST analysis:

**BEAST Log Parser**:
Reads the parameter trace file (`.log`) from each BEAST run and
calculates summary statistics including mean, and 95%
highest posterior density (HPD) intervals for all parameters. Computes effective sample
size (ESS) to assess convergence quality.

**Result Collection**:
Aggregates parsed results across all samples into a single results
table (`combined_results.csv`). Combines BEAST estimates with original
metadata for easy interpretation.

**Output** (`combined_results.csv`):
- Sample ID
- Effective sample size (ESS) for `joint`, `prior`, and `likelihood`
- Root age estimates (mean and ESS)
- Sample age estimates (mean, stdev, 95% HPD, and ESS)
- Other parameter estimates

---

**See Also**:
- [How to Run MitoDate](3.-How-to-run-MitoDate) - Execution guide
- [Requirements and Configuration](1.-Pipeline-requirements-and-configuration) - Configuration options
