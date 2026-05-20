# Pipeline Steps

For a description of the three workflow modes and how to run them, see [How to Run MitoDate](3_how_to_run.md).

The `tip_dating` and `rerun_samples` modes run these steps:

- FASTA Processing
- XML Generation
- BEAST Analysis
- Result Parsing

The `joint_tree` mode runs a separate set of steps described at the [bottom of this page](#joint-tree-workflow).

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

One XML file is generated per sample or partition combination.

## BEAST Analysis

This workflow performs Bayesian phylogenetic inference using
BEAST 1.10.4.

- Runs MCMC sampling to explore parameter space
- Estimates divergence times (root age) for samples
- Estimates substitution rates and other evolutionary parameters
- Uses hardware acceleration (BEAGLE library) for speed
- Multi-threading support scales to available CPU cores

Output files:
- `*.log` - Parameter trace file with MCMC samples
- `*.trees` - Sampled phylogenetic trees in Newick format

## Result Parsing

This workflow extracts and summarizes results from BEAST analysis:

**BEAST Log Parser**:
Reads the parameter trace file (`.log`) from each BEAST run and
calculates summary statistics, including mean and 95% highest posterior
density (HPD) intervals. It also computes effective sample size (ESS)
to help assess convergence.

**Result Collection**:
Aggregates parsed results across all samples into a single results
table and combines BEAST estimates with the original metadata.

**Output:**
- Sample ID
- Effective sample size (ESS) for `joint`, `prior`, and `likelihood`
- Root age estimates (mean and ESS)
- Sample age estimates (mean, stdev, 95% HPD, and ESS)
- Other parameter estimates

---

## Joint Tree Workflow

The `joint_tree` mode builds a final date-calibrated tree after you have
reviewed the tip-dating results from a previous `tip_dating` run.

**Merge dates into metadata**: Reads the age summary file from the earlier
run and merges the estimated sample ages into the metadata. The FASTA
alignment is relabelled accordingly.

**XML generation**: Generates a BEAST XML for the full alignment. No priors
file is needed here because all samples are now treated as dated.

**BEAST analysis**: Runs BEAST on the joint XML to produce a dated tree.

**Tree annotation**: Summarises the posterior tree sample into a single
Maximum Clade Credibility (MCC) tree using TreeAnnotator.

**Tree visualisation**: Plots the annotated MCC tree coloured by the
`Group-By` column from the metadata.

See also:
- [How to Run MitoDate](3_how_to_run.md)
- [Requirements and Configuration](1_requirements_and_configuration.md)
