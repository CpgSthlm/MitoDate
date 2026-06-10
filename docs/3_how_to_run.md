# How to Run MitoDate

## Quick Start

Copy `assets/custom.config`, fill in your paths, choose a `run_mode`,
and then run the pipeline with your preferred profile.

## Run the Pipeline

Basic command:

```bash
nextflow run main.nf -c custom.config -profile dardel
```

Or with your own profile:

```bash
nextflow run main.nf -c custom.config -profile <your_hpc_profile>
```

## Choose a Mode

Set `run_mode` in `custom.config` to one of these values:

- `tip_dating`: run the full tip-dating workflow.
- `rerun_samples`: rerun selected samples from an earlier `tip_dating` run.
- `joint_tree`: build the final dated tree from a previous age summary.

## Mode 1: `tip_dating`

Use this mode for a normal run from aligned FASTA, metadata, and priors.

Set at least these fields in `custom.config`:

```groovy
params {
    run_mode = 'tip_dating'
    sample_metadata = '/path/to/metadata.tsv'
    fasta = '/path/to/alignment.fasta'
    priors = '/path/to/priors.csv'
}
```

## Mode 2: `rerun_samples`

Use this mode when you already have FASTA files from an earlier `tip_dating`
run and only want to rerun a few samples.

Keep the same `outdir` as the earlier run and set `rerun_samples` to a
comma-separated list of sample IDs.

```groovy
params {
    run_mode = 'rerun_samples'
    sample_metadata = '/path/to/metadata.tsv'
    priors = '/path/to/priors.csv'
    rerun_samples = 'sample1,sample2,sample3'
    outdir = 'results'
}
```

## Mode 3: `joint_tree`

Use this mode after reviewing the tip-dating results and choosing the
ages you want to carry into the final tree.

Set these fields in `custom.config`:

```groovy
params {
    run_mode = 'joint_tree'
    sample_metadata = '/path/to/metadata.tsv'
    fasta = '/path/to/alignment.fasta'
    age_summary_file = 'results/05_age_summary/run_<timestamp>/Results_ageSummary.csv'
}
```

## Notes

- Input file format details are described in [Requirements](1_requirements_and_configuration.md).
- `rerun_samples` expects FASTA files to already exist in `outdir/01_fastas`.
- `joint_tree` expects `age_summary_file` from a previous `tip_dating` run.

## Error Strategy

Set `process_errorStrategy` in `custom.config` to control how Nextflow reacts when a process fails.

- `ignore` (default): keep running and continue with remaining tasks.
- `terminate`: stop the workflow immediately on first failure.
- `finish`: stop scheduling new tasks, but let running tasks complete.
- `retry`: retry failed tasks (up to process/profile retry limits).


## HPC Profiles

To create a profile for your HPC system, copy `configs/dardel.config`
and modify it with your HPC's CPU, memory, time, and partition settings:

```groovy
process {
    executor = 'slurm'
    cpus = 16
    memory = '64 GB'
    time = '48h'
    queue = 'your_partition_name'
}
```

Then run with your profile:

```bash
nextflow run main.nf -c custom.config -profile your_hpc
```

