# How to Run MitoDate

## Quick Start

Copy the configuration file from `assets/custom.config` and fill in your
parameters

## Run the Pipeline

Nextflow and Apptainer should be available in your path. Run with:

```bash
nextflow run -c custom.config -profile dardel main.nf
```

Or with your own HPC profile:

```bash
nextflow run -c custom.config -profile <your_hpc_profile> main.nf
```

## Input Data Formats

### Sample Metadata (TSV)

Tab-separated file with columns: `Sample_ID`, `Species`, `Origin`,
`Group-By`, `Calibrated_yBP` (or 'ND' for undated).

```
Sample_ID   Species                Origin      Group-By      Calibrated_yBP
YUK001      Mammuthus primigenius  Siberia     Mammoth       42000
YUK002      Mammuthus primigenius  Siberia     Mammoth       38000
```

### FASTA Sequences

Multiple sequence alignment with all sequences at same length:

```fasta
>YUK001
ATCGATCGATCGATCGATCGATCGATCGATCG
>YUK002
ATCGATCGATCGATCGATCGATCGATCGATCG
```

### Priors (CSV)

Parameter prior distributions:

```
Parameter,Distribution,Mean,StdDev,Lower,Upper
root_age,normal,50000,5000,0,100000
```

## Rerun Specific Samples

Rerun tip-dating BEAST analysis for specific samples without reprocessing:

```bash
nextflow run -c custom.config -profile dardel main.nf \
    --rerun_tip_dating "sample1,sample2,sample3"
```

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
nextflow run -c custom.config -profile your_hpc main.nf
```

---
