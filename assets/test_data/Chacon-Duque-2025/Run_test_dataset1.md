# Test Dataset: Chacón-Duque et al. (2025)

## About the Dataset

This tutorial uses a test dataset derived from the mitochondrial genomes analyzed in [Chacón-Duque et al., *Molecular Biology and Evolution* (2025)](https://doi.org/10.1093/molbev/msaf065).

The example dataset contains 10 representative undated mitogenomes, together with radiocarbon-dated sequences and outgroups used for molecular dating.

Three input files are provided:

1. **Multiple sequence alignment (FASTA):** An alignment containing outgroups, radiocarbon-dated mitogenomes, and 10 undated mitogenomes, with the variable number tandem repeat (VNTR) region removed.
2. **Sample metadata (TSV):** Sample information, including taxonomic assignments, geographic origins, and calibration ages.
3. **Tip-dating priors (CSV):** Prior distributions specifying the expected ages and associated uncertainties of undated samples.

In the following sections, we will walk through the preparation of these files, configure the pipeline step by step, and run a test analysis using MitoDate.

## 1. Set Up Your Working Environment

MitoDate is built using [Nextflow](https://www.nextflow.io/), a workflow management system designed for reproducible and scalable scientific analyses.

If you are new to Nextflow, we recommend following the [Nextflow installation guide](https://nf-co.re/docs/get_started/environment_setup/nextflow).

### 1.1 Install Nextflow

You can install Nextflow using Conda:

```bash
conda create --name nf-env bioconda::nextflow
```

Activate the environment:

```bash
conda activate nf-env
```

Verify that Nextflow is installed correctly:

```bash
nextflow info
```

Alternatively, if you are working on an HPC cluster where Nextflow is managed through an environment module system, you can load it directly:

```bash
module load nextflow
```

Note that module names and versions may vary across HPC systems.

### 1.2 Set Up a Container Engine

MitoDate uses containerized software environments to manage its dependencies. Depending on your computing environment, you can use Docker, Singularity, or Apptainer.

For additional information, see the [nf-core software dependencies guide](https://nf-co.re/docs/get_started/environment_setup/software-dependencies).

**Docker**

[Docker](https://www.docker.com/) packages software and its dependencies into isolated containers, ensuring reproducibility across computing environments. It is commonly used on local workstations, personal servers, and cloud computing platforms.

Installation instructions are available in the [Docker documentation](https://docs.docker.com/get-started/get-docker/).

**Singularity**

[SingularityCE](https://sylabs.io/singularity/) is a container platform designed for high-performance computing (HPC) environments. Unlike Docker, it can run containers without requiring administrator privileges, making it particularly suitable for shared computing clusters.

Installation instructions are available in the [SingularityCE Quick Start guide](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html).

On some HPC systems, Singularity can be loaded using:

```bash
module load singularity
```

**Apptainer**

[Apptainer](https://apptainer.org/) is an open-source container platform derived from Singularity. It provides similar functionality and is widely used in HPC environments.

If Apptainer is available through your HPC module system, load it using:

```bash
module load apptainer
```

Only one supported container engine is required.

### 1.3 Download MitoDate

Clone the MitoDate repository from GitHub:

```bash
git clone https://github.com/CpgSthlm/MitoDate.git
```

You now have the pipeline code, configuration templates, and example datasets needed for this tutorial.

## 2. Prepare the Dataset

Although the test dataset has already been prepared and is available in the repository, this section explains how the input files were generated and how they should be formatted when preparing your own data.

### 2.1 Multiple Sequence Alignment (MSA)

The input multiple sequence alignment (MSA) should include all sequences used in the molecular dating analysis.

For this test dataset, the alignment includes:

- Radiocarbon-dated mammoth mitogenomes used for tip calibration.
- African and Asian elephant mitogenomes used as outgroups.
- Ten undated mammoth mitogenomes to be analyzed.

All sequences were aligned using [MAFFT](https://mafft.cbrc.jp/alignment/software/) and saved in FASTA format.


**For non-partitioned analyses:**

The VNTR region was manually removed from the alignment to avoid potential alignment ambiguities. The resulting FASTA file is used for non-partitioned analyses:

`ChaconDuque-2025-10-VNTRremoved.fasta`

**For partitioned analyses:**

The alignment can be divided into tRNA, rRNA, coding sequences (CDS), and control regions using a corresponding GFF annotation file with coordinates matching the alignment. Since the VNTR region is annotated in the GFF file, it does not need to be manually removed from the alignment.

For partitioned analyses, use the original FASTA file containing the VNTR region:

`ChaconDuque-2025-10.fasta`

Copy the appropriate alignment file to your working directory, depending on whether you are performing a partitioned or non-partitioned analysis.

For this tutorial, we use the non-partitioned alignment as an example:

```bash
cp /path/to/MitoDate/assets/test_data/Chacon-Duque-2025/ChaconDuque-2025-10-VNTRremoved.fasta /path/to/workdir/
```

Replace `/path/to/MitoDate` and `/path/to/workdir` with the corresponding paths on your system. The same convention is used throughout this tutorial.



### 2.2 Prepare Sample Metadata

The sample metadata file provides information about each sequence included in the MSA.

A metadata template is provided in the `assets` directory. You can use this template to prepare your own dataset.

For this tutorial, the metadata has already been prepared. Take a look at few rows shown below:

```tsv
Sample_ID	Species	Origin	Group-By	Calibrated_yBP	TipDating
P033	M.tro	NAS	Mammuthus	ND	Y
P035	M.tro	NAS	Mammuthus	ND	Y
SVEP43	M.col	CNA	Mammuthus	12891	N
SVEP5681	M.col	CNA	Mammuthus	22362	N
SVEP57	M.col	CNA	Mammuthus	12891	N
SVEP72	M.col	CNA	Mammuthus	12689	N
JN673263	L.cyc	NA	Loxodonta	0	N
JN673264	L.cyc	NA	Loxodonta	0	N
KJ557424	L.cyc	NA	Loxodonta	0	N
```

Each column is described below:

- **`Sample_ID`**: A unique identifier for each sample. This must match the corresponding sequence identifier in the MSA exactly.

- **`Species`**: The species name associated with each sample. Scientific names, common names, or abbreviated names can be used. This column is primarily intended for dataset documentation and does not directly affect the molecular dating analysis.

- **`Origin`**: The geographic origin of each sample. This column is used for documentation and does not directly influence the molecular dating analysis.

- **`Group-By`**: Defines how samples are grouped for downstream analyses, including estimation of the time to the most recent common ancestor (TMRCA). In this example, sequences are grouped by genus: *Mammuthus*, *Elephas*, and *Loxodonta*.

- **`Calibrated_yBP`**: The sample age in calibrated years before present (cal yBP). Use `0` for modern samples and `ND` for undated samples. This column accepts numeric ages or `ND` only.

- **`TipDating`**: Specifies whether a sample should be included in the tip-dating procedure. Use `Y` for samples to be dated and `N` for samples that should not be tip-dated.

For example, samples `P033`, and `P035` are marked as `Y`, indicating that their ages will be estimated through molecular dating. In contrast, the radiocarbon-dated samples are marked as `N`, as their calibrated ages are already provided.

Copy the metadata file to your working directory:

```bash
cp /path/to/MitoDate/assets/test_data/Chacon-Duque-2025/ChaconDuque2025-meta-10.tsv /path/to/workdir/
```

### 2.3 Prepare Tip-Dating Priors

The next step is to define age priors for the undated samples.

MitoDate currently supports three types of prior distributions:

- **Uniform:** Assumes equal prior probability across a specified age interval.
- **Normal:** Defines a distribution using a mean and standard deviation.
- **Lognormal:** Defines an asymmetric distribution using the parameters supported by the pipeline.

The choice of prior should reflect the available chronological information and the uncertainty associated with each sample.

For example, suppose a sample is estimated to be approximately 70,000 years old, with an uncertainty of 10,000 years. You could define a uniform prior spanning 60,000–80,000 years or a normal prior centered at 70,000 years with a standard deviation of 10,000 years.

However, informative normal or lognormal priors can strongly influence posterior age estimates. We therefore recommend starting with a reasonably broad uniform prior when precise chronological constraints are unavailable, followed by sensitivity analyses to assess the influence of prior choice.

**Input format**

The first few rows of the prepared prior file are shown below:

```csv
taxa,date,prior,param1,param2,offset
P033,1050000,uniform,1000,2000000,0
P035,1200000,uniform,1000,2000000,0
P037,1390000,uniform,1000,2000000,0
```

In this example, all samples are assigned uniform priors between 1,000 and 2,000,000 years before present.

You can also specify different prior distributions for individual samples. For example:

```csv
taxa,date,prior,param1,param2,offset
P033,1050000,uniform,1000,2000000,0
P035,1200000,normal,1200000,200000,50000
P037,1390000,lognormal,1400000,400000,50000
```

The columns are defined as follows:

- **`taxa`**: The sample identifier, which must match the `Sample_ID` in the metadata file.

- **`date`**: An approximate initial age (in years before present) used to initialize the Bayesian analysis. A reasonable estimate is sufficient, as this value is not the final inferred age.

- **`prior`**: The prior distribution assigned to the sample. Supported values are `uniform`, `normal`, and `lognormal`.

- **`param1`**: The first parameter of the prior distribution. For uniform priors, this is the lower bound. For normal and lognormal priors, this represents the mean parameter.

- **`param2`**: The second parameter of the prior distribution. For uniform priors, this is the upper bound. For normal and lognormal priors, this represents the standard deviation parameter.

- **`offset`**: An offset applied to the selected prior distribution. The default value is `0`.

Make sure that the prior table includes the samples selected for tip dating and that their identifiers match those in the metadata and MSA files.

Copy the prepared priors to your working directory:

```bash
cp /path/to/MitoDate/assets/test_data/Chacon-Duque-2025/ChaconDuque-2025-priors-10.csv /path/to/workdir/
```

## 3. Prepare the Configuration File

MitoDate uses a Nextflow configuration file to define the input files, analysis parameters, and execution settings.

A configuration template is provided in the repository. Copy it to your working directory:

```bash
cp /path/to/MitoDate/assets/custom.config /path/to/workdir/ChaconDuque_2025_10_withoutPartition_tipDating.config
```

We will use this configuration file throughout the tutorial.

The pipeline supports three run modes:

1. **`tip_dating`**: Estimates the ages of undated samples using molecular dating.
2. **`rerun_samples`**: Reruns selected samples from a previous tip-dating analysis.
3. **`joint_tree`**: Constructs a jointly dated phylogenetic tree using the reviewed tip-dating results.

In this tutorial, we will first run the `tip_dating` mode.

## 4. Mode 1 — Tip Dating

Open the configuration file in your preferred text editor:

We will now configure the parameters required for the test analysis.

### 4.1 Define the Workflow Name and Run Mode

First, assign a descriptive name to the workflow and select the `tip_dating` mode:

```groovy
workflow_run_name = 'ChaconDuque_2025_10_withoutPartition'
run_mode          = 'tip_dating'
```

The workflow name is used to identify output files generated by the analysis.

The parameters `rerun_samples` and `age_summary_file` are not required for this mode and can be left at their default values.

### 4.2 Specify the Input Files

Provide the absolute paths to the three prepared input files:

**Sample metadata**

```groovy
sample_metadata = '/path/to/workdir/ChaconDuque2025-meta-10-speciesGroup.tsv'
```

**Multiple sequence alignment**

```groovy
fasta = '/path/to/workdir/Chacon-Duque-etal-2025-10representatives-VNTRremoved.fasta'
```

**Tip-dating priors**

```groovy
priors = '/path/to/workdir/Chacon-Duque-2025-priors-10representatives.csv'
```

Ensure that all paths point to existing files and that the filenames match those copied into your working directory.

Other mode-specific input parameters can be left unchanged.

### 4.3 Configure Tip-Dating Parameters

**Single-sample dating**

Enable single-sample dating:

```groovy
single_sample_dating = 'true'
```

This setting enables the single-sample dating workflow, in the following analyses the undated samples will be dated individually.


**Sequence partitioning**

MitoDate supports both partitioned and non-partitioned analyses. Configure the `split_partition` parameter according to your analysis.

**For non-partitioned analyses:**

Disable sequence partitioning:

```groovy
split_partition = 'false'
```

The entire mitochondrial alignment, excluding the VNTR region, will be analyzed as a single partition.

**For partitioned analyses:**

Enable sequence partitioning:

```groovy
split_partition = 'true'
```

In this case, specify the path to the complete mitochondrial alignment, including the VNTR region:

```groovy
fasta = '/path/to/workdir/Chacon-Duque-2025-10.fasta'
```

You also need to provide the corresponding GFF annotation file, which defines the genomic regions used for partitioning:

```groovy
gff = '/path/to/ChaconDuque-2025-NC_007596.2.liftoff.gff3'
```

The pipeline will use the GFF annotations to divide the alignment into tRNA, rRNA, coding sequences (CDS), and control regions. The VNTR region is identified through the annotation, so manual removal is not required.

**Substitution model**

Specify the nucleotide substitution model using the `substitution_model` parameter.

**For non-partitioned analyses:**

A single substitution model is applied to the entire alignment:

```groovy
substitution_model = 'HKY+G+I'
```

Here, `HKY` refers to the Hasegawa–Kishino–Yano substitution model, `+G` accounts for gamma-distributed rate variation among sites, and `+I` allows for a proportion of invariant sites.

**For partitioned analyses:**

You can specify a separate substitution model for each partition using a comma-separated list of partition names and their corresponding models:

```groovy
substitution_model = 'tRNA:HKY+G+I, rRNA:HKY+G+I, CDS:HKY+G+I, D_loop:HKY+G'
```

In this example, the `HKY+G+I` model is applied to the tRNA, rRNA, and CDS partitions, while the `HKY+G` model is used for the control region (`D_loop`).

This allows different functional regions of the mitochondrial genome to be modeled separately, accounting for potential differences in their evolutionary patterns.



**Root calibration**

The root age is calibrated using an estimated divergence time of 5.3 million years between *Loxodonta* and the *Elephas–Mammuthus* lineage.

```groovy
root_mean   = '5300000'
root_stdev  = '500000'
root_offset = '0'
```

The prior mean is set to 5.3 million years, with a standard deviation of 500,000 years.

**Population model**

Specify the coalescent population model:

```groovy
population_model = 'skygrid'
```

The Bayesian Skygrid model allows changes in effective population size over time.

**Molecular clock model**

Select a strict molecular clock:

```groovy
clock_model = 'strict'
```

The strict clock assumes a constant substitution rate across all branches of the phylogenetic tree.

### 4.4 Configure MCMC Settings

Next, specify the settings for the BEAST Markov chain Monte Carlo (MCMC) analysis.

For this tutorial, we use a short chain of 1,000 iterations to demonstrate the workflow and reduce computational time:

```groovy
chain_length = '1000'
log_every    = '100'
burnin       = '10'
```

These settings correspond to:

- `chain_length`: 1,000 MCMC iterations.
- `log_every`: Log parameters and trees every 100 iterations.
- `burnin`: Discard the first 10% of samples when summarizing the results.

**Note:** These settings are intended for testing the pipeline rather than obtaining final molecular dating estimates. For research analyses, longer MCMC chains and appropriate convergence diagnostics are required.

### 4.5 Configure HPC Execution Profiles
Now we move to the Module- (tool-)specific parameter configurations. Leave the intermediate files store in `work` folder. And change the cache folder and HPC job ID if necessary. 

MitoDate provides predefined Nextflow execution profiles for different computing environments, including the Dardel and Mjolnir HPC clusters.

The reference cluster-specific configuration files in the `configs/` directory.

Examples include:

- `dardel`: Dardel HPC cluster.
- `mjolnir`: Mjolnir HPC cluster.
- `dardel_test`: Test configuration for the Dardel cluster.
- `mjolnir_test`: Test configuration for the Mjolnir cluster.


The `mjolnir` profile, for example, loads the corresponding `configs/mjolnir.config` file and specifies the cluster's SLURM and Singularity settings, where the user define how much computation source allocated to each job. 

If you are using one of the supported HPC systems, simply select the appropriate profile when launching the pipeline.

If you are running on another HPC system, you may need to create a custom execution profile with the appropriate scheduler, resources, and container settings.

For this test tutorial, we will use the `mjolnir_test` profile.

All remaining parameters can be left at their default values.

## 5. Run the Pipeline

We are now ready to launch the test analysis.

### 5.1 Start a Persistent Terminal Session

When running the pipeline on an HPC cluster, we recommend using `tmux` or `screen` to keep the Nextflow session active if your SSH connection is interrupted.

For more information, see the [tmux cheat sheet](https://tmuxcheatsheet.com/) or the [official tmux Getting Started guide](https://github.com/tmux/tmux/wiki/Getting-Started).

To create a new tmux session:

```bash
tmux new -s MitoDate-test
```

You can detach from the session by pressing `Ctrl+B`, followed by `D`.

To reconnect later:

```bash
tmux attach -t MitoDate-test
```

### 5.2 Load Required Modules

If your HPC environment uses a module system, load the required software before launching Nextflow.

For example:

```bash
module load openjdk/20.0.0
module load nextflow
module load singularity/3.8.7
```

The available module names and versions may differ depending on your HPC system.

### 5.3 Launch MitoDate

Navigate to your working directory:

```bash
cd /path/to/workdir
```

Launch the pipeline using the prepared configuration file:

```bash
nextflow run /path/to/MitoDate/main.nf \
    -c ChaconDuque_2025_10_withoutPartition_tipDating.config \
    -profile mjolnir_test
```

Here:

- `nextflow run` launches the MitoDate workflow.
- `-c` specifies the custom configuration file.
- `-profile mjolnir_test` selects the execution profile for the Mjolnir HPC cluster.

Nextflow will automatically submit and manage the computational processes according to the selected execution profile.

For local execution, choose an appropriate container profile, such as `docker`, `singularity`, or `apptainer`, instead of `mjolnir`.

Once the workflow has completed, the results will be available in the configured output directory (`results` by default).

 

### 5.4 Monitor the Workflow

After launching the pipeline, Nextflow will display the workflow execution status in your terminal. You should see output similar to the following:

```text
N E X T F L O W   ~  version 25.10.4

Launching `/path/to/MitoDate/main.nf` [suspicious_snyder] DSL2 - revision: df02260152

        Workflow run name: ChaconDuque_2025_10-withoutPartition
        Selected run mode: tip_dating

executor >  slurm (1)
[0b/9dc387] SPLITFASTA (Chacon-Duque-2025-10-VNTRremoved.fasta) [  0%] 0 of 1
[-        ] GENERATEXML                                         -
[-        ] BEAST_TIP_DATING                                    -
[-        ] BEAST_LOG_PARSER                                    -
[-        ] COLLECT_RESULTS                                     -
```

Nextflow displays the progress of each process in the pipeline. The percentage indicates how many tasks have been completed, while the remaining processes will start as their required inputs become available.

For this analysis, the workflow consists of five main processes:

1. **`SPLITFASTA`**: Prepares the sequence alignments for individual tip-dating analyses.
2. **`GENERATEXML`**: Generates BEAST XML configuration files.
3. **`BEAST_TIP_DATING`**: Runs the Bayesian molecular dating analyses.
4. **`BEAST_LOG_PARSER`**: Extracts and summarizes results from the BEAST output files.
5. **`COLLECT_RESULTS`**: Organizes the results into the final output directory.

Once all processes have finished successfully, the terminal output will look similar to this:

```text
N E X T F L O W   ~  version 25.10.4

Launching `/path/to/MitoDate/main.nf` [suspicious_snyder] DSL2 - revision: df02260152

        Workflow run name: ChaconDuque_2025_10-withoutPartition
        Selected run mode: tip_dating

executor >  slurm (32)
[0b/9dc387] SPLITFASTA (Chacon-Duque-2025-10-VNTRremoved.fasta)                    [100%] 1 of 1 ✔
[c2/7ce3e1] GENERATEXML (Chacon-Duque-2025-10-VNTRremoved.TipDating.P033.fasta)    [100%] 10 of 10 ✔
[7a/0ad985] BEAST_TIP_DATING (Chacon-Duque-2025-10-VNTRremoved.TipDating.P033.xml) [100%] 10 of 10 ✔
[a1/f872a2] BEAST_LOG_PARSER (Chacon-Duque-2025-10-VNTRremoved.TipDating.L168)     [100%] 10 of 10 ✔
[7c/bde864] COLLECT_RESULTS (yyyy-mm-dd_xx-xx-xx)                                  [100%] 1 of 1 ✔

            Thank you for using MitoDate.

            Results are located in the folder: results.

Completed at: dd-mm-yyyy xx:xx:xx
Duration    : 7m 22s
CPU hours   : 0.1
Succeeded   : 32
```

All processes have now completed successfully! The `Succeeded` field indicates the total number of successfully completed tasks.

**Note:** Successful workflow execution does not necessarily mean that the MCMC analyses have converged. We will examine the output and convergence diagnostics in the next section.

## 6. Explore the Results

Once the workflow is complete, the output files will be organized in the `results` directory.

You will find the following subdirectories:

```text
results/
├── 01_fastas/
├── 02_xml/
├── 03_beast/
├── 04_beast_log_parser/
└── 05_age_summary/
```

Each directory contains outputs from a different stage of the analysis:

- **`01_fastas/`**: Sequence alignments prepared for individual tip-dating analyses.
- **`02_xml/`**: BEAST XML configuration files generated for each analysis.
- **`03_beast/`**: Raw BEAST outputs, including MCMC log and tree files.
- **`04_beast_log_parser/`**: Processed results and statistics extracted from the BEAST logs.
- **`05_age_summary/`**: Summary tables of estimated sample ages and associated statistics.

You are welcome to explore these directories to examine the intermediate files. For now, we will focus on the age summary.

### 6.1 Examine the Age Summary

Navigate to the `05_age_summary` directory and open the generated summary CSV file.

The first few rows of the example output are shown below:

```csv
Sample_ID,joint_ESS,prior_ESS,likelihood_ESS,age(root)_mean,age(root)_ESS,clock.rate_mean,clock.rate_ESS,age(Sample)_mean,age(Sample)_Stdev,age(Sample)_HPD_95_Lower,age(Sample)_HPD_95_Upper,age(Sample)_ESS
P048,4.145447809476067,4.870433699124905,3.609577265589051,450012.2510041753,4.144399075580542,0.25766297840984065,4.871881199717449,450000.0,0.0,450000.0,450000.0,0.0
MD228,3.226556173029218,3.186428171700787,3.3898632881092734,450020.0624216007,3.851104577532026,0.5990207586965252,3.1831237424058503,450000.0,0.0,450000.0,450000.0,0.0
P035,4.217312928656525,4.835766825269838,4.004552141112484,450059.819471998,4.746497503776188,0.346435736225648,4.828098972825187,450000.0,0.0,450000.0,450000.0,0.0
P043,3.3517581798545266,3.363186006487618,3.707546990463755,450048.80335338024,4.035149338143429,0.6273081551734815,3.3639154365410677,450000.0,0.0,450000.0,450000.0,0.0
```

The summary contains several important statistics:

- **`Sample_ID`**: The identifier of the sample being dated.
- **`joint_ESS`**: Effective sample size (ESS) of the joint posterior probability.
- **`prior_ESS`**: ESS of the prior probability.
- **`likelihood_ESS`**: ESS of the likelihood.
- **`age(root)_mean`**: Posterior mean estimate of the root age.
- **`age(root)_ESS`**: ESS of the root age estimate.
- **`clock.rate_mean`**: Posterior mean estimate of the molecular clock rate.
- **`clock.rate_ESS`**: ESS of the molecular clock rate.
- **`age(Sample)_mean`**: Posterior mean estimate of the sample age.
- **`age(Sample)_Stdev`**: Posterior standard deviation of the estimated sample age.
- **`age(Sample)_HPD_95_Lower`** and **`age(Sample)_HPD_95_Upper`**: Lower and upper bounds of the 95% highest posterior density (HPD) interval for the sample age.
- **`age(Sample)_ESS`**: ESS of the estimated sample age.

The ESS is particularly important for assessing MCMC performance. Low ESS values indicate that the chain has not produced enough effectively independent samples to reliably characterize the posterior distribution.

As a general guideline, an ESS of at least 200 is commonly recommended for key parameters.

For a more detailed explanation of BEAST output and convergence diagnostics, see the [BEAST guide to analyzing output](https://beast.community/analysing_beast_output).

The original BEAST log files are available in `03_beast/` for further inspection, for example using Tracer.

### 6.2 Evaluate the Dating Results

Remember that we used only 1000 MCMC iterations in this tutorial to demonstrate the pipeline. This is a very short chain for a Bayesian phylogenetic dating analysis.

As expected, the ESS values in our test output are extremely low, mostly between 3 and 5. These values indicate that the analyses have not adequately explored the posterior distributions.

**These results should not be interpreted as reliable molecular age estimates.**

In a full analysis, at least 200 ESS for each sample is consider converged.

If certain samples show poor convergence or unexpected age estimates, you may want to rerun them using longer MCMC chains or reconsider their prior settings where scientifically justified.

Fortunately, MitoDate allows you to rerun selected samples without repeating the entire workflow.

Let's move on to **Mode 2 — Rerun Selected Samples**!




## 7. Mode 2 — Rerun Selected Samples

After reviewing the results from the initial tip-dating analysis, you may find that some samples require additional MCMC iterations or adjustments to their prior settings.

Instead of rerunning the entire dataset, MitoDate allows you to rerun selected samples individually.

### 7.1 Select Samples for Rerunning

For example, suppose we want to rerun samples `P048` and `MD228`.

Return to the configuration file and change the following parameters:

```groovy
run_mode      = 'rerun_samples'
rerun_samples = 'P048,MD228'
```

The `rerun_samples` parameter accepts a comma-separated list of sample IDs. Make sure that these identifiers match those in your input dataset.

You can also adjust the MCMC settings, such as increasing `chain_length`, to improve convergence. If necessary, modify the corresponding entries in the tip-dating prior file, provided that any changes are supported by the available chronological evidence.

All other parameters can remain unchanged unless additional modifications are required.

### 7.2 Launch the Rerun

Save the updated configuration file and launch the pipeline:

```bash
nextflow run /path/to/MitoDate/main.nf \
    -c /path/to/workdir/ChaconDuque_2025_10_withoutPartition_rerun.config \
    -profile mjolnir_test
```

Nextflow will display output similar to the following:

```text
N E X T F L O W   ~  version 25.10.4

Launching `/path/to/MitoDate/main.nf` [goofy_raman] DSL2 - revision: df02260152

        Workflow run name: ChaconDuque_2025_10-withoutPartition
        Selected run mode: rerun_samples

        Rerunning tip-dating for samples: P048,MD228

executor >  slurm (7)
[39/76fac3] GENERATEXML (Chacon-Duque-2025-10-VNTRremoved.TipDating.MD228.fasta)    [100%] 2 of 2 ✔
[13/a79acb] BEAST_TIP_DATING (Chacon-Duque-2025-10-VNTRremoved.TipDating.MD228.xml) [100%] 2 of 2 ✔
[1a/24964f] BEAST_LOG_PARSER (Chacon-Duque-2025-10-VNTRremoved.TipDating.MD228)     [100%] 2 of 2 ✔
[fc/e5875f] COLLECT_RESULTS (yyyy-mm-dd_xx-xx-xx)                                   [100%] 1 of 1 ✔

            Thank you for using MitoDate.

            Results are located in the folder: results.

Completed at: dd-mm-yyy xx:xx:xx
Duration    : 12m 26s
Succeeded   : 7
```

Notice that only two samples were processed, rather than all ten samples from the original analysis.

The rerun mode only performs molecular dating for the selected samples.

### 7.3 Review and Combine Age Estimates

Once the rerun is complete, examine the new age estimates and their convergence diagnostics in the results directory.

You can compare the new estimates with those from the initial run and decide which results to retain.

After reviewing all samples, combine the accepted age estimates into a single CSV file. This can be done by copying the corresponding rows from the age summary files of different runs.

When preparing the combined file, make sure that:

- Each sample appears only once.
- Only results with satisfactory convergence diagnostics are retained.
- The original column names and CSV format are preserved.
- The final file includes all undated samples required for the joint-tree analysis.

This combined age summary will serve as the input for the next mode.

**For this tutorial, we provide a precomputed age summary generated using MCMC chains of 100 million iterations.** You can find the file here:

```text
assets/test_data/Chacon-Duque-2025/Results_ageSummary.csv
```

We will use this file in the next step, rather than the preliminary estimates produced by our short test runs.

## 8. Mode 3 — Joint Tree Reconstruction

Now that we have obtained molecular age estimates for the undated samples, we can move on to the final step: reconstructing a time-calibrated phylogeny!

In `joint_tree` mode, MitoDate combines the molecular age estimates of previously undated samples with the known ages of radiocarbon-dated samples.

The estimated ages and their 95% highest posterior density (HPD) intervals are used to define age constraints for the previously undated sequences. These are then incorporated into a joint BEAST analysis to reconstruct a time-calibrated phylogenetic tree.

### 8.1 Configure Joint-Tree Mode

Return to the configuration file and update the following parameters:

```groovy
run_mode = 'joint_tree'

age_summary_file = '/path/to/MitoDate/assets/test_data/Chacon-Duque-2025/Results_ageSummary.csv'
```

Here:

- **`run_mode`** specifies that the pipeline should perform joint-tree reconstruction.
- **`age_summary_file`** points to the combined age summary containing the molecular dating results from the previous analyses.

For this example, we use the precomputed age summary provided with the test dataset.

Other analysis settings, including the substitution model, clock model, and MCMC settings, can be adjusted as needed.

### 8.2 Launch the Joint-Tree Analysis

Save the configuration file and run:

```bash
nextflow run /path/to/MitoDate/main.nf \
    -c /path/to/workdir/ChaconDuque_2025_10_withoutPartition_jointTree.config \
    -profile mjolnir_test
```

Once the workflow has completed successfully, you should see output similar to the following:

```text
N E X T F L O W   ~  version 25.10.4

Launching `/path/to/MitoDate/main.nf` [peaceful_payne] DSL2 - revision: df02260152

        Workflow run name: ChaconDuque_2025_10-withoutPartition
        Selected run mode: joint_tree

        Making date-calibrated tree using estimated tip dates from
        /path/to/MitoDate/assets/test_data/Chacon-Duque-2025/Results_ageSummary.csv

executor >  slurm (5)
[07/005f66] JOINT_TREE_META (joint_tree_meta)           [100%] 1 of 1 ✔
[a8/afa1a3] JOINT_XML (joint.renamed.fasta)             [100%] 1 of 1 ✔
[91/6e147e] BEAST_JOINT_TREE (joint.renamed.xml)        [100%] 1 of 1 ✔
[7b/f6304a] TREE_ANNOTATOR (joint.renamed)              [100%] 1 of 1 ✔
[fe/7a21d8] JOINT_TREE_VISUALIZE (joint_tree_visualize) [100%] 1 of 1 ✔

            Thank you for using MitoDate.

            Results are located in the folder: results.

Completed at: dd-mm-yyyy xx:xx:xx
Duration    : 3m 27s
Succeeded   : 5
```

The joint-tree workflow consists of five processes:

1. **`JOINT_TREE_META`**: Prepares the metadata and age information for joint-tree reconstruction.
2. **`JOINT_XML`**: Generates the BEAST XML configuration for the combined dataset.
3. **`BEAST_JOINT_TREE`**: Performs Bayesian phylogenetic inference using the combined sequence data and age constraints.
4. **`TREE_ANNOTATOR`**: Summarizes the posterior tree distribution to produce an annotated phylogenetic tree.
5. **`JOINT_TREE_VISUALIZE`**: Generates a visualization of the resulting time-calibrated phylogeny.

### 8.3 Explore the Final Tree

Once the analysis is complete, navigate to the results directory to explore the generated tree files and visualization.

MitoDate produces an annotated phylogenetic tree and a corresponding visualization, allowing you to examine the evolutionary relationships and divergence times of the analyzed samples.

You can also use the output tree files for further visualization or downstream phylogenetic analyses.

**Your time-calibrated phylogeny is now ready in 06_final_tree/05_tree_plot!**


## 9. From a Test Run to a Full Analysis

Congratulations! You have now completed all three modes of the MitoDate workflow:

1. **Tip dating:** Estimate the molecular ages of undated samples.
2. **Rerun selected samples:** Refine individual age estimates when necessary.
3. **Joint tree:** Reconstruct a time-calibrated phylogeny using the estimated sample ages.

Before applying MitoDate to your own dataset, there are a few important considerations.

### 9.1 MCMC Convergence

The MCMC chains used in this tutorial are intentionally short to demonstrate the workflow without requiring extensive computational resources.

In a full analysis, the MCMC chain length should be increased until adequate convergence and mixing are achieved. For example, for the precomputed tip-dating age estimates provided in this tutorial, we used MCMC chains of 100 million iterations.

As a general guideline, an effective sample size (ESS) greater than 200 is recommended for key posterior parameters. However, ESS alone is not sufficient to establish convergence. Trace plots, posterior distributions, and consistency across independent runs should also be examined.

Note that the required chain length may vary substantially depending on dataset size, sequence divergence, model complexity, and other factors. A fixed chain length does not guarantee convergence.

### 9.2 Computational Resources

Longer MCMC chains and more complex analyses can require substantially greater computational resources than the short demonstration runs.

Before launching a full analysis, make sure that your computing environment is configured with sufficient resources, particularly:

- **CPU and memory:** Allocate appropriate resources for BEAST and other processes.
- **Walltime:** Ensure that the requested job duration is sufficient for the expected MCMC runtime.
- **Storage:** Allow enough space for intermediate files, MCMC logs, and posterior tree samples.

If you are using an HPC cluster, review your Nextflow execution profile and adjust the requested resources according to the requirements of your dataset.

You are now ready to apply MitoDate to your own dataset.

Good luck with your analyses!



