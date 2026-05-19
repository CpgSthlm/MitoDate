# MitoDate

MitoDate is a Nextflow pipeline designed to perform molecular clock dating (also known as tip dating) of ancient mitochondrial genomes (mitogenomes) using BEAST (Bayesian Evolutionary Analysis by Sampling Trees) following the methodology published in Chacón-Duque et al., MBE (2025). The pipeline automates the following processes:

- Processing and preparing multiple sequence alignments (in FASTA format)
- Generating BEAST XML configuration files with customizable parameters
- Running Bayesian phylogenetic inference using BEAST
- Parsing and collecting results from BEAST analyses
- Supporting both single and multiple sample dating approaches


### Quick Navigation

1. **[Pipeline Requirements and Configuration](docs/1_requirements_and_configuration.md)** - System requirements, dependencies, and configuration options
2. **[Pipeline Steps](docs/2_pipeline_steps.md)** - Detailed explanation of each pipeline stage
3. **[How to Run MitoDate](docs/3_how_to_run.md)** - Step-by-step guide to execute the pipeline
4. **[Data Cleanup](docs/4_data_cleanup.md)** - Managing output files and cleanup procedures
5. **[Changelogs](docs/5_changelogs.md)** - Version history and recent updates


### Support

For questions, issues, or feature requests, please refer to
the relevant documentation sections or contact:
Wenxi Li (vanssy0819@gmail.com)
Bilal Sharif (bilal.bioinfo@gmail.com)

## Citation

If you use MitoDate in your research, please cite:




