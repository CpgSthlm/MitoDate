# Data cleanup
DNAharvester produces many large files. Nextflow stores all
output files (including intermediate files) in a work directory.
The path to the work directory has to be specified in the config
file (default: `work` in the pipeline directory). Results files
are copied or symlinked from the work directory into a results directory
that is also specified in the config file. After a successful
run of DNAharvester, the work directory can be safely deleted.
Alternatively, the command `nextflow clean` can be used to
remove output from the work directory, the `.nextflow.log`
file, and clean the `.nextflow` cache directory (see
https://www.nextflow.io/docs/latest/reference/cli.html#clean
for details and options).

> The content of the work directory is used as a cache for the
output, so whenever the pipeline is re-run, results are copied/
symlinked from there unless anything has changed in the input
data, configuration, or code. If the work directory is deleted,
all its content is re-created the next time the pipeline
is executed again.
