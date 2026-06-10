process GENERATEXML {
    tag "${fasta}"
    label 'process_generate_xml'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/biopython_openpyxl_pandas:fb650661820f6788' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path(fasta)
    path(metadata)
    path(priors)

    output:
    path('*.xml')               , emit: xml
    path('versions.yml')        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in {{ name }}/bin/

    def split_partition     = task.ext.split_partition ?: params.split_partition
    def substitution_model  = task.ext.substitution_model ?: params.substitution_model
    def root_mean           = task.ext.root_mean ?: params.root_mean
    def root_stdev          = task.ext.root_stdev ?: params.root_stdev
    def chain_length        = task.ext.chain_length ?: params.chain_length
    def log_every           = task.ext.log_every ?: params.log_every
    def gff                 = task.ext.gff ?: params.gff ?: ''
    def population_model    = task.ext.population_model ?: params.population_model
    def clock_model         = task.ext.clock_model ?: params.clock_model
    def root_offset         = task.ext.root_offset ?: params.root_offset

    if (!split_partition || split_partition.toString().trim().toLowerCase() == 'false') {
        """
        main_xml_generation.py \\
            -f ${fasta} \\
            -m ${metadata} \\
            -p ${priors} \\
            --subs_model "${substitution_model}" \\
            --root_mean ${root_mean} \\
            --root_stdev ${root_stdev} \\
            --chain_length ${chain_length} \\
            --log_every ${log_every} \\
            --population_model ${population_model} \\
            --clock_model ${clock_model} \\
            --root_offset ${root_offset}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """
    } else {
        """
        main_xml_generation.py \\
            -f ${fasta} \\
            -m ${metadata} \\
            -p ${priors} \\
            --subs_model "${substitution_model}" \\
            --root_mean ${root_mean} \\
            --root_stdev ${root_stdev} \\
            --chain_length ${chain_length} \\
            --log_every ${log_every} \\
            --split_partition \\
            --annotation "${gff}" \\
            --population_model ${population_model} \\
            --clock_model ${clock_model} \\
            --root_offset ${root_offset}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """
    }
}
