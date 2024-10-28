process GENERATEXML {

    label 'process_low'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    each path(fasta)
    path(priors)
    path(gff)
    path(partition)
    val(chainlength)
    val(log_step)
    val(partition_list)
    val(nd_list)
    val(taxon_set)

    output:
    path('*.xml')               , emit: xml
    path('versions.yml')        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in {{ name }}/bin/
    """
    main.py \\
        ${fasta} \\
        ${priors} \\
        ${gff} \\
        ${partition} \\
        ${chainlength} \\
        ${log_step} \\
        ${partition_list} \\
        ${nd_list} \\
        ${taxon_set}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
