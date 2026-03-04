process SPLITFASTA {
    tag "${fasta}"
    label 'process_split_fasta'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/biopython_openpyxl_pandas:fb650661820f6788' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path(fasta)
    path(metadata)

    output:
    path('*TipDating*.fasta')             , emit: fastas
    path('versions.yml')                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in {{ name }}/bin/
    """
    singleDatingMSA.py \\
        -i ${fasta} \\
        -m ${metadata} \\
        -o ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
