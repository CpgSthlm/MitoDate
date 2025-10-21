process RENAMEFASTA {
    tag "${fasta}"
    label 'process_low'

    conda "conda-forge::biopython=1.81"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/biopython_openpyxl_pandas:fb650661820f6788' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path(fasta)
    path(metadata)

    output:
    path('*_renamed.fasta')             , emit: renamed_fasta
    path('versions.yml')                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in {{ name }}/bin/
    def output_prefix = fasta.getBaseName().replace('.fasta','').replace('.fa','')

    """
    renameFastaByMeta.py \\
        -i ${fasta} \\
        -m ${metadata} \\
        -o ${output_prefix}_renamed.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}