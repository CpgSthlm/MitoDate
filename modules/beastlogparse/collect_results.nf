process COLLECT_RESULTS {
    tag "${params.workflow_timestamp}"
    label 'process_collect_results'

    conda "conda-forge::biopython=1.81"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/biopython_openpyxl_pandas:fb650661820f6788' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path(csv)

    output:
    path('*ageSummary.csv')             , emit: age_summary
    path('versions.yml')                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in {{ name }}/bin/

    """
    CollectResults.py \\
        ${csv} \\
        -o Results_ageSummary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}