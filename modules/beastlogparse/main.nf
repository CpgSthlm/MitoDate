process BEAST_LOG_PARSER {
    tag "${beast_logs.baseName}"
    label 'process_medium'

    conda "conda-forge::python=3.8.3, "
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/numpy_pandas_pathlib:4292f4f6e675421b' :
        'community.wave.seqera.io/library/numpy_pandas_pathlib:09d913c976849beb' }"

    input:
    path(beast_logs)

    output:
    path('*.csv')               , emit: parsed_results
    path('*.log')               , emit: parsed_log
    path('versions.yml')        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in {{ name }}/bin/

    def burnin      = task.ext.burnin ?: params.burnin
    def output_name = beast_logs.baseName.replaceAll('.log', '')

    """
    ParsingLogs.py \\
        ${beast_logs} \\
        -b ${burnin} \\
        -o ${output_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
