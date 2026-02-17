process BEAST {
    tag "${xml}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/beast%3A1.10.4--hdfd78af_2' :
        'community.wave.seqera.io/library/beast:1.10.4--2b8ba8c6c4be979a' }"

    input:
    path(xml)

    output:
    path('*.log')               , emit: beast_out_log
    path ('*.trees')            , emit: beast_out_trees
    path('versions.yml')        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Run BEAST
    beast -beagle -overwrite ${xml}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beast: \$(echo \$(beast -version))
    END_VERSIONS
    """
}