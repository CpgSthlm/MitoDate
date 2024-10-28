process BEAST {
    label 'process_high'

    input:
    each path(xml)

    output:
    path('*')                   , emit: beast_output
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