process JOINT_TREE_VISUALIZE {
    tag "joint_tree_visualize"
    label 'process_joint_tree_visualize'

    conda "conda-forge::biopython=1.81 conda-forge::matplotlib=3.7"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/matplotlib_pip_biopython:df1c4b0b702add4f' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path(tree)
    path(joint_meta)

    output:
    path('joint_tree.png')      , emit: tree_png
    path('joint_tree.pdf')      , emit: tree_pdf
    path('versions.yml')        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in bin/
    def figsize_w = task.ext.figsize_w ?: params.figsize_w ?: 12
    def figsize_h = task.ext.figsize_h ?: params.figsize_h ?: 12

    """
    joint_tree_visualize.py \\
        -i ${tree} \\
        -m ${joint_meta} \\
        -o joint_tree \\
        --figsize ${figsize_w} ${figsize_h}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
