process JOINT_TREE_META {
    tag "joint_tree_meta"
    label 'process_joint_tree_meta'

    conda "conda-forge::biopython=1.81"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/biopython_openpyxl_pandas:fb650661820f6788' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path(age_summary)
    path(metadata)
    path(fasta)

    output:
    path('joint_meta.tsv')          , emit: joint_meta
    path('joint.renamed.fasta')     , emit: joint_fasta
    path('versions.yml')            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in bin/
    """
    JointTreeMeta.py \\
        -i ${age_summary} \\
        -m ${metadata} \\
        -f ${fasta} \\
        -o joint_meta.tsv \\
        -F joint.renamed.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
