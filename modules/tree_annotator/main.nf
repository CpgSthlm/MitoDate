process TREE_ANNOTATOR {
    tag "${trees.baseName}"
    label 'process_tree_annotator'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/beast%3A1.10.4--hdfd78af_2' :
        'community.wave.seqera.io/library/beast:1.10.4--2b8ba8c6c4be979a' }"

    input:
    path(trees)

    output:
    path('joint.tree.out')      , emit: annotated_tree
    path('versions.yml')        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def burnin       = task.ext.burnin       ?: params.burnin
    def chain_length = task.ext.chain_length ?: params.chain_length
    // treeannotator -burnin expects number of states to discard
    def burnin_states = (burnin.toFloat() / 100.0 * chain_length.toLong()).toLong()

    """
    treeannotator \\
        -burnin ${burnin_states} \\
        ${trees} \\
        joint.tree.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        treeannotator: \$(treeannotator 2>&1 | grep -o 'v[0-9][^ ]*' | head -1 || echo 'BEAST 1.10.4')
    END_VERSIONS
    """
}
