process SNPAD {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(input)

    output:
    path("priors*"),   emit: priors
    path("errors*"),   emit: errors

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    snpAD \
        --cpus 30 \
        --priors_out priors_${prefix}.txt \
        --errors_out errors_${prefix}.txt \
        $input
    """
}