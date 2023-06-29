process SNPAD {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("priors*"),   emit: priors
    tuple val(meta), path("errors*"),   emit: errors

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
        ${prefix}_chr21_mapped.snpAD
    """
}