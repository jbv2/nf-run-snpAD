process BCFTOOLS_STATS {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.id}.stats"), emit: stats

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools stats -s - $vcf > ${prefix}.stats
    """
}