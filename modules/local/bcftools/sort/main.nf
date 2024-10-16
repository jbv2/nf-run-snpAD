process BCFTOOLS_SORT {
    tag "$meta.id"
    label 'process_high_memory'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.sorted.vcf.gz"), emit: sorted_vcf

    script:
    def memory = task.memory.toMega()  // Use memory in megabytes
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools sort \\
        --max-mem ${memory}M \\
        -Oz -o ${prefix}.sorted.vcf.gz \\
        ${vcf}
    """
}
