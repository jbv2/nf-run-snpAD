process BCFTOOLS_ANNOTATE {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf)
    path(tbi)

    output:
    tuple val(meta), path("${meta.id}.annotated.vcf.gz"), emit: annotated_vcf

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools annotate \\
        --annotations $args \\
        --columns ID \\
        --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' \\
        --threads ${task.cpus} \\
        -Oz -o ${prefix}.annotated.vcf.gz \\
        ${vcf}
    """
}
