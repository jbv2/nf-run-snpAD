process BCFTOOLS_CONCAT {
    tag "$meta.id"
    label 'process_medium'
    scratch true

    input:
    tuple val(meta), path(vcfs)

    output:
    tuple val(meta), path("*.concat.vcf.gz"), emit: concat_vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Concatenation of VCF files
    """
    bcftools concat \\
        ${vcfs} \\
        --threads ${task.cpus} \\
        -Oz -o ${prefix}.concat.vcf.gz
    """
}
