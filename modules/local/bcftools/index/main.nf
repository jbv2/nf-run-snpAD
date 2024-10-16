process BCFTOOLS_INDEX {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}*.tbi"), emit: index_file

    script:
    """
    bcftools index \\
        --tbi ${vcf} \\
        --threads ${task.cpus}
    """
}
