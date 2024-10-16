process BCFTOOLS_FILTER {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf)
    path(fasta_ref)

    output:
    tuple val(meta), path("${meta.id}.filtered.vcf.gz"), emit: filtered_vcf

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools filter \\
        -i 'FMT/DP > 3 && FMT/GQ > 29' \\
        --set-GTs . \\
        --threads ${task.cpus} \\
        ${vcf} \\
    | bcftools view \\
    --exclude-uncalled \\
    --threads ${task.cpus} \\
    | bcftools norm \\
    --check-ref s \\
        --fasta-ref $fasta_ref \\
        --threads ${task.cpus} \\
        -Oz -o ${prefix}.filtered.vcf.gz \\
    """
}
