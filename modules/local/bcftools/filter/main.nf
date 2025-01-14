process BCFTOOLS_FILTER {
    tag "$meta.id"
    label 'process_low'
    scratch true

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
        -i 'FMT/DP > 1 && FMT/GQ > 29' \\
        --set-GTs . \\
        --threads ${task.cpus} \\
        ${vcf} \\
    bcftools filter \\
        -i 'GT="0/1" && ((ALT="A" && FORMAT/A/FORMAT/DP >= 0.2) \\
        || (ALT="C" && FORMAT/C/FORMAT/DP >= 0.2) \\
        || (ALT="G" && FORMAT/G/FORMAT/DP >= 0.2) \\
        || (ALT="T" && FORMAT/T/FORMAT/DP >= 0.2)) \\
        || GT="0/0" \\
        || GT="1/1"'
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
