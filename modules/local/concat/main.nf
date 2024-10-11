process CONCAT {
    tag "$meta.id"
    label 'process_high_memory'

    input:
    tuple val(meta), path(vcfs)
    path(fasta_ref)

    output:
    tuple val(meta), path("*.sorted.vcf.gz*"),   emit: sorted_vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    memory=\$(echo "$task.memory" | sed "s# GB#G#")
    bcftools concat \\
    $vcfs \\
    --threads ${task.cpus} \\
    | bcftools filter \\
    -i 'FMT/GQ>30' \\
    --set-GTs . \\
    --threads ${task.cpus} \\
    |  bcftools view \\
    --exclude-uncalled \\
    --threads ${task.cpus} \\
    | bcftools sort \\
    --max-mem \${memory} \\
        -Oz \\
        -o ${prefix}.sorted_tmp.vcf.gz \\
 && bcftools index \\
     --tbi ${prefix}.sorted_tmp.vcf.gz \\
     --threads ${task.cpus} \\
 && bcftools annotate \\
    --annotations $args \\
    --columns ID \\
    --set-id +'%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' \\
    --threads ${task.cpus} \\
    ${prefix}.sorted_tmp.vcf.gz \\
| bcftools norm \\
--check-ref s \\
        --fasta-ref $fasta_ref \\
        --threads ${task.cpus} \\
    -Oz -o ${prefix}.sorted.vcf.gz \\
 && bcftools index \\
    --tbi ${prefix}.sorted.vcf.gz \\
    --threads 8
    """
}
