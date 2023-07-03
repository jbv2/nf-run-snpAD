process CONCAT {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf)
    val(rsID)

    output:
    tuple val(meta), path("*.sorted.vcf.gz*"),   emit: sorted_vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools concat \\
    \$(ls ${prefix}*.vcf | sort -V ) \\
     --threads 8 \\
     | bcftools filter \\
    -i 'FMT/GQ>20' \\
    --set-GTs . \\
    --threads 8 \\
     | bcftools sort \\
        -Oz \\
        -o ${prefix}.sorted_tmp.vcf.gz \\
 && bcftools index \\
     --tbi ${prefix}.sorted_tmp.vcf.gz \\
     --threads 8 \\
 && bcftools annotate \\
    --annotations $rsID \\
    --columns ID \\
    --set-id +'%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' \\
    --threads 8 \\
    ${prefix}.sorted_tmp.vcf.gz \\
    -Oz -o ${prefix}.sorted.vcf.gz \\
 && bcftools index \\
    --tbi ${prefix}.sorted.vcf.gz \\
    --threads 8
    """
}