process CALL {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(priors)
    tuple val(meta), path(errors)
    tuple val(meta), path(inputs)
    val(ref_fasta_fai)

    output:
    tuple val(meta), path("*.vcf"),   emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for i in `seq 1 3`;
    do
    snpADCall \
    --name ${prefix} \
    --error $errors \
    --priors $priors \
    ${prefix}_chr\${i}_mapped.snpAD > ${prefix}_chr\${i}.vcf \
    && bcftools reheader \
     --fai $ref_fasta_fai \
     ${prefix}_chr\${i}.vcf \
     --threads 8 \
     -o ${prefix}_chr\${i}_tmp \
    && mv ${prefix}_chr\${i}_tmp ${prefix}_chr\${i}.vcf ;
    done
    """
}