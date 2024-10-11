process CALL {
    tag "$meta.id"
    label 'process_single'
    scratch true

    input:
    path(priors)
    path(errors)
    tuple val(meta), path(inputs)
    val(ch_i)
    val(ref_fasta_fai)

    output:
    tuple val(meta), path("*.vcf"),   emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    snpADCall \
    --name ${prefix} \
    --error $errors \
    --priors $priors \
    ${prefix}_chr$ch_i"_mapped.snpAD" > ${prefix}_chr$ch_i".vcf" \
    && bcftools reheader \
     --fai $ref_fasta_fai \
     ${prefix}_chr$ch_i".vcf" \
     --threads 8 \
     -o ${prefix}_chr$ch_i"_tmp" \
    && mv ${prefix}_chr$ch_i"_tmp" ${prefix}_chr$ch_i".vcf"
    """
}