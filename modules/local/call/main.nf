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
    def accesibility = task.ext.acc ?: "${params.accesibility}"
    """
    # Conditionally set poplistname based on accesibility
    if [ '$accesibility' == 'true' ]; then
    snpADCall \
    --name ${prefix} \
    --error $errors \
    --priors $priors \
    --set_zero_priors 0.0000000003125 \
    ${prefix}_chr$ch_i"_mapped_strict1kgp.snpAD" > ${prefix}_chr$ch_i".vcf" \
    && bcftools reheader \
     --fai $ref_fasta_fai \
     ${prefix}_chr$ch_i".vcf" \
     --threads 8 \
     -o ${prefix}_chr$ch_i"_tmp" \
    && mv ${prefix}_chr$ch_i"_tmp" ${prefix}_chr$ch_i".vcf"
    else
    snpADCall \
    --name ${prefix} \
    --error $errors \
    --priors $priors \
    --set_zero_priors 0.0000000003125 \
    ${prefix}_chr$ch_i"_mapped.snpAD" > ${prefix}_chr$ch_i".vcf" \
    && bcftools reheader \
     --fai $ref_fasta_fai \
     ${prefix}_chr$ch_i".vcf" \
     --threads 8 \
     -o ${prefix}_chr$ch_i"_tmp" \
    && mv ${prefix}_chr$ch_i"_tmp" ${prefix}_chr$ch_i".vcf"
    fi 
    """
}