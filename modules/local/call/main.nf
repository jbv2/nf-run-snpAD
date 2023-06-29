process CALL {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(priors)
    tuple val(meta), path(errors)
    tuple val(meta), path(inputs)

    output:
    tuple val(meta), path("*.vcf"),   emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for i in `seq 1 22` X Y ;
    do
    snpADCall \
    --name ${prefix} \
    --error $errors \
    --priors $priors \
    ${prefix}_chr\${i}_mapped.snpAD > ${prefix}_chr\${i}.vcf
    """
}