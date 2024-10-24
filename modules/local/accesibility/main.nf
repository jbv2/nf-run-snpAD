process ACCESIBILITY {
    tag "$meta.id"
    label 'process_single'
    scratch false

    input:
    tuple val(meta), path(input)
    val(ch_i)
    path(accesibility_bed)

    output:
    tuple val(meta), path("*_mapped_strict1kgp.snpAD"), val(ch_i),   emit: strict_inputs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    intersectbed.pl \
        ${prefix}_chr$ch_i"_mapped.snpAD" \
        ${accesibility_bed}/$ch_i".bed" > ${prefix}_chr$ch_i"_mapped_strict1kgp.snpAD"
    """
}