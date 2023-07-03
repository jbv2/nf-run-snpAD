process MAPABILITY {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(input)
    val(ch_i)
    path(map_bed)

    output:
    tuple val(meta), path("*.snpAD"), val(ch_i),   emit: inputs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    intersectbed.pl \
        ${prefix}_chr$ch_i".snpAD" \
        ${map_bed}/$ch_i".bed" > ${prefix}_chr$ch_i"_mapped.snpAD"
    """
}