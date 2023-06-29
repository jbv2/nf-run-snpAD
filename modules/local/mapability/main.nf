process MAPABILITY {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(input)
    path(map_bed)

    output:
    tuple val(meta), path("*.snpAD"),   emit: inputs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for i in `seq 1 22` X Y ;
    do
    intersectbed.pl ${prefix}_chr\${i}".snpAD" ${map_bed}/\${i}.bed > ${prefix}_chr\${i}_mapped.snpAD ;
    done
    """
}