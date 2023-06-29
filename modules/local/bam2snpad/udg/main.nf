process BAM2SNPAD_UDG {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(bai)
    path(ref_fasta)

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
    Bam2snpAD $bam --bam_index $bai --region \${i} --map_qual 25 --fasta $ref_fasta > ${prefix}_chr\${i}".snpAD";
    done
    """
}