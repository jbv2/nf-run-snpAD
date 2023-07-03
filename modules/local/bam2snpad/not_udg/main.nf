process BAM2SNPAD_NOT_UDG {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(bam)
    path(bai)
    path(ref_fasta)
    val(ch_i)

    output:
    tuple val(meta), path("*.snpAD"), val(ch_i),  emit: inputs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    Bam2snpAD \
        --region $ch_i \
        --bam_index $bai  \
        --map_qual 25 \
        --fasta $ref_fasta \
        $bam >${prefix}_chr$ch_i".snpAD" 
    """
}