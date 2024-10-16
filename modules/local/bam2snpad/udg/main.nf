process BAM2SNPAD_UDG {
    tag "$meta.id"
    label 'process_single'
    scratch true

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
        --bam_index $bai \
        --region $ch_i \
        --map_qual 30 \
        --fasta $ref_fasta \
        --offset 31 \
        $bam >${prefix}_chr$ch_i".snpAD" 
    """
}