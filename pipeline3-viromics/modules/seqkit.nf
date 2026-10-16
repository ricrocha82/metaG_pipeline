process SPLIT_VOTUS {
    // iPHoP runtime grows with the number of input viruses; chunks run as
    // independent jobs and their predictions are concatenated afterwards.
    label 'process_low'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/seqkit/2.7.0/seqkit_2.7.0.sif"

    input:
    path votus

    output:
    path "chunks/*.fna", emit: chunks

    script:
    """
    # split2 keeps the input extension: chunks/vOTUs.part_001.fna, ...
    seqkit split2 -s ${params.iphop_chunk_size} -O chunks -j ${task.cpus} ${votus}
    """
}
