process VOTU_CLUSTERING {
    // CheckV's recommended ANI clustering (blastn -> anicalc -> aniclust;
    // scripts copied from the CheckV repo into bin/). aniclust uses the
    // longest sequence of each cluster as the representative (vOTU).
    // The CheckV image provides blastn, Biopython and numpy.
    label 'process_high'
    publishDir "${params.outdir}/06_votus", mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/CheckV/0.8.1/CheckV-0.8.1.sif"

    input:
    path viral_contigs

    output:
    path "vOTUs.fna",          emit: votus
    path "vOTU_clusters.tsv",  emit: clusters
    path "viral_contigs_ani.tsv", emit: ani

    script:
    """
    makeblastdb -in ${viral_contigs} -dbtype nucl -out blastdb
    blastn \\
        -query ${viral_contigs} \\
        -db blastdb \\
        -outfmt '6 std qlen slen' \\
        -max_target_seqs 10000 \\
        -num_threads ${task.cpus} \\
        -out viral_contigs_blast.tsv

    anicalc.py -i viral_contigs_blast.tsv -o viral_contigs_ani.tsv
    aniclust.py \\
        --fna ${viral_contigs} \\
        --ani viral_contigs_ani.tsv \\
        --out vOTU_clusters.tsv \\
        --min_ani ${params.votu_min_ani} \\
        --min_tcov ${params.votu_min_tcov} \\
        --min_qcov 0

    # Representatives = first column of the cluster table
    cut -f1 vOTU_clusters.tsv > representatives.txt
    awk 'NR==FNR { keep[\$1]; next } /^>/ { p = (substr(\$1, 2) in keep) } p' \\
        representatives.txt ${viral_contigs} > vOTUs.fna
    """
}
