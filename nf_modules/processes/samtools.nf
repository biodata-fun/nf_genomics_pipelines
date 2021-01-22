/*
 *  Module containing processes for running different Samtools commands
 *
 * This module relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

process RUN_SAMTOOLS_SORT {
    /*
    This process will sort an alignment file

    Returns
    -------
    Path to sorted BAM file
    */

    input:
    path bam_file

    output:
        path "${bam_file.baseName}.sorted.bam"

    script:
    """
    samtools sort ${bam_file} -O BAM -o ${bam_file.baseName}.sorted.bam
    """
}

process RUN_SAMTOOLS_MERGE {
    /*
    This process will run samtools merge
    to merge multiple sorted BAM files

    Returns
    -------
    Path to merged BAM file
    */

    input:
        path bams

    output:
    path "out.merged.bam"

    script:
    """
    samtools merge out.merged.bam ${bams} -O BAM
    """
}