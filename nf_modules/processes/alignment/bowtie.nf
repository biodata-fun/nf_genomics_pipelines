/*
 *  Module containing processes for running Bowtie2
 *
 * This module relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

process RUN_BOWTIE {
    /*
    This process will align a pair of FASTQ using Bowtie2

    Returns
    -------
    Path to BAM with alignments
    */

    input:
        tuple val(sampleId), file(R1), file(R2)
        val(reference)

    output:
    path "${sampleId}.bam"

    script:
    """
    bowtie2 -x $reference -1 $R1 -2 $R2 | samtools view -bS - > ${sampleId}.bam
    """
}