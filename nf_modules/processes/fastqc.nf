/*
 *  Module containing a process for running FASTQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
 *  with a FASTQ (https://en.wikipedia.org/wiki/FASTQ_format) file
 *
 * This module relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

process FASTQC {
    /*
    This process that runs FASTQC 

    Returns
    -------
    Path to splitted VCF
    */

    publishDir "results/", mode: 'move', overwrite: true
    executor 'local'

    input:
      tuple val(sampleId), file(R1), file(R2)

    output:
    file "${sampleId}*"
    
    script:
    """
    fastqc $R1 $R2
    """
}