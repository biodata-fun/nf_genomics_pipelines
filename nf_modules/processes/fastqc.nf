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

    executor 'local'

    input:
      tuple val(sampleId), file(reads)
    
    //output:
    //path 'out.bcftoolsnorm.vcf.gz'
    script:
    """
    echo fastqc --sample $sampleId --reads $reads 
    """
}