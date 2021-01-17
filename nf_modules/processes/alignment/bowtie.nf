/*
 *  Module containing processes for executing the different
 *  bcftools commands on a VCF
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

    executor 'local'

    input:
      tuple val(sampleId), file(reads)
      val(reference)

//    output:
//    path 'out.bcftoolsnorm.vcf.gz'

    script:
    """
    echo bowtie2 -x $reference -1 $reads[0] 
    """
}