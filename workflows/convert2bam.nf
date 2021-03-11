#!/usr/bin/env nextflow

/* 
 * Workflow to convert file/s in the CRAM format to the BAM format
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * More information on how to run the workflow can be found in:
 *
 * https://github.com/biodata-fun/nf_genomics_pipelines/wiki/cov_per_base.nf
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

// params defaults
params.help = false

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to convert files in the CRAM format to the BAM format'
    log.info '--------------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow -C convert2bam.config run convert2bam.nf --input \'/path/to/cram/*.cram\''
    log.info ''
    log.info 'Options:'
    log.info '    --help        Show this message and exit.'
    log.info '    --file        File containing the paths to the alignment files. The format of the file will be:'
    log.info '                          path/to/file1.bam\n'
    log.info '              path/to/file2.bam\n'
    log.info '              ...'
    log.info '                          Note. Each bam file needs an index.'
    log.info '    --genome      Absolute path for file with genome as required by BEDTools.'
    log.info '                  The genome file spec is defined by BEDTools and has the following format:'
    log.info '                   <seq_name><\t><chrom_size>'
    log.info '    --window      Number of bases for each of the BEDTools generated windows.'
    log.info '    --outprefix   Prefix for output file.'
    log.info ''
    exit 1
}

log.info 'Starting the conversion.....'

files = Channel.fromPath(params.input)

process convert2bam {
    /*
    Process to run samtools view to convert the .cram file to .bam format
        
	It will also create an index for the converted BAM file
        
	Returns
    -------
    Path to a BAM file
    */
	tag { cram_f }

    publishDir "converted_bam", mode: 'copy', overwrite: true

    input:
    file cram_f from files

    output:
    file "${cram_f}.bam" into out_bam

    """
    samtools view -b -o ${cram_f}.bam ${cram_f}
    """
}
