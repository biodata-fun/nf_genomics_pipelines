/* 
 * Workflow to run the BWA short read aligner. See http://bio-bwa.sourceforge.net/bwa.shtml for more 
 * details
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

// params defaults
params.help = false
params.threads = 1
params.queue = 'production-rh74'

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to run BWA on FASTQ file/s'
    log.info '-------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run_bwa.nf --dir DIR --ref REF_FILE --format BAM'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--dir DIR    Path to the dir containing the FASTQ file pairs to be aligned.'
    log.info '  --ref REF_FILE  Path to the reference FASTA file.'
    log.info '  --format FMT Output format sam, bam, cram.'
    log.info ''
    exit 1
}

log.info 'Starting the analysis.....'

Channel.fromFilePairs(params.dir+"/*_{1,2}.fastq.gz").set{fastqChannel}

process run_bwa {
    /*
    Process to run the BWA read mapper

    Returns
    -------
    A BAM file with the aligned reads
    */

    memory '5 GB'
    executor 'lsf'
    queue "${params.queue}"

    publishDir "alignments/", mode: 'copy', overwrite: true

    input:
        set val(runId), file(reads) from fastqChannel
    
    output:
        file "${runId}.${params.format}"


    script:
    if (params.format == 'sam')
    """
    bwa mem ${params.ref} ${reads} -o ${runId}.${params.format}
    """
    else if (params.format == 'bam')
    """
    bwa mem ${params.ref} ${reads} |samtools view -b -o ${runId}.${params.format}
    """
    else if (params.format == 'cram')
    """
    bwa mem ${params.ref} ${reads} |samtools view -C -T ${params.ref} > ${runId}.${params.format}
    """
}