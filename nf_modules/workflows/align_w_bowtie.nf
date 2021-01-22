/* 
 * Workflow to align a FASTQ with Bowtie
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

nextflow.enable.dsl=2

include { FASTQC } from '../processes/fastqc.nf' 
include { RUN_BOWTIE } from '../processes/alignment/bowtie.nf'
include { RUN_SAMTOOLS_SORT } from '../processes/samtools.nf'
include { RUN_SAMTOOLS_MERGE } from '../processes/samtools.nf'
include { SAVE_FILE } from '../processes/utils.nf'


// params defaults
params.help = false
params.threads = 1
params.outdir = './results'

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to run a BOWTIE aligment on multiple FASTQ files'
    log.info '---------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow align_w_bowtie.nf --dir DIR'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--DIR DIR    DIR containing the FASTQ files.'
    log.info '	--ref REF   File with reference genome.'
    log.info '	--prefix PREFIX   Prefix for merged BAM file.'
    log.info ''
    exit 1
}

log.info 'Starting alignment.....'

workflow  {

    Channel
    .fromFilePairs(params.dir+"*_R{1,2}_001.fastq.gz", 
    flat:true, 
    checkIfExists:true)
    .set { f_ch}

    main:
        FASTQC( f_ch)
        RUN_BOWTIE( f_ch, params.ref)
        RUN_SAMTOOLS_SORT(RUN_BOWTIE.out)
        all_sorted_bams = RUN_SAMTOOLS_SORT.out.collect()
        RUN_SAMTOOLS_MERGE(all_sorted_bams)
        SAVE_FILE(RUN_SAMTOOLS_MERGE.out, params.prefix)
}
