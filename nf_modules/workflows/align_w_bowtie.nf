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
//include { SAVE_FILE } from '../processes/utils.nf'


// params defaults
params.help = false
params.threads = 1
params.outdir = './results'

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to run a BOWTIE aligment on a FASTQ file'
    log.info '-------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow align_w_bowtie.nf --dir DIR'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--DIR DIR    DIR containing the FASTQ files.'
    log.info '	--ref REF   File with reference genome.'
    log.info ''
    exit 1
}

log.info 'Starting alignment.....'

workflow  {
    samples_ch = Channel.fromFilePairs(params.dir+"*_R{1,2}*.fastq.gz")
    main:
        FASTQC( samples_ch)
        RUN_BOWTIE( samples_ch, params.ref)
//     SELECT_VARIANTS(ALLELIC_PRIMITIVES.out, params.vt, params.threads)
   //     RUN_BCFTOOLS_SORT(SELECT_VARIANTS.out)
   //     RUN_VT_UNIQ(RUN_BCFTOOLS_SORT.out)
   //    SAVE_FILE(RUN_VT_UNIQ.out)
}

