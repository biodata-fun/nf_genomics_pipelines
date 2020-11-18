/* 
 * Workflow to calculate the MD5 on a set of files for which the
 * paths are passed in a file. 
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 */

 // params defaults
params.help = false
params.threads = 1
params.queue = 'production-rh74'

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to run md5sum on a set of file/s'
    log.info '-----------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run_md5.nf --list FILE'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--list FILE    File with the paths for the files being analysed.'
    log.info ''
    exit 1
}

log.info 'Starting the analysis.....'

listSeq=params.list

Channel.fromPath(listSeq)
        .splitCsv()
        .map { row ->
            def path = file("${row[0]}")
            return path
        }
         .set { file_list}

process md5_for_path {
    /*
    * Process to run 'md5sum' on each of the files

    * Returns
    * -------
    * 
    */
    memory '512 MB'
    executor 'lsf'
    queue "${params.queue}"
    maxForks 20

    publishDir "result_md5", mode: 'move', overwrite: true

    input:
    file x from file_list

    output:
    file "${x}.md5"

    script:
    """
    md5sum ${x} > ${x}.md5
    """
}