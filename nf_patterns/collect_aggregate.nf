/*
* Script to show how a process emits different files and these
* are collected and aggregated by a downstream process
*
* This script relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
*
* USAGE: nextflow run collect_aggregate.nf
*/


files = Channel.from 'fileA', 'fileB', 'fileC'

process generate_files {
    """
    Create the text files
    """
    input:
        val(f) from files
    
    output:
        file("${f}.txt") into aggr_ch
    
    script:
    """
    echo "hello" > ${f}.txt
    """
}

allf_ch = aggr_ch.collect()

process print_merged{
    """
    Process receiving the collected files
    and creating a single file with the 
    contents of each of the collected files
    """
    input:
    file(f) from allf_ch

    script:
    """
    cat $f > merged.txt
    """
}