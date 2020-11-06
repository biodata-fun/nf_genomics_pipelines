/*
* This script shows how to build 2 different command lines depending on
* a parameter named mode. It illustrates how to implement conditional 
* execution in nextflow
*
* USAGE: nextflow run conditional_execution.nf --mode SE
*
*/

process cond_execution {
    script:
    if ( params.mode == 'SE')
    """
    echo fastqc -q read1
    """
    else if ( params.mode == 'PE')
    """
    echo fastqc -q read1 read2
    """
}

process cond_execution_w_output {
    /*
    This process will build 2 different command lines
    and will generate 2 different output files depending
    on params.mode
    */

    if (params.mode == 'SE'){
    output:
    file "out.se.txt"
    } else if ( params.mode == 'PE'){
    output:
    file "out.pe.txt"
    }
    
    script:
    if ( params.mode == 'SE')
    """
    touch out.se.txt
    """
    else if ( params.mode == 'PE')
    """
    touch out.pe.txt
    """
}