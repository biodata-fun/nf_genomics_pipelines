/* 
* This script shows how to create a channel using a file containing file paths
* USAGE: nextflow run ch_from_paths_inafile.nf --list file_paths.txt
*/

listSeq=params.list

Channel.fromPath(listSeq)
        .splitCsv()
        .map { row ->
            def path = file("${row[0]}")
            return path
        }
        .set { file_list}


process foo_per_line {
    /*
    * Process that process each of lines in 
    */
    
    input:
    file x from file_list

    script:
    """
    cat $x
    """
}