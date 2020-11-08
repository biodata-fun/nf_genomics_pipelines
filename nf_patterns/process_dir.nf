/*
*  This script will glob all *.txt files in dir that will be
*  processed with foo_glob
*/

my_ch = Channel.fromPath( '/path/to/*.txt' )

process foo_glob {
    /*
    * Process that process each of the *.txt files
    * in a dir
    */
    input:
    file(afile) from my_ch

    script:
    """
    echo this is file ${afile}
    """
}