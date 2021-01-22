process SAVE_FILE {
    /*
    This process will split the multiallelic variants by using BCFTools

    Returns
    -------
    Path to splitted VCF
    */
    
    publishDir "results/", mode: 'move', overwrite: true

    executor 'local'

    input:
	path(afile)
    val(prefix)

    output:
    path "${prefix}.merged.bam"

    """
    echo $afile $prefix
    mv ${afile} ${prefix}.merged.bam
    """
}