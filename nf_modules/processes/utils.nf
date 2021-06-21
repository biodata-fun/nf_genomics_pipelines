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

process SAVE_VCF_FILE {
    /*
    This process will save a certain VCF
    Returns
    -------
    Path to saved VCF
    */
    
    publishDir "results/", mode: 'copy', overwrite: true

    executor 'local'

    input:
	path vcf

    output:
    path "out.norm.vcf.gz"

    """
    ls $vcf
    """
}