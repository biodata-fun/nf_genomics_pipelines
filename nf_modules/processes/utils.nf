process SAVE_FILE {
    /*
    This process will split the multiallelic variants by using BCFTools

    Returns
    -------
    Path to splitted VCF
    */
    
    publishDir "results/", mode: 'copy', overwrite: true

    executor 'local'

    input:
	path(vcf)
    val(vt)

    output:
    path "out.norm.${vt}.vcf.gz"

    """
    mv $vcf out.norm.${vt}.vcf.gz
    """
}