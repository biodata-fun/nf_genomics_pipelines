/* 
 *  Module containing processes for executing the different
 *  bcftools commands on a VCF
 *
 * This module relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

process SPLIT_MULTIALLELIC {
    /*
    This process will split the multiallelic variants by using BCFTools

    Returns
    -------
    Path to splitted VCF
    */
    input:
    path vcf
    val(threads)

    output:
    path 'out.bcftoolsnorm.vcf.gz'

    """
    bcftools norm -m -any ${vcf} -o out.bcftoolsnorm.vcf.gz -Oz --threads ${threads}
    """
}

process GET_HEADER {
    /*
    Process is used to get the header from a VCF file

    Returns
    -------
    Path to a text file with a header
    */
    input:
    path vcf

    output:
    path 'header.txt'

    """
    bcftools view -h ${vcf} > header.txt
    """
}

process REHEADER {
    /*
    Process to replace the header of a VCF file

    Parameters
    ----------
    header.txt : file with header
    vcf : vcf that will be altered

    Output
    ------
    gzipped file
    */
    input:
    path(header)
    path(vcf)

    output:
    path 'reheaded.vcf.gz'

    """
    bcftools reheader -h ${header} -o 'reheaded.vcf.gz' ${vcf}
    """
}

process SELECT_VARIANTS {
    /*
    Process to select the desired variants type (snps/indels)
    */
    input:
    path(vcf)
    val(vt)
    val(threads)

    output:
    path "out.${vt}.vcf.gz"

    """
    bcftools view $vcf -v ${vt} -o out.${vt}.vcf.gz -Oz --threads ${threads}
    """
}

process RUN_BCFTOOLS_SORT {
    /*
    Process to run bcftools sort
    */
    input:
    path(vcf)

    output:
    path "out.sort.vcf.gz"

    """
	mkdir tmpdir
	bcftools sort -T tmpdir/ ${vcf} -o out.sort.vcf.gz -Oz
    """
}

process EXC_NON_VTS {
    /*
    This process will select the variants on a VCF

    Returns
    -------
    Path to VCF containing only variants
    */
    input:
    path(vcf)
    val(threads)

    output:
    file "out.onlyvariants.vcf.gz"

    """
    bcftools view -c1 ${vcf} -o out.onlyvariants.vcf.gz --threads ${threads} -Oz
    """
}

process INTERSECTION_CALLSET {
    /*
    Process to find the intersection between a call set and the Gold
    standard call set

    Parameters
    ----------
    vcf : path to VCF file used for training
    vt : type of variant ('snps','indels')
    true_cs : path to gold standard call set
    true_cs_ix : path to index for 'true_cs'

    Output
    ------
    3 VCFs containing the False Positives in FP.vcf.gz, the False Negatives in FN.vcf.gz
    and the True Positives in TP.vcf.gz
    */
    input:
    path(vcf)
    val(vt)
    path(true_cs)
    path(true_cs_ix)

    output:
    path 'FP.vcf.gz', emit: fp_vcf
    path 'FN.vcf.gz', emit: fn_vcf 
    path 'TP.vcf.gz', emit: tp_vcf

    """
    tabix ${vcf}
    bcftools isec -c ${vt}  -p 'dir/' ${vcf} ${true_cs}
    bgzip -c dir/0000.vcf > FP.vcf.gz
    bgzip -c dir/0001.vcf > FN.vcf.gz
    bgzip -c dir/0002.vcf > TP.vcf.gz
    """
}

process BCFT_QUERY {
    /*
    Process to run 'bcftools query' on a VCF file

    Output
    ------
    .tsv file compressed with gzip containing the annotations
    that have been 
    */
    input:
    path(vcf)
    val(annotations)

    output:
    path("${vcf.baseName}.tsv.gz")

    """
    bcftools query -H -f '${annotations}' ${vcf} | bgzip -c > ${vcf.baseName}.tsv.gz
    """
}

process SELECT_REGION {
    /*
    Process to fetch a certain region from the VCF using 'bcftools view'

    Parameters
    ----------
    vcf : path to VCF file
    vcf_ix : path to tabix index for this VCF file
    region : region to fetch. For example: chr20

    Output
    ------
    gzipped VCF file with a particular region
    */
    input:
    path(vcf)
    path(vcf_ix)
    val(region)

    output:
    path("sub.${region}.vcf.gz")

    """
    bcftools view -r ${region} ${vcf} -o sub.${region}.vcf.gz -Oz
    """
}

process DROP_GTPS {
    /*
    Process to drop the genotype information from a VCF

    Parameters
    ----------
    vcf : path to the VCF file
    vt : select/exclude comma-separated list of variant types: snps,indels,mnps,ref,bnd,other

    Output
    ------
    gzipped VCF file without genotypes
    */
    input:
    path(vcf)
    val(vt)
    val(threads)

    output:
    path("out.${vt}.noGTPS.vcf.gz")

    """
    bcftools view -G -v ${vt} ${vcf} -o out.${vt}.noGTPS.vcf.gz -Oz --threads ${threads}
    """
}