#!/usr/bin/env nextflow

/* 
 * Script for splitting a VCF containing multiple chromosomes into 1 VCF per chromosome.
 * It will only consider the biallelic sites
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

// params defaults
params.help = false
params.prefix = "out"
params.outdir = "outdir"
params.cpus = 1

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to split a VCF into chromosomes'
    log.info '----------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow split_into_chros.nf --vcf callset.vcf.gz --chros chr1,chr2 --prefix callset --outdir dirname'
    log.info ''
    log.info 'Options:'
    log.info '  --help  Show this message and exit.'
    log.info '  --vcf VCF	VCF that will be splited.'
    log.info '	--prefix PREFIX		  String used for output files. Default: out'
    log.info '  --outdir DIRNAME          Namd of output dir.'
    log.info '  --chros LIST_OF_CHROS	  List of single chromosome-VCFs to generate.'
    log.info '  --cpus INT	  Number of CPUs to use. Default 1.'
    log.info ''
    exit 1
}

if( !params.chros) {
  exit 1, "Undefined --chros parameter. I need a comma-separated string with chromosomes: chr1,chr2, ..."
}
if( !params.vcf) {
  exit 1, "Undefined --vcf parameter. I need the path to a VCF file"
}

vcfFile = file(params.vcf)
if( !vcfFile.exists() ) {
  exit 1, "The specified VCF file does not exist: ${params.vcf}"
}

chr_str = params.chros.toString()
chr_list = Channel.from( chr_str.split(','))

process splitVCF {
	/*
	Function to split a multi-chromosome VCF into single chromosome VCF

	Returns
	-------
	Returns 2 files per chromosome:
		1) A VCF format file for each splitted chromosome
		2) A tabix index for that VCF
	*/
	
  maxForks = "${params.cpus}"
 	publishDir params.outdir, mode:'move'

	input:
	val chr from chr_list
	path vcf from params.vcf

	output:
	path("${params.prefix}.${chr}.vcf.gz*")

	"""
	bcftools index -t ${vcf}
	bcftools view -r ${chr} ${vcf} -o ${params.prefix}.${chr}.vcf.gz -O z
	bcftools index -t ${params.prefix}.${chr}.vcf.gz
	"""
}
