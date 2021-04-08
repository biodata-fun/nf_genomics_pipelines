/* 
 * Workflow to calculate the coverage per base from multiple alignment files
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * More information on how to run the workflow can be found in:
 *
 * https://github.com/biodata-fun/nf_genomics_pipelines/wiki/cov_per_base.nf
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

// params defaults
params.help = false
params.queue = 'production-rh74'
params.region = false

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to calculate the coverage per base for multiple alignment files'
    log.info '------------------------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow -C cov_per_base.config run cov_per_base.nf --file <aln_paths.txt> --genome <genome.txt> --window <int>'
    log.info ''
    log.info 'Options:'
    log.info '    --help	Show this message and exit.'
    log.info '    --file	File containing the paths to the alignment files. The format of the file will be:'
    log.info '				path/to/file1.bam\n'
    log.info '              path/to/file2.bam\n'
    log.info '              ...'
    log.info '				Note. Each bam file needs an index.'
    log.info '    --genome      Absolute path for file with genome as required by BEDTools.'
    log.info '                  The genome file spec is defined by BEDTools and has the following format:'
    log.info '                   <seq_name><\t><chrom_size>'
    log.info '    --window      Number of bases for each of the BEDTools generated windows.'
    log.info '    --outprefix   Prefix for output file.'
    log.info ''
    exit 1
}

log.info 'Starting the analysis.....'

bamPaths = file(params.file)
genomeFile = file(params.genome)

process make_windows {
    /*
    Process to create genomic windows of a certain width (in bases)
    */
 
    memory '500 MB'

    input:
    file genome from genomeFile

    output:
    stdout ivals_ch
 
    """
    bedtools makewindows -g $genome -w ${params.window} |awk -F'\\t' '{print \$1":"\$2"-"\$3}' -
    """
}

ivals_ch.splitText().map{it -> it.trim()}.set{monoival_ch}

process get_cov {
	/*
	Process to run SAMTools depth on params.pos_file and get a
	pos file that will be used later
	*/
	tag "Depth for ival: $ival"
	
	memory '500 MB'

	input:
	val ival from monoival_ch
    file url from bamPaths

	output:
	file("*.cov.gz") into cov_chunks

	exec:
	def match = (ival =~ /(.*):(\d+)-(\d+)/)
	match.matches()
	def chrom = match.group(1)
	def start = match.group(2)
	def toadd=10-start.length()
    def nstart='0'*toadd+start
	
	script:
	"""
    samtools depth -a -d 0 -r ${ival} -f $url |bgzip -c > ${chrom}.${nstart}.cov.gz
    """
}

sorted_covchunks = cov_chunks.collect().sort { a,b ->
     return a.baseName <=> b.baseName
}

process merge_chunks {
	/*
	This process will collect and merge each of the genomic
	chunks generated above
	*/
	memory '500 MB'

	input:
		file(cov_f) from sorted_covchunks

	output:
		file("merged.cov.gz") into merged_file

	script:
	"""
	zcat $cov_f |gzip -c > merged.cov.gz
	"""
}

process aggregate_depth {
	/*
	This process will aggregate across all sample-level coverages
	to produce an aggregated number
	*/
	publishDir "results", mode: 'copy', overwrite: true

	memory '500 MB'

	input:
        file(merged_file) from merged_file

	output:
		file("out.cov.gz") into agg_file

	"""
	sum_covs.py --ifile ${merged_file} --prefix out
	"""
}